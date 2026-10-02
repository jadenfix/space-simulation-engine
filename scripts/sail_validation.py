#!/usr/bin/env python3
"""Bounded sail-validation tools. Standard library only; never submits GPU jobs.

The numerical screens deliberately cannot promote chamber or flight claims.
See docs/SAIL_VALIDATION_CAMPAIGN.md for assumptions and remaining evidence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import time

C = 299792458.0
QE = 1.602176634e-19
EPS0 = 8.8541878128e-12
ME = 9.1093837015e-31
MP = 1.67262192369e-27


def strict_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value):
        raise ValueError(f"nonfinite JSON number: {value}")

    return json.loads(text, object_pairs_hook=pairs, parse_constant=constant)


def shape(obj, required, optional=()):
    if not isinstance(obj, dict):
        raise ValueError("expected a JSON object")
    missing = set(required) - obj.keys()
    unknown = obj.keys() - set(required) - set(optional)
    if missing or unknown:
        raise ValueError(f"missing keys: {sorted(missing)}; unknown keys: {sorted(unknown)}")


def number(value, name, *, minimum=None, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def integer(value, name, *, minimum=1):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def digest(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def nonempty(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def sha256_file(path):
    hasher = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            hasher.update(block)
    return hasher.hexdigest()


def serialize(report):
    return json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n'


def estimate_compute(cfg):
    """Sizing arithmetic, not a solver benchmark or an accuracy guarantee."""
    shape(cfg, ('schema', 'solver', 'cells', 'lengths_m', 'particles_per_cell_total',
                'bytes_per_particle', 'particle_copies', 'field_bytes_per_cell',
                'number_density_m3', 'electron_temperature_ev', 'max_particle_speed_m_s',
                'throughput', 'gpus', 'parallel_efficiency'),
          ('steps', 'duration_s', 'dt_s', 'cfl_safety', 'plasma_frequency_dt_factor',
           'particle_crossing_fraction', 'device_memory_gib', 'budget_gpu_hours', 'extra_memory_bytes'))
    if cfg['schema'] != 'spacewind.compute-input.v1':
        raise ValueError('unsupported compute-input schema')
    if cfg['solver'] not in ('explicit_em', 'electrostatic'):
        raise ValueError('solver must be explicit_em or electrostatic')
    if not isinstance(cfg['cells'], list) or len(cfg['cells']) not in (2, 3):
        raise ValueError('cells must have two or three dimensions')
    if not isinstance(cfg['lengths_m'], list) or len(cfg['lengths_m']) != len(cfg['cells']):
        raise ValueError('lengths_m must match the cell dimensions')
    cells = [integer(v, 'cells') for v in cfg['cells']]
    lengths = [number(v, 'lengths_m', positive=True) for v in cfg['lengths_m']]
    dx = [length / count for length, count in zip(lengths, cells)]
    density = number(cfg['number_density_m3'], 'number_density_m3', positive=True)
    temperature = number(cfg['electron_temperature_ev'], 'electron_temperature_ev', positive=True)
    vmax = number(cfg['max_particle_speed_m_s'], 'max_particle_speed_m_s', positive=True)
    if vmax >= C:
        raise ValueError('max_particle_speed_m_s must be subluminal')
    factors = {key: number(cfg.get(key, default), key, positive=True) for key, default in
               [('cfl_safety', .9), ('plasma_frequency_dt_factor', .2), ('particle_crossing_fraction', .25)]}
    if any(v > 1 for v in factors.values()):
        raise ValueError('sampling/safety factors must be <= 1')
    omega_pe = math.sqrt(density * QE * QE / (EPS0 * ME))
    limits = {
        'electron_plasma_sampling_s': factors['plasma_frequency_dt_factor'] / omega_pe,
        'particle_crossing_sampling_s': factors['particle_crossing_fraction'] * min(dx) / vmax,
    }
    if cfg['solver'] == 'explicit_em':
        limits['em_cfl_s'] = factors['cfl_safety'] / (C * math.sqrt(sum(1 / d**2 for d in dx)))
    suggested_dt = min(limits.values())
    dt = number(cfg.get('dt_s', suggested_dt), 'dt_s', positive=True)
    if dt > suggested_dt * (1 + 1e-12):
        raise ValueError('dt_s exceeds a declared timestep ceiling')
    if ('steps' in cfg) == ('duration_s' in cfg):
        raise ValueError('provide exactly one of steps or duration_s')
    steps = (integer(cfg['steps'], 'steps') if 'steps' in cfg else
             math.ceil(number(cfg['duration_s'], 'duration_s', positive=True) / dt))
    ppc = integer(cfg['particles_per_cell_total'], 'particles_per_cell_total')
    bpp = integer(cfg['bytes_per_particle'], 'bytes_per_particle')
    copies = integer(cfg['particle_copies'], 'particle_copies')
    field_bytes = integer(cfg['field_bytes_per_cell'], 'field_bytes_per_cell', minimum=0)
    extra = integer(cfg.get('extra_memory_bytes', 0), 'extra_memory_bytes', minimum=0)
    gpus = integer(cfg['gpus'], 'gpus')
    efficiency = number(cfg['parallel_efficiency'], 'parallel_efficiency', positive=True)
    if efficiency > 1:
        raise ValueError('parallel_efficiency must be <= 1')
    throughput = cfg['throughput']
    shape(throughput, ('kind', 'particle_updates_per_gpu_second'), ('benchmark_receipt_sha256', 'solver_revision'))
    if throughput['kind'] not in ('assumed', 'user_supplied_measurement'):
        raise ValueError('throughput kind must be assumed or user_supplied_measurement')
    if throughput['kind'] == 'user_supplied_measurement':
        digest(throughput.get('benchmark_receipt_sha256'), 'benchmark_receipt_sha256')
        nonempty(throughput.get('solver_revision'), 'solver_revision')
    rate = number(throughput['particle_updates_per_gpu_second'], 'throughput', positive=True)
    ncell = math.prod(cells)
    nparticle = ncell * ppc
    raw_bytes = nparticle * bpp
    memory = raw_bytes * copies + ncell * field_bytes + extra
    gpu_hours_ideal = nparticle * steps / rate / 3600
    gpu_hours = gpu_hours_ideal / efficiency
    memory_fits = None
    minimum_gpus = None
    if 'device_memory_gib' in cfg:
        capacity = number(cfg['device_memory_gib'], 'device_memory_gib', positive=True) * 2**30
        minimum_gpus = math.ceil(memory / capacity)
        memory_fits = memory <= gpus * capacity
    within_budget = None
    if 'budget_gpu_hours' in cfg:
        within_budget = gpu_hours <= number(cfg['budget_gpu_hours'], 'budget_gpu_hours', positive=True)
    report = {
        'schema': 'spacewind.compute-report.v1', 'inputs': cfg,
        'cell_count': ncell, 'macroparticle_count': nparticle, 'particle_updates': nparticle * steps,
        'steps': steps, 'dt_s': dt, 'covered_physical_time_s': steps * dt,
        'timestep_ceilings_s': limits, 'limiting_timestep_ceiling': min(limits, key=limits.get),
        'cell_widths_m': dx, 'debye_length_m': math.sqrt(EPS0 * temperature / (density * QE)),
        'electron_inertial_length_m': C / omega_pe,
        'proton_inertial_length_m': C / math.sqrt(density * QE * QE / (EPS0 * MP)),
        'raw_particle_gib': raw_bytes / 2**30, 'planned_memory_gib': memory / 2**30,
        'ideal_gpu_hours': gpu_hours_ideal, 'allocated_gpu_hours': gpu_hours,
        'elapsed_hours_at_declared_parallel_efficiency': gpu_hours / gpus,
        'memory_minimum_gpu_count': minimum_gpus, 'aggregate_memory_fits': memory_fits,
        'within_declared_gpu_hour_budget': within_budget,
        'direct_dft_forward_cell_mode_contributions': ncell**2,
        'hardware_validation_established': False,
        'warnings': [
            'Throughput is an assumption or caller-supplied measurement, not independently verified here.',
            'Effective throughput must include deposition, field solves, synchronization and diagnostics.',
            'Memory includes only declared copies/arrays; allocator, load imbalance and output costs may add overhead.',
            'Aggregate memory fit does not establish feasible domain decomposition.',
            'Timestep ceilings are screening rules, not convergence or stability proofs for a complete device.',
            'max_particle_speed_m_s must cover voltage-accelerated electrons, not just the incident solar wind.',
            'This mesh does not establish wire/sheath resolution or a sufficient domain size/settling time.',
        ],
    }
    serialize(report)  # Reject overflow/nonfinite derived outputs before returning.
    return report


def check_convergence(cfg):
    """Three-level scalar Richardson/GCI screen; no physical claim promotion."""
    shape(cfg, ('schema', 'observable', 'units', 'solver_revision', 'geometry_sha256',
                'case_sha256', 'levels', 'absolute_tolerance', 'relative_tolerance', 'minimum_order'),
          ('maximum_order', 'gci_safety_factor'))
    if cfg['schema'] != 'spacewind.convergence-input.v1':
        raise ValueError('unsupported convergence schema')
    for key in ('observable', 'units', 'solver_revision'):
        nonempty(cfg[key], key)
    for key in ('geometry_sha256', 'case_sha256'):
        digest(cfg[key], key)
    levels = cfg['levels']
    if not isinstance(levels, list) or len(levels) != 3:
        raise ValueError('exactly three levels, coarse to fine, are required')
    h, f, se, run_hashes = [], [], [], []
    for level in levels:
        shape(level, ('h', 'value', 'standard_error', 'run_sha256'))
        h.append(number(level['h'], 'h', positive=True))
        f.append(number(level['value'], 'value'))
        se.append(number(level['standard_error'], 'standard_error', minimum=0))
        run_hashes.append(digest(level['run_sha256'], 'run_sha256'))
    if len(set(run_hashes)) != 3:
        raise ValueError('each level must identify a distinct run')
    if not h[0] > h[1] > h[2]:
        raise ValueError('resolution h must decrease from coarse to fine')
    ratio = h[0] / h[1]
    if not math.isclose(ratio, h[1] / h[2], rel_tol=1e-9):
        raise ValueError('this screen requires a constant refinement ratio')
    atol = number(cfg['absolute_tolerance'], 'absolute_tolerance', minimum=0)
    rtol = number(cfg['relative_tolerance'], 'relative_tolerance', minimum=0)
    if atol == rtol == 0:
        raise ValueError('at least one tolerance must be positive')
    pmin = number(cfg['minimum_order'], 'minimum_order', positive=True)
    pmax = number(cfg.get('maximum_order', 6), 'maximum_order', positive=True)
    safety = number(cfg.get('gci_safety_factor', 1.25), 'gci_safety_factor', minimum=1.25)
    if pmin > pmax:
        raise ValueError('minimum_order exceeds maximum_order')
    d1, d2 = f[1] - f[0], f[2] - f[1]
    blockers = []
    # Sum of standard errors is conservative without assuming cross-level independence.
    if abs(d1) <= 3 * (se[0] + se[1]) or abs(d2) <= 3 * (se[1] + se[2]):
        blockers.append('refinement_signal_unresolved_or_identical_outputs')
    if d1 * d2 <= 0:
        blockers.append('not_monotone_convergence')
    order = None
    error = None
    extrapolated = None
    tolerance = max(atol, rtol * abs(f[2]))
    if not blockers:
        order = math.log(abs(d1 / d2)) / math.log(ratio)
        if not pmin <= order <= pmax:
            blockers.append('observed_order_outside_declared_range')
        else:
            denominator = math.expm1(order * math.log(ratio))
            error = safety * abs(d2) / denominator + 3 * se[2]
            extrapolated = f[2] + d2 / denominator
            if error > tolerance:
                blockers.append('estimated_error_exceeds_declared_tolerance')
    report = {
        'schema': 'spacewind.convergence-report.v1', 'inputs': cfg,
        'numerical_screen_passed': not blockers, 'blockers': blockers,
        'observed_order': order, 'fine_grid_gci_plus_3se': error,
        'declared_tolerance': tolerance, 'richardson_extrapolated_value': extrapolated,
        'physical_validation_established': False, 'flight_validation_established': False,
        'nonclaim': 'One scalar asymptotic numerical screen only. Run IDs are caller supplied, not authenticated. '
                    'Requires separate particle/domain/boundary studies, correlated-sample uncertainty, '
                    'independent implementation and held-out measurements. Exact/flat solutions need an analytic oracle.',
    }
    serialize(report)
    return report


def read_config(path):
    cfg = {}
    for line_no, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.split('#', 1)[0].strip()
        if not line or line.startswith(';'):
            continue
        if '=' not in line:
            raise ValueError(f'line {line_no}: expected key=value')
        key, value = (part.strip() for part in line.split('=', 1))
        if not key or key in cfg:
            raise ValueError(f'line {line_no}: empty or duplicate key {key}')
        cfg[key] = value
    return cfg


def write_manifest(path, report):
    # All paths are inside a newly created run directory owned by this invocation.
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(serialize(report))
    temporary.replace(path)


def sweep(args):
    engine = Path(args.engine).resolve(strict=True)
    source = Path(args.config).resolve(strict=True)
    base = read_config(source)
    duration = number(args.duration, 'duration', positive=True)
    step = number(args.step, 'step', positive=True)
    timeout = number(args.timeout_seconds, 'timeout_seconds', positive=True)
    cds = [number(x, 'cd', minimum=0) for x in args.cd]
    cls = list(dict.fromkeys([0.0] + [number(x, 'cl') for x in args.cl]))
    cases = [(cd, cl, shear, True) for cd, cl, shear in itertools.product(cds, cls, (False, True))]
    cases += [(0.0, 0.0, shear, False) for shear in (False, True)]
    steps_per_case = math.ceil(duration / step)
    max_steps = integer(args.max_total_steps, 'max_total_steps')
    # One extra clipped final step can be needed after floating-point accumulation.
    planned_steps = (steps_per_case + 1) * len(cases)
    if len(cases) > 128 or planned_steps > max_steps:
        raise ValueError('sweep exceeds the case/step budget; reduce the design or explicitly raise max-total-steps')
    target = Path(args.output).resolve()
    target.mkdir(parents=True, exist_ok=False)
    (target/'source.cfg').write_bytes(source.read_bytes())
    report = {
        'schema': 'spacewind.sail-sweep.v1', 'base_config_sha256': sha256_file(source),
        'engine_binary_sha256': sha256_file(engine), 'planned_total_steps': planned_steps,
        'planned_case_count': len(cases), 'completed': False,
        'duration_s': duration, 'step_s': step, 'cases': [],
        'hardware_validation_established': False, 'mission_feasibility_established': False,
        'nonclaim': 'Sensitivity experiment on assumed electric-sail forces, not kinetic or chamber validation. '
                    'No-shear cases retain the analytic Parker environment/turbulence; they are NOT uniform-flow proofs. '
                    'Radiation, magnetic-sail and Lorentz paths are disabled; gravity is inherited or explicitly enabled. '
                    'No mission-success threshold or spacecraft power model is supplied.',
    }
    manifest = target/'summary.json'
    write_manifest(manifest, report)
    for index, (cd, cl, shear, enabled) in enumerate(cases):
        case = target/f'case_{index:03d}'
        case.mkdir()
        cfg = dict(base)
        cfg.update(integrator='rk4', duration_s=format(duration, '.17g'), step_s=format(step, '.17g'),
                   min_step_s=format(step, '.17g'), max_step_s=format(step, '.17g'),
                   output_every_steps=str(max(1, math.ceil(steps_per_case / 1000))),
                   electric_sail_cd=format(cd, '.17g'), electric_sail_cl=format(cl, '.17g'),
                   shear_enabled=str(shear).lower(), enable_electric_sail=str(enabled).lower(),
                   enable_radiation_pressure='false', enable_magnetic_sail='false', enable_lorentz_force='false')
        if args.enable_gravity:
            cfg.update(enable_gravity='true', enable_1pn='true')
        config_path = case/'input.cfg'
        config_path.write_text(''.join(f'{k} = {v}\n' for k, v in sorted(cfg.items())))
        row = {'id': case.name, 'cd': cd, 'cl': cl, 'synthetic_shear': shear,
               'electric_sail_enabled': enabled, 'config_sha256': sha256_file(config_path),
               'completed': False}
        start = time.perf_counter()
        try:
            run = subprocess.run([str(engine), 'simulate', str(config_path), str(case/'trajectory.csv'),
                                  str(case/'receipt.json')], capture_output=True, text=True, timeout=timeout)
            (case/'stdout.log').write_text(run.stdout)
            (case/'stderr.log').write_text(run.stderr)
            row['returncode'] = run.returncode
            if run.returncode != 0:
                raise ValueError('engine returned nonzero; see stderr.log')
            receipt = strict_json((case/'receipt.json').read_text())
            if not isinstance(receipt, dict) or receipt.get('completed') is not True:
                raise ValueError('missing completed receipt')
            validation = receipt.get('field_sail_validation')
            if not isinstance(validation, dict) or validation.get('schema') != 'spacewind.field-sail-screen.v1':
                raise ValueError('missing supported field-sail validation section')
            if validation.get('physical_validation_established') is not False or validation.get('flight_validation_established') is not False:
                raise ValueError('unexpected physical/flight claim in phenomenological sweep')
            accepted_steps = integer(receipt.get('accepted_steps'), 'accepted_steps')
            if accepted_steps > steps_per_case + 1:
                raise ValueError('native engine exceeded the planned per-case step budget')
            with (case/'trajectory.csv').open() as stream:
                reader = csv.DictReader(stream)
                first = next(reader)
                last = first
                for record in reader:
                    last = record
            initial = [number(float(first[k]), k) for k in ('vx_m_s', 'vy_m_s', 'vz_m_s')]
            final = [number(float(last[k]), k) for k in ('vx_m_s', 'vy_m_s', 'vz_m_s')]
            row.update(completed=True, accepted_steps=receipt['accepted_steps'],
                       delta_velocity_norm_m_s=math.sqrt(sum((b-a)**2 for a,b in zip(initial, final))),
                       final_speed_m_s=number(receipt['final_speed_m_s'], 'final_speed_m_s'),
                       field_sail_validation=receipt['field_sail_validation'],
                       receipt_sha256=sha256_file(case/'receipt.json'),
                       trajectory_sha256=sha256_file(case/'trajectory.csv'))
        except (OSError, ValueError, KeyError, StopIteration, subprocess.TimeoutExpired) as error:
            row['error'] = str(error)
        row['elapsed_wall_s'] = time.perf_counter() - start
        report['cases'].append(row)
        write_manifest(manifest, report)
    report['completed'] = all(row['completed'] for row in report['cases'])
    write_manifest(manifest, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ('compute', 'convergence'):
        item = sub.add_parser(command)
        item.add_argument('input', type=Path)
        item.add_argument('--output', type=Path, help='new report file; existing files are not overwritten')
    item = sub.add_parser('sweep')
    item.add_argument('--engine', default='build/spacewind')
    item.add_argument('--config', default='configs/dynamic_soaring_shear.cfg')
    item.add_argument('--output', required=True, help='new directory; existing directories are not overwritten')
    item.add_argument('--duration', type=float, default=600)
    item.add_argument('--step', type=float, default=5)
    item.add_argument('--cd', type=float, nargs='+', default=[.4])
    item.add_argument('--cl', type=float, nargs='+', default=[0, .3, 1.2])
    item.add_argument('--timeout-seconds', type=float, default=30)
    item.add_argument('--max-total-steps', type=int, default=1_000_000)
    item.add_argument('--enable-gravity', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.command == 'sweep':
            report = sweep(args)
            print(serialize({'completed': report['completed'], 'cases': len(report['cases']),
                             'summary': str(Path(args.output)/'summary.json')}), end='')
            return 0 if report['completed'] else 1
        cfg = strict_json(args.input.read_text())
        report = estimate_compute(cfg) if args.command == 'compute' else check_convergence(cfg)
        text = serialize(report)
        if args.output:
            with args.output.open('x') as stream:
                stream.write(text)
        else:
            print(text, end='')
        if args.command == 'convergence':
            return 0 if report['numerical_screen_passed'] else 3
        return 3 if report['within_declared_gpu_hour_budget'] is False or report['aggregate_memory_fits'] is False else 0
    except (OSError, ValueError, TypeError, OverflowError, ZeroDivisionError) as error:
        print(f'validation error: {error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
