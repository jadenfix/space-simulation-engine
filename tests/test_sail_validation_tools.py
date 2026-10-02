"""Independent numerical fixtures and real-CLI integration; no synthetic flight evidence."""
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('sail_validation', ROOT/'scripts/sail_validation.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
ENGINE = Path(os.environ.get('SPACEWIND_TEST_BINARY', ROOT/'build/spacewind')).resolve()
TOOL = ROOT/'scripts/sail_validation.py'


def compute_fixture():
    return json.loads((ROOT/'configs/validation/compute_3d_example.json').read_text())


def convergence_fixture():
    # Manufactured observable f(h)=10+0.2*h^2, not a measurement of any sail.
    return {
        'schema': 'spacewind.convergence-input.v1', 'observable': 'manufactured_force_x',
        'units': 'N', 'solver_revision': 'manufactured-test-only',
        'geometry_sha256': hashlib.sha256(b'manufactured geometry').hexdigest(),
        'case_sha256': hashlib.sha256(b'manufactured case').hexdigest(),
        'levels': [{'h': h, 'value': 10+.2*h*h, 'standard_error': 0,
                    'run_sha256': hashlib.sha256(f'manufactured h={h}'.encode()).hexdigest()}
                   for h in [1, .5, .25]],
        'absolute_tolerance': 1e-5, 'relative_tolerance': .005, 'minimum_order': 1.8,
    }


class ComputeTests(unittest.TestCase):
    def test_reference_arithmetic(self):
        r = V.estimate_compute(compute_fixture())
        self.assertEqual(r['cell_count'], 256**3)
        self.assertEqual(r['macroparticle_count'], 1073741824)
        self.assertEqual(r['raw_particle_gib'], 64)
        self.assertEqual(r['planned_memory_gib'], 132)
        self.assertAlmostEqual(r['dt_s'], .9*(1000/256)/(V.C*math.sqrt(3)), places=18)
        self.assertEqual(r['steps'], 14769942)
        self.assertAlmostEqual(r['allocated_gpu_hours'], 44053.06795403946, places=5)
        self.assertAlmostEqual(r['debye_length_m'], 10.5131815908, places=8)
        self.assertEqual(r['direct_dft_forward_cell_mode_contributions'], (256**3)**2)
        self.assertFalse(r['hardware_validation_established'])

    def test_work_and_parallel_scaling(self):
        cfg = compute_fixture()
        del cfg['duration_s']
        cfg['steps'] = 1000
        a = V.estimate_compute(cfg)
        cfg['steps'] *= 2
        b = V.estimate_compute(cfg)
        self.assertEqual(b['particle_updates'], 2*a['particle_updates'])
        self.assertAlmostEqual(b['allocated_gpu_hours'], 2*a['allocated_gpu_hours'])
        cfg['gpus'] *= 2
        c = V.estimate_compute(cfg)
        self.assertEqual(c['allocated_gpu_hours'], b['allocated_gpu_hours'])
        self.assertEqual(c['elapsed_hours_at_declared_parallel_efficiency'], b['elapsed_hours_at_declared_parallel_efficiency']/2)
        cfg['parallel_efficiency'] = .5
        d = V.estimate_compute(cfg)
        self.assertEqual(d['allocated_gpu_hours'], 2*c['allocated_gpu_hours'])

    def test_total_ppc_and_2d(self):
        cfg = compute_fixture()
        cfg.update(cells=[128, 128], lengths_m=[100, 100], particles_per_cell_total=32)
        r = V.estimate_compute(cfg)
        self.assertEqual(r['macroparticle_count'], 128**2*32)
        self.assertAlmostEqual(r['timestep_ceilings_s']['em_cfl_s'], .9*(100/128)/(V.C*math.sqrt(2)))

    def test_electrostatic_still_resolves_fast_particles(self):
        cfg = compute_fixture()
        cfg['solver'] = 'electrostatic'
        r = V.estimate_compute(cfg)
        self.assertNotIn('em_cfl_s', r['timestep_ceilings_s'])
        self.assertEqual(r['dt_s'], .25*(1000/256)/1e8)
        self.assertIn('voltage-accelerated', ' '.join(r['warnings']))

    def test_sampling_ceiling_is_enforced(self):
        cfg = compute_fixture()
        cfg['dt_s'] = 1e-3
        with self.assertRaises(ValueError):
            V.estimate_compute(cfg)

    def test_memory_and_compute_budgets(self):
        cfg = compute_fixture()
        cfg.update(gpus=1, device_memory_gib=80, budget_gpu_hours=1)
        r = V.estimate_compute(cfg)
        self.assertFalse(r['aggregate_memory_fits'])
        self.assertFalse(r['within_declared_gpu_hour_budget'])
        self.assertEqual(r['memory_minimum_gpu_count'], 2)

    def test_input_boundaries(self):
        for key, value in [('cells', [0, 1, 1]), ('cells', [True, 4, 4]), ('cells', [4]),
                           ('lengths_m', [1, 2]), ('lengths_m', [0, 2, 3]),
                           ('duration_s', 0), ('duration_s', float('nan')),
                           ('gpus', 0), ('gpus', 1.5), ('particle_copies', 0),
                           ('parallel_efficiency', 1.1), ('field_bytes_per_cell', -1),
                           ('max_particle_speed_m_s', V.C), ('number_density_m3', float('inf')),
                           ('particles_per_cell_total', True), ('solver', 'fake'), ('cfl_safety', 1.1)]:
            with self.subTest(key=key, value=value):
                cfg = compute_fixture()
                cfg[key] = value
                with self.assertRaises(ValueError):
                    V.estimate_compute(cfg)

    def test_required_unknown_and_duration_exclusivity(self):
        for change in ('missing', 'unknown', 'both'):
            cfg = compute_fixture()
            if change == 'missing':
                del cfg['throughput']
            elif change == 'unknown':
                cfg['particles_per_species'] = 64
            else:
                cfg['steps'] = 10
            with self.assertRaises(ValueError):
                V.estimate_compute(cfg)

    def test_measured_rate_requires_declared_provenance(self):
        cfg = compute_fixture()
        cfg['throughput']['kind'] = 'user_supplied_measurement'
        with self.assertRaises(ValueError):
            V.estimate_compute(cfg)
        cfg['throughput'].update(benchmark_receipt_sha256='a'*64, solver_revision='example-only')
        r = V.estimate_compute(cfg)
        self.assertIn('not independently verified', r['warnings'][0])

    def test_strict_json_rejects_nonfinite_and_duplicates(self):
        for text in ['{"a":NaN}', '{"a":Infinity}', '{"a":1,"a":2}']:
            with self.assertRaises(ValueError):
                V.strict_json(text)


class ConvergenceTests(unittest.TestCase):
    def test_manufactured_second_order(self):
        r = V.check_convergence(convergence_fixture())
        self.assertTrue(r['numerical_screen_passed'])
        self.assertAlmostEqual(r['observed_order'], 2)
        self.assertAlmostEqual(r['richardson_extrapolated_value'], 10)
        self.assertAlmostEqual(r['fine_grid_gci_plus_3se'], .015625)
        self.assertFalse(r['physical_validation_established'])
        self.assertFalse(r['flight_validation_established'])

    def test_flat_data_does_not_manufacture_convergence(self):
        cfg = convergence_fixture()
        for level in cfg['levels']:
            level['value'] = 10
        r = V.check_convergence(cfg)
        self.assertFalse(r['numerical_screen_passed'])
        self.assertIsNone(r['observed_order'])

    def test_noise_dominated_is_blocked(self):
        cfg = convergence_fixture()
        for level in cfg['levels']:
            level['standard_error'] = .1
        self.assertFalse(V.check_convergence(cfg)['numerical_screen_passed'])

    def test_nonmonotone_and_divergent_are_blocked(self):
        for values in [(10, 12, 11), (10, 11, 13)]:
            cfg = convergence_fixture()
            for level, value in zip(cfg['levels'], values):
                level['value'] = value
            self.assertFalse(V.check_convergence(cfg)['numerical_screen_passed'])

    def test_tight_error_budget_is_blocked(self):
        cfg = convergence_fixture()
        cfg.update(absolute_tolerance=1e-10, relative_tolerance=1e-10)
        r = V.check_convergence(cfg)
        self.assertIn('estimated_error_exceeds_declared_tolerance', r['blockers'])

    def test_absolute_tolerance_for_small_observable(self):
        cfg = convergence_fixture()
        for level in cfg['levels']:
            level['value'] = .2*level['h']**2
        cfg.update(absolute_tolerance=.02, relative_tolerance=0)
        self.assertTrue(V.check_convergence(cfg)['numerical_screen_passed'])

    def test_bad_level_structure_and_provenance(self):
        variants = []
        a = convergence_fixture(); a['levels'][1]['h'] = .6; variants.append(a)
        a = convergence_fixture(); a['levels'][1]['h'] = 2; variants.append(a)
        a = convergence_fixture(); a['levels'][1]['standard_error'] = -1; variants.append(a)
        a = convergence_fixture(); a['levels'][1]['run_sha256'] = a['levels'][0]['run_sha256']; variants.append(a)
        a = convergence_fixture(); a['geometry_sha256'] = 'missing'; variants.append(a)
        a = convergence_fixture(); a['levels'].pop(); variants.append(a)
        a = convergence_fixture(); a['gci_safety_factor'] = 0; variants.append(a)
        for cfg in variants:
            with self.assertRaises(ValueError):
                V.check_convergence(cfg)


class CommandTests(unittest.TestCase):
    def invoke(self, *args):
        return subprocess.run([sys.executable, str(TOOL), *map(str, args)], cwd=ROOT,
                              capture_output=True, text=True, timeout=45)

    def test_compute_budget_exit_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.json'
            cfg = compute_fixture(); cfg['budget_gpu_hours'] = 1
            path.write_text(json.dumps(cfg))
            run = self.invoke('compute', path)
            self.assertEqual(run.returncode, 3, run.stderr)
            self.assertFalse(json.loads(run.stdout)['within_declared_gpu_hour_budget'])

    def test_report_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/'report.json'
            target.write_text('KEEP')
            run = self.invoke('compute', ROOT/'configs/validation/compute_3d_example.json', '--output', target)
            self.assertEqual(run.returncode, 2)
            self.assertEqual(target.read_text(), 'KEEP')

    def test_invalid_input_exits_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad.json'; path.write_text('{"schema":NaN}')
            run = self.invoke('compute', path)
            self.assertEqual(run.returncode, 2)
            self.assertNotIn('Traceback', run.stderr)

    def test_convergence_rejection_exit_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.json'
            cfg = convergence_fixture()
            for level in cfg['levels']: level['value'] = 0
            path.write_text(json.dumps(cfg))
            run = self.invoke('convergence', path)
            self.assertEqual(run.returncode, 3, run.stderr)
            self.assertFalse(json.loads(run.stdout)['physical_validation_established'])


class NativeIntegrationTests(unittest.TestCase):
    invoke = CommandTests.invoke
    @classmethod
    def setUpClass(cls):
        if not ENGINE.is_file():
            raise RuntimeError('Build the engine first: make validation-tools')

    def test_all_existing_configs_emit_unvalidated_receipts(self):
        with tempfile.TemporaryDirectory() as tmp:
            for source in sorted((ROOT/'configs').glob('*.cfg')):
                with self.subTest(config=source.name):
                    cfg = V.read_config(source)
                    cfg.update(duration_s='1', step_s='.1', max_step_s='.1', min_step_s='.1')
                    path = Path(tmp)/source.name
                    path.write_text(''.join(f'{k}={v}\n' for k, v in cfg.items()))
                    receipt = Path(tmp)/'receipt.json'
                    run = subprocess.run([str(ENGINE), 'simulate', str(path), str(Path(tmp)/'trace.csv'), str(receipt)],
                                         capture_output=True, text=True, timeout=30)
                    self.assertEqual(run.returncode, 0, run.stderr)
                    r = V.strict_json(receipt.read_text())
                    self.assertTrue(r['completed'])
                    v = r['field_sail_validation']
                    self.assertFalse(v['physical_validation_established'])
                    self.assertFalse(v['flight_validation_established'])
                    self.assertFalse(v['coupled_device_power_budget'])
                    if source.name == 'dynamic_soaring_shear.cfg':
                        self.assertEqual(v['electric']['status'], 'requires_additional_accounting')
                        self.assertAlmostEqual(v['electric']['outgoing_kinetic_power_lower_bound_ratio'], 1.8)

    def test_paired_sweep_and_negative_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/'sweep'
            before = V.sha256_file(ROOT/'configs/dynamic_soaring_shear.cfg')
            run = self.invoke('sweep', '--engine', ENGINE, '--output', target, '--duration', 1,
                              '--step', .1, '--cl', 0, 1.2)
            self.assertEqual(run.returncode, 0, run.stderr)
            report = V.strict_json((target/'summary.json').read_text())
            self.assertTrue(report['completed'])
            self.assertEqual(len(report['cases']), 6)
            self.assertEqual({r['synthetic_shear'] for r in report['cases']}, {False, True})
            self.assertEqual(sum(not r['electric_sail_enabled'] for r in report['cases']), 2)
            self.assertTrue(all(r['receipt_sha256'] for r in report['cases']))
            self.assertEqual(before, V.sha256_file(ROOT/'configs/dynamic_soaring_shear.cfg'))
            self.assertFalse(report['mission_feasibility_established'])
            run = self.invoke('sweep', '--engine', ENGINE, '--output', target)
            self.assertEqual(run.returncode, 2)

    def test_step_budget_prevents_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/'blocked'
            run = self.invoke('sweep', '--engine', ENGINE, '--output', target, '--max-total-steps', 1)
            self.assertEqual(run.returncode, 2)
            self.assertFalse(target.exists())

    def test_config_parser_rejects_ambiguous_overrides(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad.cfg'; path.write_text('mass_kg=1\nmass_kg=2\n')
            with self.assertRaises(ValueError): V.read_config(path)

    def test_failed_engine_is_not_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp)/'fake_engine'
            fake.write_text(f'#!{sys.executable}\nraise SystemExit(1)\n')
            fake.chmod(0o755)
            target = Path(tmp)/'failed'
            run = self.invoke('sweep', '--engine', fake, '--output', target, '--cl', 0)
            self.assertEqual(run.returncode, 1, run.stderr)
            r = json.loads((target/'summary.json').read_text())
            self.assertFalse(r['completed'])
            self.assertTrue(all(not x['completed'] for x in r['cases']))
            self.assertFalse(r['hardware_validation_established'])

    def test_malformed_receipt_fails_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp)/'fake_engine'
            fake.write_text(f'#!{sys.executable}\nimport pathlib,sys\npathlib.Path(sys.argv[4]).write_text("[]")\n')
            fake.chmod(0o755)
            target = Path(tmp)/'malformed'
            run = self.invoke('sweep', '--engine', fake, '--output', target, '--cl', 0)
            self.assertEqual(run.returncode, 1, run.stderr)
            self.assertNotIn('Traceback', run.stderr)
            r = json.loads((target/'summary.json').read_text())
            self.assertFalse(r['completed'])
            self.assertTrue(all('error' in x for x in r['cases']))

    def test_timeout_retains_failed_case_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp)/'fake_engine'
            fake.write_text(f'#!{sys.executable}\nimport time\ntime.sleep(10)\n')
            fake.chmod(0o755)
            target = Path(tmp)/'timeout'
            run = self.invoke('sweep', '--engine', fake, '--output', target, '--cl', 0,
                              '--timeout-seconds', .05)
            self.assertEqual(run.returncode, 1, run.stderr)
            r = json.loads((target/'summary.json').read_text())
            self.assertEqual(len(r['cases']), r['planned_case_count'])
            self.assertFalse(r['completed'])
            self.assertTrue(all('timed out' in x['error'] for x in r['cases']))
            self.assertFalse((target/'summary.json.tmp').exists())

    def test_identical_config_has_identical_native_trace(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = V.read_config(ROOT/'configs/dynamic_soaring_shear.cfg')
            cfg.update(duration_s='2', step_s='.1', max_step_s='.1', min_step_s='.1')
            path = Path(tmp)/'input.cfg'
            path.write_text(''.join(f'{k}={v}\n' for k,v in cfg.items()))
            outputs = []
            for i in range(2):
                trace = Path(tmp)/f'trace{i}.csv'
                receipt = Path(tmp)/f'receipt{i}.json'
                run = subprocess.run([str(ENGINE), 'simulate', str(path), str(trace), str(receipt)],
                                     capture_output=True, text=True, timeout=30)
                self.assertEqual(run.returncode, 0, run.stderr)
                outputs.append((trace.read_bytes(), receipt.read_bytes()))
            self.assertEqual(outputs[0], outputs[1])


if __name__ == '__main__':
    unittest.main()
