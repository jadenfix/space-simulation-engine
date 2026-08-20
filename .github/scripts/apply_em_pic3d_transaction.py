from __future__ import annotations

import json
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    target = Path(path)
    text = target.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise SystemExit(
            f"{path}: expected exactly one replacement anchor, found {count}"
        )
    target.write_text(text.replace(old, new, 1), encoding="utf-8")


replace_once(
    "assurance/include/spacewind/em_pic3d.h",
    """    int field_work_passes;
    int momentum_passes;
    int passes;
} swa_em_pic3d_result;""",
    """    int field_work_passes;
    int momentum_passes;
    int state_committed;
    int passes;
} swa_em_pic3d_result;""",
)

replace_once(
    "assurance/src/em_pic3d_step.inc",
    """    if (!out->finite) {
        goto cleanup;
    }
    if (state->particle_count > 0U) {
        memcpy(
            state->particles, particle_final,
            state->particle_count * sizeof(*particle_final)
        );
    }
    memcpy(state->ex_V_m, ex_final, count * sizeof(double));
    memcpy(state->ey_V_m, ey_final, count * sizeof(double));
    memcpy(state->ez_V_m, ez_final, count * sizeof(double));
    memcpy(state->bx_T, bx_final, count * sizeof(double));
    memcpy(state->by_T, by_final, count * sizeof(double));
    memcpy(state->bz_T, bz_final, count * sizeof(double));
    ok = 1;""",
    """    if (!out->finite) {
        goto cleanup;
    }
    out->state_committed = 0;
    if (out->passes) {
        if (state->particle_count > 0U) {
            memcpy(
                state->particles, particle_final,
                state->particle_count * sizeof(*particle_final)
            );
        }
        memcpy(state->ex_V_m, ex_final, count * sizeof(double));
        memcpy(state->ey_V_m, ey_final, count * sizeof(double));
        memcpy(state->ez_V_m, ez_final, count * sizeof(double));
        memcpy(state->bx_T, bx_final, count * sizeof(double));
        memcpy(state->by_T, by_final, count * sizeof(double));
        memcpy(state->bz_T, bz_final, count * sizeof(double));
        out->state_committed = 1;
    }
    ok = 1;""",
)

replace_once(
    "assurance/src/em_pic3d_step.inc",
    """        \"  \\\"local_current_passes\\\": %s,\\n\"
        \"  \\\"spectral_oracle_passes\\\": %s,\\n\"
        \"  \\\"passes\\\": %s,\\n\"""",
    """        \"  \\\"local_current_passes\\\": %s,\\n\"
        \"  \\\"spectral_oracle_passes\\\": %s,\\n\"
        \"  \\\"state_committed\\\": %s,\\n\"
        \"  \\\"passes\\\": %s,\\n\"""",
)

replace_once(
    "assurance/src/em_pic3d_step.inc",
    """        result->local_current_passes ? \"true\" : \"false\",
        result->spectral_oracle_passes ? \"true\" : \"false\",
        result->passes ? \"true\" : \"false\"""",
    """        result->local_current_passes ? \"true\" : \"false\",
        result->spectral_oracle_passes ? \"true\" : \"false\",
        result->state_committed ? \"true\" : \"false\",
        result->passes ? \"true\" : \"false\"""",
)

test_insert = r'''
static void capture_field_state(
    const swa_em_pic3d *state,
    double *snapshot
) {
    const size_t bytes = state->cell_count * sizeof(double);
    memcpy(snapshot, state->ex_V_m, bytes);
    memcpy(snapshot + state->cell_count, state->ey_V_m, bytes);
    memcpy(snapshot + 2U * state->cell_count, state->ez_V_m, bytes);
    memcpy(snapshot + 3U * state->cell_count, state->bx_T, bytes);
    memcpy(snapshot + 4U * state->cell_count, state->by_T, bytes);
    memcpy(snapshot + 5U * state->cell_count, state->bz_T, bytes);
}

static int field_state_matches_snapshot(
    const swa_em_pic3d *state,
    const double *snapshot
) {
    const size_t bytes = state->cell_count * sizeof(double);
    return memcmp(state->ex_V_m, snapshot, bytes) == 0 &&
           memcmp(state->ey_V_m, snapshot + state->cell_count, bytes) == 0 &&
           memcmp(state->ez_V_m, snapshot + 2U * state->cell_count, bytes) == 0 &&
           memcmp(state->bx_T, snapshot + 3U * state->cell_count, bytes) == 0 &&
           memcmp(state->by_T, snapshot + 4U * state->cell_count, bytes) == 0 &&
           memcmp(state->bz_T, snapshot + 5U * state->cell_count, bytes) == 0;
}

static void test_transactional_step_acceptance(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    double *snapshot = NULL;
    const size_t rejected_cells = 4U;
    const double rejected_dx = 1.0 / (double)rejected_cells;
    const double rejected_dt = 1.05 * rejected_dx /
        (SWA_C * sqrt(3.0));
    size_t index;

    memset(&state, 0, sizeof(state));
    check_true(
        swa_em_pic3d_init(
            &state,
            rejected_cells, rejected_cells, rejected_cells, 0U,
            1.0, 1.0, 1.0, rejected_dt
        ),
        "transactional rejected-step grid initializes"
    );
    fill_divergence_free_fields(&state);
    snapshot = calloc(6U * state.cell_count, sizeof(double));
    check_true(snapshot != NULL,
               "allocate rejected-step field snapshot");
    if (snapshot != NULL) {
        capture_field_state(&state, snapshot);
        limits.maximum_courant = 0.8;
        check_true(
            swa_em_pic3d_step(&state, &limits, &result),
            "rejected transactional step returns an audit receipt"
        );
        check_true(!result.passes && !result.state_committed,
                   "failed scientific gates do not commit candidate state");
        check_true(field_state_matches_snapshot(&state, snapshot),
                   "rejected 3D3V step preserves every field bit");
    }
    free(snapshot);
    swa_em_pic3d_destroy(&state);

    memset(&state, 0, sizeof(state));
    limits = permissive_limits();
    check_true(
        swa_em_pic3d_init(
            &state, 17U, 1U, 1U, 0U,
            1.0, 1000.0, 1000.0,
            0.20 * (1.0 / 17.0) / SWA_C
        ),
        "transactional accepted-step grid initializes"
    );
    initialize_plane_wave(&state, 0, 2U, 1.0, 1);
    snapshot = calloc(6U * state.cell_count, sizeof(double));
    check_true(snapshot != NULL,
               "allocate accepted-step field snapshot");
    if (snapshot != NULL) {
        capture_field_state(&state, snapshot);
        check_true(
            swa_em_pic3d_step(&state, &limits, &result),
            "accepted transactional step executes"
        );
        check_true(result.passes && result.state_committed,
                   "passing 3D3V step commits exactly one candidate state");
        check_true(!field_state_matches_snapshot(&state, snapshot),
                   "accepted 3D3V step advances the field state");
        for (index = 0U; index < state.cell_count; ++index) {
            check_true(isfinite(state.ex_V_m[index]) &&
                       isfinite(state.ey_V_m[index]) &&
                       isfinite(state.ez_V_m[index]) &&
                       isfinite(state.bx_T[index]) &&
                       isfinite(state.by_T[index]) &&
                       isfinite(state.bz_T[index]),
                       "committed transactional state remains finite");
        }
    }
    free(snapshot);
    swa_em_pic3d_destroy(&state);
}

'''
replace_once(
    "assurance/tests/test_em_pic3d.c",
    "\nint main(void) {",
    "\n" + test_insert + "int main(void) {",
)
replace_once(
    "assurance/tests/test_em_pic3d.c",
    """    test_coupled_3d3v_reference();
    test_fail_closed_gates();""",
    """    test_coupled_3d3v_reference();
    test_fail_closed_gates();
    test_transactional_step_acceptance();""",
)

doc_anchor = "## Fail-closed gates\n"
doc_insert = """## Transactional step acceptance

The timestep is evaluated into candidate particle and field buffers. The public
state is committed only when every declared acceptance gate passes. A finite
rejected step still returns its complete diagnostic receipt, but leaves all
particle and field bytes unchanged. Internal errors likewise leave the state
unchanged. The receipt exposes `state_committed` so callers cannot confuse a
computed candidate with an accepted physical trajectory point.

This prevents a failed CFL, continuity, Gauss, magnetic-divergence, current,
nonlinear, energy, momentum, or subluminal gate from contaminating subsequent
steps.

"""
replace_once(
    "docs/EM_PIC3D_REFERENCE.md",
    doc_anchor,
    doc_insert + doc_anchor,
)

claims_path = Path("claims/em-pic3d-claims.json")
claims = json.loads(claims_path.read_text(encoding="utf-8"))
claims["version"] = "0.2.0"
target = next(item for item in claims["claims"] if item["id"] == "SW-EM3D-0004")
suffix = (
    " Rejected finite timesteps are additionally required to preserve "
    "the complete public particle and field state bit-for-bit; only a "
    "fully passing receipt may commit the candidate state."
)
if suffix.strip() not in target["statement"]:
    target["statement"] += suffix
claims_path.write_text(json.dumps(claims, indent=2) + "\n", encoding="utf-8")
