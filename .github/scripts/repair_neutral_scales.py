from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    target = Path(path)
    text = target.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise SystemExit(
            f"{path}: expected one match, found {count}: {old!r}"
        )
    target.write_text(text.replace(old, new, 1), encoding="utf-8")


support = Path("assurance/src/em_pic3d_support.inc")
support_text = support.read_text(encoding="utf-8")
marker = "static double pic3d_absolute_particle_charge_C("
if marker not in support_text:
    support_text += r'''

static double pic3d_absolute_particle_charge_C(
    const swa_em_pic3d_particle *particles,
    size_t particle_count
) {
    swa_kahan total = {0.0, 0.0};
    size_t p;
    if (particles == NULL && particle_count != 0U) {
        return NAN;
    }
    for (p = 0U; p < particle_count; ++p) {
        swa_kahan_add(
            &total,
            fabs(particles[p].charge_C * particles[p].macro_weight)
        );
    }
    return total.sum;
}

static double pic3d_particle_momentum_scale_Ns(
    const swa_em_pic3d_particle *particles,
    size_t particle_count
) {
    swa_kahan total = {0.0, 0.0};
    size_t p;
    if (particles == NULL && particle_count != 0U) {
        return NAN;
    }
    for (p = 0U; p < particle_count; ++p) {
        const double speed2 = swa_vdot(
            particles[p].velocity_mps,
            particles[p].velocity_mps
        );
        const double gamma = 1.0 / sqrt(
            1.0 - speed2 / (SWA_C * SWA_C)
        );
        swa_kahan_add(
            &total,
            particles[p].macro_weight * particles[p].mass_kg *
            gamma * swa_vnorm(particles[p].velocity_mps)
        );
    }
    return total.sum;
}
'''
    support.write_text(support_text, encoding="utf-8")

replace_once(
    "assurance/src/em_pic3d_diagnostics.inc",
    r'''        const double gamma = 1.0 / sqrt(
            1.0 - speed2 / (SWA_C * SWA_C)
        );
        swa_kahan_add(
            &energy,
            particles[p].macro_weight * particles[p].mass_kg *
            SWA_C * SWA_C * (gamma - 1.0)
        );''',
    r'''        const double beta2 = speed2 / (SWA_C * SWA_C);
        const double inverse_gamma = sqrt(1.0 - beta2);
        const double gamma_minus_one =
            beta2 / (inverse_gamma * (1.0 + inverse_gamma));
        swa_kahan_add(
            &energy,
            particles[p].macro_weight * particles[p].mass_kg *
            SWA_C * SWA_C * gamma_minus_one
        );''',
)

replace_once(
    "assurance/src/em_pic3d_spectral.inc",
    r'''    swa_kahan squared = {0.0, 0.0};
    swa_kahan mean = {0.0, 0.0};
    double scale = DBL_MIN;
    size_t i;
    size_t j;
    size_t k;
    *maximum = 0.0;
''',
    r'''    swa_kahan squared = {0.0, 0.0};
    swa_kahan mean = {0.0, 0.0};
    const double inverse_length = sqrt(
        1.0 / (dx * dx) +
        1.0 / (dy * dy) +
        1.0 / (dz * dz)
    );
    const double total_volume =
        state->length_x_m * state->length_y_m * state->length_z_m;
    const double absolute_charge_rate =
        pic3d_absolute_particle_charge_C(
            state->particles, state->particle_count
        ) / (total_volume * state->dt_s);
    double maximum_current = 0.0;
    double scale;
    size_t index;
    size_t i;
    size_t j;
    size_t k;
    *maximum = 0.0;
    for (index = 0U; index < state->cell_count; ++index) {
        maximum_current = pic3d_max4(
            maximum_current,
            fabs(jx[index]), fabs(jy[index]), fabs(jz[index])
        );
    }
    scale = pic3d_max3(
        absolute_charge_rate,
        maximum_current * inverse_length,
        DBL_MIN
    );
''',
)

replace_once(
    "assurance/src/em_pic3d_step.inc",
    r'''    charge_scale = pic3d_max4(
        fabs(out->initial_charge_C),
        fabs(out->final_charge_C),
        charge_l1_initial,
        pic3d_max2(charge_l1_final, DBL_MIN)
    );''',
    r'''    charge_scale = pic3d_max4(
        fabs(out->initial_charge_C),
        fabs(out->final_charge_C),
        pic3d_absolute_particle_charge_C(
            particle_initial, state->particle_count
        ),
        pic3d_max2(
            pic3d_max2(charge_l1_initial, charge_l1_final),
            DBL_MIN
        )
    );''',
)

replace_once(
    "assurance/src/em_pic3d_step.inc",
    r'''    momentum_scale = pic3d_max4(
        swa_vnorm(swa_vadd(
            out->initial_particle_momentum_Ns,
            out->initial_field_momentum_Ns
        )),
        swa_vnorm(swa_vadd(
            out->final_particle_momentum_Ns,
            out->final_field_momentum_Ns
        )),
        swa_vnorm(swa_vsub(
            out->final_particle_momentum_Ns,
            out->initial_particle_momentum_Ns
        )),
        pic3d_max2(
            swa_vnorm(swa_vsub(
                out->final_field_momentum_Ns,
                out->initial_field_momentum_Ns
            )),
            DBL_MIN
        )
    );''',
    r'''    {
        const double particle_momentum_scale = pic3d_max2(
            pic3d_particle_momentum_scale_Ns(
                particle_initial, state->particle_count
            ),
            pic3d_particle_momentum_scale_Ns(
                particle_final, state->particle_count
            )
        );
        const double field_momentum_scale = pic3d_max2(
            out->initial_field_energy_J / SWA_C,
            out->final_field_energy_J / SWA_C
        );
        const double exchange_momentum_scale = pic3d_max4(
            particle_momentum_scale,
            field_momentum_scale,
            swa_vnorm(swa_vsub(
                out->final_particle_momentum_Ns,
                out->initial_particle_momentum_Ns
            )),
            swa_vnorm(swa_vsub(
                out->final_field_momentum_Ns,
                out->initial_field_momentum_Ns
            ))
        );
        momentum_scale = pic3d_max4(
            swa_vnorm(swa_vadd(
                out->initial_particle_momentum_Ns,
                out->initial_field_momentum_Ns
            )),
            swa_vnorm(swa_vadd(
                out->final_particle_momentum_Ns,
                out->final_field_momentum_Ns
            )),
            exchange_momentum_scale,
            DBL_MIN
        );
    }''',
)

replace_once(
    "assurance/tests/test_em_pic3d_plasma_diagnostic.c",
    "    const double target_omega = 2.0 * SWA_PI / (64.0 * dt);",
    "    const double target_omega = 2.0 * SWA_PI / (1024.0 * dt);",
)
