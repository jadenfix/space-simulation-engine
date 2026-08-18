#include "spacewind/field_ledger.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static size_t cell_index(size_t i, size_t j, size_t nx) {
    return j * nx + i;
}
static size_t xface_index(size_t i, size_t j, size_t nx) {
    return j * (nx + 1U) + i;
}
static size_t yface_index(size_t i, size_t j, size_t nx) {
    return j * nx + i;
}
static size_t vertex_index(size_t i, size_t j, size_t nx) {
    return j * (nx + 1U) + i;
}
static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}
static int near(double a, double b, double relative, double absolute) {
    return fabs(a - b) <=
           fmax(absolute, relative * fmax(fabs(a), fabs(b)));
}
static int near_vec(swa_vec3 a, swa_vec3 b,
                    double relative, double absolute) {
    return near(a.x, b.x, relative, absolute) &&
           near(a.y, b.y, relative, absolute) &&
           near(a.z, b.z, relative, absolute);
}

static void test_charge_continuity(
    swa_charge_continuity_result *receipt_result
) {
    enum {
        NX = 7,
        NY = 5,
        CELLS = NX * NY,
        XFACES = (NX + 1) * NY,
        YFACES = NX * (NY + 1)
    };
    const double dx = 0.4;
    const double dy = 0.3;
    const double dt = 0.02;
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double jx[XFACES];
    double jy[YFACES];
    swa_charge_continuity_input input;
    swa_charge_continuity_result result;
    swa_charge_continuity_result corrupted;
    size_t i;
    size_t j;

    for (j = 0U; j < NY; ++j) {
        for (i = 0U; i <= NX; ++i) {
            jx[xface_index(i, j, NX)] =
                1.0e-3 * (0.7 * (double)i -
                          0.2 * (double)j +
                          0.03 * (double)(i * j));
        }
    }
    for (j = 0U; j <= NY; ++j) {
        for (i = 0U; i < NX; ++i) {
            jy[yface_index(i, j, NX)] =
                1.0e-3 * (-0.4 * (double)j +
                          0.1 * (double)i -
                          0.02 * (double)(i * j));
        }
    }
    for (j = 0U; j < NY; ++j) {
        for (i = 0U; i < NX; ++i) {
            const size_t c = cell_index(i, j, NX);
            const double divergence =
                (jx[xface_index(i + 1U, j, NX)] -
                 jx[xface_index(i, j, NX)]) / dx +
                (jy[yface_index(i, j + 1U, NX)] -
                 jy[yface_index(i, j, NX)]) / dy;
            rho_initial[c] =
                2.0e-5 + 1.0e-7 * (double)i -
                2.0e-7 * (double)j;
            rho_final[c] = rho_initial[c] - dt * divergence;
        }
    }

    input.nx = NX;
    input.ny = NY;
    input.dx_m = dx;
    input.dy_m = dy;
    input.dt_s = dt;
    input.rho_initial_C_m3 = rho_initial;
    input.rho_final_C_m3 = rho_final;
    input.jx_face_A_m2 = jx;
    input.jy_face_A_m2 = jy;

    check_true(
        swa_audit_charge_continuity(
            &input, 1.0e-12, 1.0e-12, 1.0e-12, &result
        ),
        "continuity audit executes"
    );
    check_true(result.finite, "continuity outputs finite");
    check_true(result.local_passes, "manufactured local continuity passes");
    check_true(result.global_passes, "manufactured global continuity passes");
    check_true(result.passes, "manufactured continuity passes");
    check_true(
        fabs(result.divergence_theorem_error_C_m) < 1.0e-15,
        "discrete divergence theorem closes"
    );
    *receipt_result = result;

    rho_final[cell_index(3U, 2U, NX)] += dt * 1.0e-4;
    check_true(
        swa_audit_charge_continuity(
            &input, 1.0e-12, 1.0e-12, 1.0e-12, &corrupted
        ),
        "corrupted continuity audit executes"
    );
    check_true(!corrupted.local_passes,
               "local charge corruption detected");
    check_true(!corrupted.global_passes,
               "global charge corruption detected");
    check_true(!corrupted.passes,
               "corrupted continuity rejected");
    rho_final[cell_index(3U, 2U, NX)] -= dt * 1.0e-4;

    input.dt_s = 0.0;
    check_true(
        !swa_audit_charge_continuity(
            &input, 1.0e-12, 1.0e-12, 1.0e-12, &corrupted
        ),
        "zero continuity timestep rejected"
    );
}

static void test_gauss_law(
    swa_field_constraint_result *receipt_result
) {
    enum {
        NX = 8,
        NY = 6,
        CELLS = NX * NY,
        XFACES = (NX + 1) * NY,
        YFACES = NX * (NY + 1)
    };
    const double dx = 0.25;
    const double dy = 0.35;
    double ex[XFACES];
    double ey[YFACES];
    double rho[CELLS];
    swa_staggered_scalar_constraint_input input;
    swa_field_constraint_result result;
    swa_field_constraint_result corrupted;
    size_t i;
    size_t j;

    for (j = 0U; j < NY; ++j) {
        for (i = 0U; i <= NX; ++i) {
            ex[xface_index(i, j, NX)] =
                1.0e3 * (0.2 * (double)i +
                         0.03 * (double)j);
        }
    }
    for (j = 0U; j <= NY; ++j) {
        for (i = 0U; i < NX; ++i) {
            ey[yface_index(i, j, NX)] =
                1.0e3 * (-0.1 * (double)j +
                         0.02 * (double)i);
        }
    }
    for (j = 0U; j < NY; ++j) {
        for (i = 0U; i < NX; ++i) {
            const double divergence =
                (ex[xface_index(i + 1U, j, NX)] -
                 ex[xface_index(i, j, NX)]) / dx +
                (ey[yface_index(i, j + 1U, NX)] -
                 ey[yface_index(i, j, NX)]) / dy;
            rho[cell_index(i, j, NX)] = SWA_EPS0 * divergence;
        }
    }
    input.nx = NX;
    input.ny = NY;
    input.dx_m = dx;
    input.dy_m = dy;
    input.x_face_field = ex;
    input.y_face_field = ey;
    input.cell_source = rho;
    check_true(
        swa_audit_gauss_law(
            &input, 1.0e-12, 1.0e-18, &result
        ),
        "Gauss audit executes"
    );
    check_true(result.passes, "manufactured Gauss law passes");
    *receipt_result = result;

    rho[cell_index(4U, 3U, NX)] += 1.0e-8;
    check_true(
        swa_audit_gauss_law(
            &input, 1.0e-12, 1.0e-18, &corrupted
        ),
        "corrupted Gauss audit executes"
    );
    check_true(!corrupted.passes, "Gauss-law corruption detected");
}

static void test_magnetic_divergence(
    swa_field_constraint_result *receipt_result
) {
    enum {
        NX = 9,
        NY = 7,
        VERTICES = (NX + 1) * (NY + 1),
        XFACES = (NX + 1) * NY,
        YFACES = NX * (NY + 1)
    };
    const double dx = 0.2;
    const double dy = 0.15;
    double psi[VERTICES];
    double bx[XFACES];
    double by[YFACES];
    swa_staggered_scalar_constraint_input input;
    swa_field_constraint_result result;
    swa_field_constraint_result corrupted;
    size_t i;
    size_t j;

    for (j = 0U; j <= NY; ++j) {
        for (i = 0U; i <= NX; ++i) {
            psi[vertex_index(i, j, NX)] =
                1.0e-8 * (
                    sin(0.3 * (double)i) +
                    cos(0.2 * (double)j) +
                    0.1 * (double)(i * j)
                );
        }
    }
    for (j = 0U; j < NY; ++j) {
        for (i = 0U; i <= NX; ++i) {
            bx[xface_index(i, j, NX)] =
                (psi[vertex_index(i, j + 1U, NX)] -
                 psi[vertex_index(i, j, NX)]) / dy;
        }
    }
    for (j = 0U; j <= NY; ++j) {
        for (i = 0U; i < NX; ++i) {
            by[yface_index(i, j, NX)] =
                -(psi[vertex_index(i + 1U, j, NX)] -
                  psi[vertex_index(i, j, NX)]) / dx;
        }
    }
    input.nx = NX;
    input.ny = NY;
    input.dx_m = dx;
    input.dy_m = dy;
    input.x_face_field = bx;
    input.y_face_field = by;
    input.cell_source = NULL;
    check_true(
        swa_audit_magnetic_divergence(
            &input, 1.0e-11, 1.0e-18, &result
        ),
        "magnetic divergence audit executes"
    );
    check_true(result.passes,
               "discrete curl construction is divergence free");
    *receipt_result = result;

    bx[xface_index(5U, 3U, NX)] += 1.0e-7;
    check_true(
        swa_audit_magnetic_divergence(
            &input, 1.0e-11, 1.0e-18, &corrupted
        ),
        "corrupted magnetic divergence audit executes"
    );
    check_true(!corrupted.passes,
               "magnetic monopole corruption detected");
}

static void test_integral_ledgers(
    swa_poynting_ledger_result *energy_receipt,
    swa_field_momentum_result *momentum_receipt
) {
    swa_poynting_ledger energy;
    swa_poynting_ledger_result energy_result;
    swa_poynting_ledger_result corrupted_energy;
    swa_field_momentum_ledger momentum;
    swa_field_momentum_result momentum_result;
    swa_field_momentum_result corrupted_momentum;

    energy.field_energy_initial_J = 100.0;
    energy.field_energy_final_J = 80.0;
    energy.particle_work_J = 12.0;
    energy.outward_poynting_energy_J = 10.0;
    energy.impressed_source_energy_J = 2.0;
    check_true(
        swa_audit_poynting_ledger(
            &energy, 1.0e-14, 1.0e-13, &energy_result
        ),
        "Poynting ledger audit executes"
    );
    check_true(energy_result.passes, "Poynting ledger closes");
    *energy_receipt = energy_result;
    energy.field_energy_final_J += 0.1;
    check_true(
        swa_audit_poynting_ledger(
            &energy, 1.0e-14, 1.0e-13, &corrupted_energy
        ),
        "corrupted Poynting audit executes"
    );
    check_true(!corrupted_energy.passes,
               "untracked field energy detected");

    momentum.field_momentum_initial_Ns = swa_v3(1.0, 2.0, 3.0);
    momentum.mechanical_impulse_Ns = swa_v3(0.5, -0.2, 0.1);
    momentum.outward_maxwell_impulse_Ns = swa_v3(0.1, 0.3, -0.2);
    momentum.impressed_external_impulse_Ns = swa_v3(0.2, 0.4, 0.5);
    momentum.field_momentum_final_Ns = swa_v3(0.6, 2.3, 3.6);
    check_true(
        swa_audit_field_momentum(
            &momentum, 1.0e-14, 1.0e-13, &momentum_result
        ),
        "field momentum audit executes"
    );
    check_true(momentum_result.passes,
               "field and mechanical momentum close");
    *momentum_receipt = momentum_result;
    momentum.field_momentum_final_Ns.z += 1.0e-3;
    check_true(
        swa_audit_field_momentum(
            &momentum, 1.0e-14, 1.0e-13, &corrupted_momentum
        ),
        "corrupted momentum audit executes"
    );
    check_true(!corrupted_momentum.passes,
               "untracked electromagnetic impulse detected");
}

static void test_local_electromagnetism(void) {
    size_t i;
    swa_maxwell_traction_result traction;
    swa_maxwell_traction_result rotated_traction;
    swa_local_em_result local;

    check_true(
        swa_maxwell_traction(
            swa_v3(1.0e5, 0.0, 0.0),
            swa_v3(0.0, 0.0, 0.0),
            swa_v3(1.0, 0.0, 0.0),
            &traction
        ),
        "electric normal Maxwell traction executes"
    );
    check_true(
        near(
            traction.traction_N_m2.x,
            0.5 * SWA_EPS0 * 1.0e10,
            1.0e-14,
            1.0e-14
        ),
        "electric normal pressure matches analytic value"
    );
    check_true(
        !swa_maxwell_traction(
            swa_v3(1.0, 0.0, 0.0),
            swa_v3(0.0, 1.0, 0.0),
            swa_v3(0.0, 0.0, 0.0),
            &traction
        ),
        "zero surface normal rejected"
    );

    for (i = 1U; i <= 2000U; ++i) {
        const swa_vec3 electric = swa_v3(
            2.0e5 * (swa_halton(i, 2U) - 0.5),
            2.0e5 * (swa_halton(i, 3U) - 0.5),
            2.0e5 * (swa_halton(i, 5U) - 0.5)
        );
        const swa_vec3 magnetic = swa_v3(
            2.0e-6 * (swa_halton(i, 7U) - 0.5),
            2.0e-6 * (swa_halton(i, 11U) - 0.5),
            2.0e-6 * (swa_halton(i, 13U) - 0.5)
        );
        const swa_vec3 normal = swa_v3(
            0.1 + swa_halton(i, 17U),
            0.2 + swa_halton(i, 19U),
            0.3 + swa_halton(i, 23U)
        );
        const double angle = 2.0 * SWA_PI * swa_halton(i, 29U);
        const swa_vec3 rotated_electric = swa_rotate_z(electric, angle);
        const swa_vec3 rotated_magnetic = swa_rotate_z(magnetic, angle);
        const swa_vec3 rotated_normal = swa_rotate_z(normal, angle);

        check_true(
            swa_local_electromagnetic_state(
                electric, magnetic, &local
            ),
            "local electromagnetic state executes"
        );
        check_true(local.finite,
                   "local electromagnetic state finite");
        check_true(local.energy_density_J_m3 >= 0.0,
                   "electromagnetic energy density nonnegative");
        check_true(local.dominant_energy_ratio <= 1.0 + 1.0e-12,
                   "Poynting flux obeys dominant-energy bound");
        check_true(
            near_vec(
                local.momentum_density_Ns_m3,
                swa_vscale(local.poynting_flux_W_m2,
                           1.0 / (SWA_C * SWA_C)),
                1.0e-14,
                1.0e-30
            ),
            "field momentum equals Poynting flux over c squared"
        );
        check_true(
            swa_maxwell_traction(
                electric, magnetic, normal, &traction
            ),
            "Maxwell traction executes"
        );
        check_true(
            swa_maxwell_traction(
                rotated_electric, rotated_magnetic,
                rotated_normal, &rotated_traction
            ),
            "rotated Maxwell traction executes"
        );
        check_true(
            near_vec(
                swa_rotate_z(traction.traction_N_m2, angle),
                rotated_traction.traction_N_m2,
                2.0e-12,
                1.0e-18
            ),
            "Maxwell traction rotates covariantly"
        );
    }
}

int main(void) {
    swa_charge_continuity_result continuity;
    swa_field_constraint_result gauss;
    swa_field_constraint_result magnetic_divergence;
    swa_poynting_ledger_result energy;
    swa_field_momentum_result momentum;
    FILE *fp;

    memset(&continuity, 0, sizeof(continuity));
    memset(&gauss, 0, sizeof(gauss));
    memset(&magnetic_divergence, 0, sizeof(magnetic_divergence));
    memset(&energy, 0, sizeof(energy));
    memset(&momentum, 0, sizeof(momentum));

    test_charge_continuity(&continuity);
    test_gauss_law(&gauss);
    test_magnetic_divergence(&magnetic_divergence);
    test_integral_ledgers(&energy, &momentum);
    test_local_electromagnetism();

    fp = fopen("output/field_ledger_receipt.json", "w");
    check_true(fp != NULL, "open field ledger receipt");
    if (fp != NULL) {
        check_true(
            swa_write_field_ledger_receipt(
                fp, &continuity, &gauss, &magnetic_divergence,
                &energy, &momentum
            ),
            "write field ledger receipt"
        );
        (void)fclose(fp);
    }

    printf(
        "spacewind field ledger assurance: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
