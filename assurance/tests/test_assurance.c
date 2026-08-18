#include "spacewind/assurance.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks=0,failures=0;
static uint64_t transcript_hash=UINT64_C(14695981039346656037);

static void mix_u64(uint64_t x) {
    transcript_hash^=x;
    transcript_hash*=UINT64_C(1099511628211);
}
static void check_true(int condition,const char *name) {
    ++checks; mix_u64((uint64_t)checks); mix_u64(swa_fnv1a64(name,strlen(name)));
    if (!condition) {++failures;fprintf(stderr,"FAIL: %s\n",name);}
}
static int near(double a,double b,double rel,double abs_tol) {
    const double scale=fmax(fabs(a),fabs(b));
    return fabs(a-b)<=fmax(abs_tol,rel*scale);
}

static swa_plasma_state nominal_plasma(void) {
    swa_plasma_state s;
    memset(&s,0,sizeof(s));
    s.ion_density_m3=5.0e6; s.electron_density_m3=5.0e6;
    s.ion_temperature_K=1.0e5; s.electron_temperature_K=1.5e5;
    s.ion_mass_kg=SWA_MP; s.ion_charge_C=SWA_QE;
    s.neutral_density_m3=1.0; s.collision_cross_section_m2=1.0e-19;
    s.electrode_radius_m=25.0e-6; s.interaction_length_m=100.0;
    s.electrode_voltage_V=20000.0; s.circuit_capacitance_F=1.0e-6;
    s.available_power_W=1000.0; s.collection_current_A=0.02;
    s.bulk_velocity_mps=swa_v3(4.0e5,0.0,0.0);
    s.magnetic_field_T=swa_v3(0.0,0.0,5.0e-9);
    return s;
}

static void test_dimensions(void) {
    const swa_dimension mass=swa_dim(1,0,0,0,0,0,0);
    const swa_dimension length=swa_dim(0,1,0,0,0,0,0);
    const swa_dimension time=swa_dim(0,0,1,0,0,0,0);
    const swa_dimension velocity=swa_dim_div(length,time);
    const swa_dimension acceleration=swa_dim_div(velocity,time);
    const swa_dimension force=swa_dim_mul(mass,acceleration);
    const swa_dimension energy=swa_dim_mul(force,length);
    swa_quantity q;
    check_true(swa_dim_equal(force,swa_dim(1,1,-2,0,0,0,0)),"force dimension");
    check_true(swa_dim_equal(energy,swa_dim(1,2,-2,0,0,0,0)),"energy dimension");
    check_true(swa_quantity_add(swa_quantity_make(2.0,force),swa_quantity_make(3.0,force),&q)&&q.value==5.0,"same-dimension add");
    check_true(!swa_quantity_add(swa_quantity_make(2.0,force),swa_quantity_make(3.0,energy),&q),"reject mixed-dimension add");
    check_true(swa_dim_equal(swa_quantity_mul(swa_quantity_make(2.0,mass),swa_quantity_make(3.0,acceleration)).dimension,force),"quantity multiplication dimensions");
}

static void test_intervals(void) {
    size_t i;
    for (i=1;i<=2000;++i) {
        const double a=-10.0+20.0*swa_halton(i,2U);
        const double b=-10.0+20.0*swa_halton(i,3U);
        const double c=0.1+9.9*swa_halton(i,5U);
        const double d=0.1+9.9*swa_halton(i,7U);
        const swa_interval x=swa_interval_make(a,b);
        const swa_interval y=swa_interval_make(c,d);
        const double sx=x.lo+(x.hi-x.lo)*swa_halton(i,11U);
        const double sy=y.lo+(y.hi-y.lo)*swa_halton(i,13U);
        swa_interval q,root;
        check_true(swa_interval_contains(swa_interval_add(x,y),sx+sy),"interval add containment");
        check_true(swa_interval_contains(swa_interval_sub(x,y),sx-sy),"interval sub containment");
        check_true(swa_interval_contains(swa_interval_mul(x,y),sx*sy),"interval mul containment");
        check_true(swa_interval_div(x,y,&q)&&swa_interval_contains(q,sx/sy),"interval div containment");
        check_true(swa_interval_sqrt(y,&root)&&swa_interval_contains(root,sqrt(sy)),"interval sqrt containment");
    }
    { swa_interval q; check_true(!swa_interval_div(swa_interval_make(1,2),swa_interval_make(-1,1),&q),"division rejects zero-containing denominator"); }
}

static void test_numerics(void) {
    size_t i;
    swa_kahan k={0.0,0.0};
    swa_online_stats stats;
    swa_convergence_result c;
    double naive=0.0;
    swa_stats_reset(&stats);
    for (i=0;i<100000;++i) { naive+=1e-6; swa_kahan_add(&k,1e-6); swa_stats_push(&stats,(double)i); }
    check_true(fabs(k.sum-0.1)<=fabs(naive-0.1)+1e-15,"Kahan no worse than naive");
    check_true(stats.count==100000,"online stats count");
    check_true(near(stats.mean,49999.5,1e-15,1e-12),"online stats mean");
    check_true(swa_richardson_three_level(1.0+0.25,1.0+0.0625,1.0+0.015625,2.0,2.0,0.05,&c),"Richardson call");
    check_true(c.passes,"second-order convergence passes");
    check_true(near(c.extrapolated_value,1.0,1e-12,1e-12),"Richardson extrapolation");
    check_true(swa_richardson_three_level(1.4,1.1,1.3,2.0,2.0,0.1,&c)&&!c.passes,"nonmonotone convergence rejected");
}

static void test_ledgers(void) {
    swa_ledger l;
    swa_ledger_result r;
    memset(&l,0,sizeof(l));
    l.mass_initial_kg=10.0;l.mass_in_kg=2.0;l.mass_out_kg=1.0;l.mass_final_kg=11.0;
    l.charge_initial_C=1.0;l.charge_in_C=0.5;l.charge_out_C=0.25;l.charge_final_C=1.25;
    l.momentum_initial_Ns=swa_v3(1,2,3);l.impulse_external_Ns=swa_v3(2,-1,1);
    l.momentum_in_Ns=swa_v3(0.5,0.25,0);l.momentum_out_Ns=swa_v3(0.25,0,0.5);
    l.momentum_final_Ns=swa_v3(3.25,1.25,3.5);
    l.energy_initial_J=100;l.work_source_J=20;l.heat_in_J=5;l.energy_in_J=10;l.energy_out_J=15;l.energy_final_J=120;
    check_true(swa_audit_ledger(&l,1e-12,1e-12,&r)&&r.passes,"closed multi-ledger passes");
    l.energy_final_J+=0.1;
    check_true(swa_audit_ledger(&l,1e-12,1e-12,&r)&&!r.passes,"energy leak detected");
    l.energy_final_J-=0.1;l.momentum_final_Ns.z+=1e-3;
    check_true(swa_audit_ledger(&l,1e-12,1e-12,&r)&&!r.passes,"momentum leak detected");
}

static void test_covariance_and_scaling(void) {
    size_t i;
    for (i=1;i<=1000;++i) {
        const swa_vec3 a=swa_v3(-2+4*swa_halton(i,2),-2+4*swa_halton(i,3),-2+4*swa_halton(i,5));
        const swa_vec3 b=swa_v3(-2+4*swa_halton(i,7),-2+4*swa_halton(i,11),-2+4*swa_halton(i,13));
        const double theta=2*SWA_PI*swa_halton(i,17);
        const swa_vec3 lhs=swa_rotate_z(swa_vcross(a,b),theta);
        const swa_vec3 rhs=swa_vcross(swa_rotate_z(a,theta),swa_rotate_z(b,theta));
        check_true(swa_vnorm(swa_vsub(lhs,rhs))<1e-12,"rotation covariance of cross product");
        check_true(near(swa_vdot(a,b),swa_vdot(swa_rotate_z(a,theta),swa_rotate_z(b,theta)),1e-12,1e-12),"rotation invariance of dot product");
        {
            const double r1=0.2*SWA_AU+5*SWA_AU*swa_halton(i,19);
            const double r2=0.2*SWA_AU+5*SWA_AU*swa_halton(i,23);
            const double p1=SWA_SOLAR_IRRADIANCE_1AU*(SWA_AU/r1)*(SWA_AU/r1)/SWA_C;
            const double p2=SWA_SOLAR_IRRADIANCE_1AU*(SWA_AU/r2)*(SWA_AU/r2)/SWA_C;
            check_true(near(p1/p2,(r2/r1)*(r2/r1),1e-13,1e-13),"photon pressure inverse-square metamorphism");
            check_true(near((SWA_G*1.98847e30/(r1*r1))/(SWA_G*1.98847e30/(r2*r2)),(r2/r1)*(r2/r1),1e-13,1e-13),"gravity inverse-square metamorphism");
        }
    }
}

static void test_plasma_similarity(void) {
    swa_plasma_state s=nominal_plasma();
    swa_plasma_regime r1,r2;
    swa_similarity_signature a,b;
    swa_similarity_result result;
    check_true(swa_compute_plasma_regime(&s,&r1)&&r1.finite,"nominal plasma regime finite");
    check_true(r1.debye_length_m>0&&r1.ion_inertial_length_m>0&&r1.fast_magnetosonic_speed_mps>0,"plasma scales positive");
    swa_make_similarity_signature(&r1,3,&a);
    b=a;
    check_true(swa_compare_similarity(&a,&b,3,&result)&&result.inside&&result.weighted_log_distance==0.0,"identical similarity accepted");
    s.electron_temperature_K*=100.0;
    check_true(swa_compute_plasma_regime(&s,&r2),"mismatched regime computes");
    swa_make_similarity_signature(&r2,3,&b);
    check_true(swa_compare_similarity(&a,&b,3,&result)&&!result.inside,"large similarity mismatch rejected");
    b=a;b.evidence_level=1;
    check_true(swa_compare_similarity(&a,&b,3,&result)&&!result.inside,"insufficient evidence rejected");
}

static void test_cycles(void) {
    swa_cycle_audit a;
    swa_cycle_result r;
    FILE *fp;
    memset(&a,0,sizeof(a));
    a.start.position_m=swa_v3(1e6,2e6,0);a.end.position_m=a.start.position_m;
    a.start.velocity_mps=swa_v3(1000,2000,0);a.end.velocity_mps=a.start.velocity_mps;
    a.start.stored_energy_J=100;a.end.stored_energy_J=110;
    a.start.controller_phase_rad=0.1;a.end.controller_phase_rad=0.1+2*SWA_PI;
    a.craft_work_J=20;a.source_flow_work_J=25;a.relative_dissipation_J=5;
    a.actuator_energy_J=5;a.thermal_loss_J=5;a.numerical_residual_J=0;
    a.uncertain_net_gain_J=swa_interval_make(8,12);
    a.position_scale_m=1e6;a.velocity_scale_mps=1000;a.energy_scale_J=100;a.phase_scale_rad=1;
    a.evidence_level=2;a.required_evidence_level=4;
    check_true(swa_audit_cycle(&a,1e-12,1e-12,&r),"cycle audit call");
    check_true(r.numerical_cycle_passes,"closed numerical cycle passes");
    check_true(!r.physical_promotion_passes,"low-evidence cycle not promoted");
    a.evidence_level=4;
    check_true(swa_audit_cycle(&a,1e-12,1e-12,&r)&&r.physical_promotion_passes,"robust high-evidence cycle promotes");
    a.uncertain_net_gain_J=swa_interval_make(-1,12);
    check_true(swa_audit_cycle(&a,1e-12,1e-12,&r)&&!r.physical_promotion_passes,"uncertainty crossing zero blocks promotion");
    a.uncertain_net_gain_J=swa_interval_make(8,12);a.end.position_m.x+=100;
    check_true(swa_audit_cycle(&a,1e-8,1e-12,&r)&&!r.numerical_cycle_passes,"state nonclosure detected");
    a.end.position_m=a.start.position_m;
    fp=fopen("output/cycle_receipt.json","w");
    check_true(fp!=NULL,"open cycle receipt");
    if (fp!=NULL) { swa_audit_cycle(&a,1e-12,1e-12,&r);check_true(swa_write_cycle_receipt(fp,&a,&r),"write cycle receipt");fclose(fp); }
}

static void test_faults(void) {
    swa_fault_state s;
    swa_fault_result r;
    memset(&s,0,sizeof(s));
    s.temperature_K=300;s.temperature_limit_K=500;s.bus_voltage_V=100;s.minimum_bus_voltage_V=50;
    s.maximum_angular_rate_rad_s=1.0;
    check_true(swa_assess_faults(&s,&r)&&!r.safe_mode_required,"nominal fault state stays active");
    s.arc_energy_J=10;
    check_true(swa_assess_faults(&s,&r)&&r.arc_detected&&r.safe_mode_required,"arc forces safe mode");
    s.arc_energy_J=0;s.temperature_K=600;
    check_true(swa_assess_faults(&s,&r)&&r.thermal_trip&&r.safe_mode_required,"thermal trip forces safe mode");
    s.temperature_K=300;s.bus_voltage_V=20;
    check_true(swa_assess_faults(&s,&r)&&r.undervoltage&&r.safe_mode_required,"undervoltage forces safe mode");
    s.bus_voltage_V=100;s.angular_rate_rad_s=2;
    check_true(swa_assess_faults(&s,&r)&&r.attitude_trip&&r.safe_mode_required,"attitude trip forces safe mode");
    s.angular_rate_rad_s=0;s.dropout_fraction=0.5;s.tether_cut_fraction=0.1;
    check_true(swa_assess_faults(&s,&r)&&r.safe_mode_required,"combined comms and tether degradation forces safe mode");
}

static void test_determinism(void) {
    size_t i;
    uint64_t h1=UINT64_C(14695981039346656037),h2=UINT64_C(14695981039346656037);
    for (i=1;i<=1000;++i) {
        const double values[3]={swa_halton(i,2),swa_halton(i,3),swa_halton(i,5)};
        h1^=swa_fnv1a64(values,sizeof(values));h1*=UINT64_C(1099511628211);
    }
    for (i=1;i<=1000;++i) {
        const double values[3]={swa_halton(i,2),swa_halton(i,3),swa_halton(i,5)};
        h2^=swa_fnv1a64(values,sizeof(values));h2*=UINT64_C(1099511628211);
    }
    check_true(h1==h2,"deterministic low-discrepancy transcript");
}

int main(void) {
    FILE *fp;
    test_dimensions();test_intervals();test_numerics();test_ledgers();
    test_covariance_and_scaling();test_plasma_similarity();test_cycles();test_faults();test_determinism();
    fp=fopen("output/assurance_receipt.json","w");
    if (fp!=NULL) { swa_write_assurance_receipt(fp,checks,failures,transcript_hash);fclose(fp); }
    printf("spacewind assurance: %zu checks, %zu failures, hash=%016llx\n",
           checks,failures,(unsigned long long)transcript_hash);
    return failures==0?EXIT_SUCCESS:EXIT_FAILURE;
}
