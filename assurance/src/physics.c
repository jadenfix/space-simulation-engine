#include "spacewind/assurance.h"

#include <math.h>
#include <string.h>

static double max2(double a,double b){return a>b?a:b;}
static double clamp01(double x){return x<0.0?0.0:(x>1.0?1.0:x);}
static double safe_scale(double x,double floor_value){return max2(fabs(x),floor_value);}
static double positive_log_ratio(double a,double b){
    if (!(a>0.0) || !(b>0.0)) return INFINITY;
    return fabs(log(a/b));
}

int swa_compute_plasma_regime(const swa_plasma_state *s, swa_plasma_regime *o) {
    double rho, bmag, speed, thermal_pressure, magnetic_pressure, collision_frequency;
    double flow_time, circuit_time, required_power;
    if (s==NULL || o==NULL || !(s->ion_density_m3>0.0) || !(s->electron_density_m3>0.0) ||
        !(s->ion_temperature_K>0.0) || !(s->electron_temperature_K>0.0) ||
        !(s->ion_mass_kg>0.0) || s->ion_charge_C==0.0 || !(s->interaction_length_m>0.0)) return 0;
    memset(o,0,sizeof(*o));
    rho=s->ion_density_m3*s->ion_mass_kg;
    bmag=swa_vnorm(s->magnetic_field_T);
    speed=swa_vnorm(s->bulk_velocity_mps);
    o->debye_length_m=sqrt(SWA_EPS0*SWA_KB*s->electron_temperature_K/
        (s->electron_density_m3*s->ion_charge_C*s->ion_charge_C));
    o->ion_plasma_frequency_rad_s=sqrt(s->ion_density_m3*s->ion_charge_C*s->ion_charge_C/
        (SWA_EPS0*s->ion_mass_kg));
    o->electron_plasma_frequency_rad_s=sqrt(s->electron_density_m3*SWA_QE*SWA_QE/(SWA_EPS0*SWA_ME));
    o->ion_inertial_length_m=SWA_C/o->ion_plasma_frequency_rad_s;
    o->ion_gyrofrequency_rad_s=fabs(s->ion_charge_C)*bmag/s->ion_mass_kg;
    o->ion_thermal_speed_mps=sqrt(2.0*SWA_KB*s->ion_temperature_K/s->ion_mass_kg);
    o->ion_gyroradius_m=(o->ion_gyrofrequency_rad_s>0.0)?o->ion_thermal_speed_mps/o->ion_gyrofrequency_rad_s:INFINITY;
    o->sound_speed_mps=sqrt(SWA_KB*(s->electron_temperature_K+3.0*s->ion_temperature_K)/s->ion_mass_kg);
    o->alfven_speed_mps=(rho>0.0)?bmag/sqrt(SWA_MU0*rho):INFINITY;
    o->fast_magnetosonic_speed_mps=hypot(o->sound_speed_mps,o->alfven_speed_mps);
    thermal_pressure=SWA_KB*(s->ion_density_m3*s->ion_temperature_K+s->electron_density_m3*s->electron_temperature_K);
    magnetic_pressure=bmag*bmag/(2.0*SWA_MU0);
    o->plasma_beta=(magnetic_pressure>0.0)?thermal_pressure/magnetic_pressure:INFINITY;
    o->sonic_mach=speed/o->sound_speed_mps;
    o->alfven_mach=(o->alfven_speed_mps>0.0)?speed/o->alfven_speed_mps:0.0;
    o->fast_mach=speed/o->fast_magnetosonic_speed_mps;
    collision_frequency=s->neutral_density_m3*s->collision_cross_section_m2*max2(speed,o->ion_thermal_speed_mps);
    o->knudsen=(collision_frequency>0.0)?speed/(collision_frequency*s->interaction_length_m):INFINITY;
    o->normalized_voltage=fabs(s->ion_charge_C*s->electrode_voltage_V)/(SWA_KB*s->electron_temperature_K);
    o->radius_to_debye=s->electrode_radius_m/o->debye_length_m;
    o->gyroradius_to_length=o->ion_gyroradius_m/s->interaction_length_m;
    o->inertial_to_length=o->ion_inertial_length_m/s->interaction_length_m;
    flow_time=s->interaction_length_m/max2(speed,1.0);
    circuit_time=s->circuit_capacitance_F*fabs(s->electrode_voltage_V)/max2(fabs(s->collection_current_A),1e-30);
    o->circuit_to_flow_time=circuit_time/flow_time;
    required_power=fabs(s->electrode_voltage_V*s->collection_current_A);
    o->power_margin=s->available_power_W/max2(required_power,1e-30);
    o->finite=isfinite(o->debye_length_m)&&isfinite(o->ion_plasma_frequency_rad_s)&&
        isfinite(o->sound_speed_mps)&&isfinite(o->fast_magnetosonic_speed_mps)&&
        isfinite(o->normalized_voltage)&&isfinite(o->radius_to_debye)&&
        isfinite(o->inertial_to_length)&&isfinite(o->circuit_to_flow_time)&&isfinite(o->power_margin);
    return 1;
}

void swa_make_similarity_signature(const swa_plasma_regime *r,uint32_t evidence_level,
                                   swa_similarity_signature *o) {
    static const double default_tolerance[SWA_SIMILARITY_COUNT]={
        0.25,0.25,0.25,0.25,0.35,0.35,0.35,0.35,0.35,0.35,0.50,0.35};
    if (r==NULL||o==NULL) return;
    o->value[0]=r->normalized_voltage;
    o->value[1]=r->radius_to_debye;
    o->value[2]=r->gyroradius_to_length;
    o->value[3]=r->inertial_to_length;
    o->value[4]=r->sonic_mach;
    o->value[5]=r->alfven_mach;
    o->value[6]=r->fast_mach;
    o->value[7]=r->plasma_beta;
    o->value[8]=r->knudsen;
    o->value[9]=r->circuit_to_flow_time;
    o->value[10]=r->power_margin;
    o->value[11]=r->debye_length_m;
    memcpy(o->tolerance_log,default_tolerance,sizeof(default_tolerance));
    o->evidence_level=evidence_level;
}

int swa_compare_similarity(const swa_similarity_signature *target,
                           const swa_similarity_signature *candidate,
                           uint32_t required_evidence_level,
                           swa_similarity_result *out) {
    size_t i;
    double sum=0.0, worst=-1.0;
    if (target==NULL||candidate==NULL||out==NULL) return 0;
    memset(out,0,sizeof(*out));
    out->inside=1;
    for (i=0;i<SWA_SIMILARITY_COUNT;++i) {
        const double mismatch=positive_log_ratio(target->value[i],candidate->value[i]);
        const double tolerance=max2(target->tolerance_log[i],1e-12);
        const double normalized=mismatch/tolerance;
        if (!isfinite(normalized)||normalized>1.0) out->inside=0;
        sum+=normalized*normalized;
        if (mismatch>worst) {worst=mismatch;out->worst_index=i;}
    }
    out->weighted_log_distance=sqrt(sum/(double)SWA_SIMILARITY_COUNT);
    out->maximum_log_mismatch=worst;
    out->minimum_evidence_level=target->evidence_level<candidate->evidence_level?
        target->evidence_level:candidate->evidence_level;
    if (out->minimum_evidence_level<required_evidence_level) out->inside=0;
    return 1;
}

int swa_audit_cycle(const swa_cycle_audit *a,double state_tolerance,
                    double identity_tolerance,swa_cycle_result *o) {
    double position_error,velocity_error,energy_state_error,phase_error;
    double moving_residual,energy_residual,source_scale,energy_scale;
    if (a==NULL||o==NULL||state_tolerance<0.0||identity_tolerance<0.0) return 0;
    memset(o,0,sizeof(*o));
    position_error=swa_vnorm(swa_vsub(a->end.position_m,a->start.position_m))/safe_scale(a->position_scale_m,1.0);
    velocity_error=swa_vnorm(swa_vsub(a->end.velocity_mps,a->start.velocity_mps))/safe_scale(a->velocity_scale_mps,1.0);
    energy_state_error=fabs(a->end.stored_energy_J-a->start.stored_energy_J)/safe_scale(a->energy_scale_J,1.0);
    phase_error=fabs(remainder(a->end.controller_phase_rad-a->start.controller_phase_rad,2.0*SWA_PI))/safe_scale(a->phase_scale_rad,1.0);
    o->normalized_state_closure=max2(max2(position_error,velocity_error),max2(energy_state_error,phase_error));
    moving_residual=a->source_flow_work_J-(a->craft_work_J+a->relative_dissipation_J);
    source_scale=max2(max2(fabs(a->source_flow_work_J),fabs(a->craft_work_J)),max2(fabs(a->relative_dissipation_J),1.0));
    o->moving_medium_identity_relative=fabs(moving_residual)/source_scale;
    energy_residual=a->craft_work_J-a->actuator_energy_J-a->thermal_loss_J-a->numerical_residual_J-
        (a->end.stored_energy_J-a->start.stored_energy_J);
    energy_scale=max2(max2(fabs(a->craft_work_J),fabs(a->actuator_energy_J)),max2(fabs(a->thermal_loss_J),1.0));
    o->energy_identity_relative=fabs(energy_residual)/energy_scale;
    o->net_gain_J=a->craft_work_J-a->actuator_energy_J-a->thermal_loss_J;
    o->net_gain_interval_J=a->uncertain_net_gain_J;
    o->state_closed=o->normalized_state_closure<=state_tolerance;
    o->identity_closed=o->moving_medium_identity_relative<=identity_tolerance&&
                       o->energy_identity_relative<=identity_tolerance;
    o->evidence_sufficient=a->evidence_level>=a->required_evidence_level;
    o->robust_positive_gain=a->uncertain_net_gain_J.lo>0.0;
    o->numerical_cycle_passes=o->state_closed&&o->identity_closed;
    o->physical_promotion_passes=o->numerical_cycle_passes&&o->evidence_sufficient&&o->robust_positive_gain;
    return 1;
}

int swa_assess_faults(const swa_fault_state *s,swa_fault_result *o) {
    double severity=0.0;
    if (s==NULL||o==NULL) return 0;
    memset(o,0,sizeof(*o));
    o->sensor_degraded=fabs(s->sensor_bias_fraction)>0.05;
    o->communications_degraded=s->dropout_fraction>0.20;
    o->tether_degraded=s->tether_cut_fraction>0.02;
    o->arc_detected=s->arc_energy_J>1.0;
    o->thermal_trip=s->temperature_K>=s->temperature_limit_K;
    o->undervoltage=s->bus_voltage_V<=s->minimum_bus_voltage_V;
    o->attitude_trip=fabs(s->angular_rate_rad_s)>=s->maximum_angular_rate_rad_s;
    severity=max2(severity,clamp01(fabs(s->sensor_bias_fraction)/0.20));
    severity=max2(severity,clamp01(s->dropout_fraction/0.50));
    severity=max2(severity,clamp01(s->tether_cut_fraction/0.20));
    severity=max2(severity,clamp01(s->arc_energy_J/100.0));
    if (s->temperature_limit_K>0.0) severity=max2(severity,clamp01(s->temperature_K/s->temperature_limit_K));
    if (s->minimum_bus_voltage_V>0.0) severity=max2(severity,clamp01(1.0-s->bus_voltage_V/s->minimum_bus_voltage_V));
    if (s->maximum_angular_rate_rad_s>0.0) severity=max2(severity,clamp01(fabs(s->angular_rate_rad_s)/s->maximum_angular_rate_rad_s));
    o->severity=severity;
    o->safe_mode_required=o->arc_detected||o->thermal_trip||o->undervoltage||o->attitude_trip||
                          (o->tether_degraded&&o->communications_degraded);
    return 1;
}

int swa_write_cycle_receipt(FILE *fp,const swa_cycle_audit *a,const swa_cycle_result *r) {
    if (fp==NULL||a==NULL||r==NULL) return 0;
    return fprintf(fp,
        "{\n  \"schema\": \"spacewind.cycle-audit/v1\",\n"
        "  \"state_closure\": %.17g,\n  \"moving_medium_identity\": %.17g,\n"
        "  \"energy_identity\": %.17g,\n  \"net_gain_J\": %.17g,\n"
        "  \"net_gain_interval_J\": [%.17g, %.17g],\n"
        "  \"evidence_level\": %u,\n  \"required_evidence_level\": %u,\n"
        "  \"numerical_cycle_passes\": %s,\n  \"physical_promotion_passes\": %s,\n"
        "  \"nonclaim\": \"a numerical closed cycle is not chamber or flight evidence\"\n}\n",
        r->normalized_state_closure,r->moving_medium_identity_relative,r->energy_identity_relative,
        r->net_gain_J,r->net_gain_interval_J.lo,r->net_gain_interval_J.hi,
        a->evidence_level,a->required_evidence_level,
        r->numerical_cycle_passes?"true":"false",r->physical_promotion_passes?"true":"false")>0;
}
