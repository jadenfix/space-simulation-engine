#include "spacewind/assurance.h"

#include <math.h>
#include <string.h>

static double swa_down(double x) { return nextafter(x, -INFINITY); }
static double swa_up(double x) { return nextafter(x, INFINITY); }
static double swa_max2(double a, double b) { return a > b ? a : b; }
static double swa_max4(double a, double b, double c, double d) {
    return swa_max2(swa_max2(a, b), swa_max2(c, d));
}
static int swa_all_finite3(swa_vec3 v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

swa_vec3 swa_v3(double x, double y, double z) { swa_vec3 v = {x, y, z}; return v; }
swa_vec3 swa_vadd(swa_vec3 a, swa_vec3 b) { return swa_v3(a.x+b.x, a.y+b.y, a.z+b.z); }
swa_vec3 swa_vsub(swa_vec3 a, swa_vec3 b) { return swa_v3(a.x-b.x, a.y-b.y, a.z-b.z); }
swa_vec3 swa_vscale(swa_vec3 a, double s) { return swa_v3(a.x*s, a.y*s, a.z*s); }
double swa_vdot(swa_vec3 a, swa_vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
swa_vec3 swa_vcross(swa_vec3 a, swa_vec3 b) {
    return swa_v3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}
double swa_vnorm(swa_vec3 a) { return sqrt(swa_vdot(a, a)); }
swa_vec3 swa_rotate_z(swa_vec3 a, double angle_rad) {
    const double c = cos(angle_rad);
    const double s = sin(angle_rad);
    return swa_v3(c*a.x-s*a.y, s*a.x+c*a.y, a.z);
}

swa_dimension swa_dim(int8_t mass, int8_t length, int8_t time, int8_t current,
                      int8_t temperature, int8_t amount, int8_t luminous) {
    swa_dimension d = {{mass, length, time, current, temperature, amount, luminous}};
    return d;
}
int swa_dim_equal(swa_dimension a, swa_dimension b) {
    return memcmp(a.exponent, b.exponent, sizeof(a.exponent)) == 0;
}
swa_dimension swa_dim_mul(swa_dimension a, swa_dimension b) {
    swa_dimension out;
    size_t i;
    for (i=0; i<SWA_DIM_COUNT; ++i) out.exponent[i]=(int8_t)(a.exponent[i]+b.exponent[i]);
    return out;
}
swa_dimension swa_dim_div(swa_dimension a, swa_dimension b) {
    swa_dimension out;
    size_t i;
    for (i=0; i<SWA_DIM_COUNT; ++i) out.exponent[i]=(int8_t)(a.exponent[i]-b.exponent[i]);
    return out;
}
swa_quantity swa_quantity_make(double value, swa_dimension dimension) {
    swa_quantity q = {value, dimension}; return q;
}
int swa_quantity_add(swa_quantity a, swa_quantity b, swa_quantity *out) {
    if (out==NULL || !swa_dim_equal(a.dimension,b.dimension)) return 0;
    *out=swa_quantity_make(a.value+b.value,a.dimension); return isfinite(out->value);
}
swa_quantity swa_quantity_mul(swa_quantity a, swa_quantity b) {
    return swa_quantity_make(a.value*b.value,swa_dim_mul(a.dimension,b.dimension));
}
swa_quantity swa_quantity_div(swa_quantity a, swa_quantity b) {
    return swa_quantity_make(a.value/b.value,swa_dim_div(a.dimension,b.dimension));
}

swa_interval swa_interval_make(double a, double b) {
    swa_interval out;
    if (a<=b) { out.lo=a; out.hi=b; } else { out.lo=b; out.hi=a; }
    return out;
}
swa_interval swa_interval_add(swa_interval a, swa_interval b) {
    return swa_interval_make(swa_down(a.lo+b.lo),swa_up(a.hi+b.hi));
}
swa_interval swa_interval_sub(swa_interval a, swa_interval b) {
    return swa_interval_make(swa_down(a.lo-b.hi),swa_up(a.hi-b.lo));
}
swa_interval swa_interval_mul(swa_interval a, swa_interval b) {
    const double p0=a.lo*b.lo, p1=a.lo*b.hi, p2=a.hi*b.lo, p3=a.hi*b.hi;
    const double lo=fmin(fmin(p0,p1),fmin(p2,p3));
    const double hi=fmax(fmax(p0,p1),fmax(p2,p3));
    return swa_interval_make(swa_down(lo),swa_up(hi));
}
int swa_interval_div(swa_interval a, swa_interval b, swa_interval *out) {
    swa_interval reciprocal;
    if (out==NULL || (b.lo<=0.0 && b.hi>=0.0)) return 0;
    reciprocal=swa_interval_make(swa_down(1.0/b.hi),swa_up(1.0/b.lo));
    *out=swa_interval_mul(a,reciprocal); return 1;
}
int swa_interval_sqrt(swa_interval a, swa_interval *out) {
    if (out==NULL || a.lo<0.0) return 0;
    *out=swa_interval_make(swa_down(sqrt(a.lo)),swa_up(sqrt(a.hi))); return 1;
}
int swa_interval_contains(swa_interval a, double x) { return x>=a.lo && x<=a.hi; }
double swa_interval_width(swa_interval a) { return a.hi-a.lo; }

void swa_kahan_add(swa_kahan *accumulator, double value) {
    const double y=value-accumulator->correction;
    const double t=accumulator->sum+y;
    accumulator->correction=(t-accumulator->sum)-y;
    accumulator->sum=t;
}
void swa_stats_reset(swa_online_stats *stats) {
    if (stats==NULL) return;
    stats->count=0; stats->mean=0.0; stats->m2=0.0; stats->minimum=INFINITY; stats->maximum=-INFINITY;
}
void swa_stats_push(swa_online_stats *stats, double value) {
    double delta, delta2;
    if (stats==NULL || !isfinite(value)) return;
    stats->count++;
    delta=value-stats->mean;
    stats->mean+=delta/(double)stats->count;
    delta2=value-stats->mean;
    stats->m2+=delta*delta2;
    if (value<stats->minimum) stats->minimum=value;
    if (value>stats->maximum) stats->maximum=value;
}
double swa_stats_variance(const swa_online_stats *stats) {
    return (stats!=NULL && stats->count>1) ? stats->m2/(double)(stats->count-1) : 0.0;
}
double swa_halton(size_t index, unsigned base) {
    double f=1.0, r=0.0;
    size_t i=index;
    if (base<2U) return 0.0;
    while (i>0U) { f/=(double)base; r+=f*(double)(i%(size_t)base); i/=(size_t)base; }
    return r;
}
uint64_t swa_fnv1a64(const void *data, size_t size) {
    const unsigned char *p=(const unsigned char *)data;
    uint64_t h=UINT64_C(14695981039346656037);
    size_t i;
    for (i=0;i<size;++i) { h^=(uint64_t)p[i]; h*=UINT64_C(1099511628211); }
    return h;
}

int swa_audit_ledger(const swa_ledger *ledger, double relative_tolerance,
                     double absolute_tolerance, swa_ledger_result *out) {
    double mass_scale, charge_scale, momentum_scale, energy_scale;
    swa_vec3 expected_momentum;
    if (ledger==NULL || out==NULL || relative_tolerance<0.0 || absolute_tolerance<0.0) return 0;
    memset(out,0,sizeof(*out));
    out->mass_residual_kg=ledger->mass_final_kg-(ledger->mass_initial_kg+ledger->mass_in_kg-ledger->mass_out_kg);
    out->charge_residual_C=ledger->charge_final_C-(ledger->charge_initial_C+ledger->charge_in_C-ledger->charge_out_C);
    expected_momentum=swa_vadd(ledger->momentum_initial_Ns,
        swa_vadd(ledger->impulse_external_Ns,swa_vsub(ledger->momentum_in_Ns,ledger->momentum_out_Ns)));
    out->momentum_residual_Ns=swa_vsub(ledger->momentum_final_Ns,expected_momentum);
    out->momentum_residual_norm_Ns=swa_vnorm(out->momentum_residual_Ns);
    out->energy_residual_J=ledger->energy_final_J-(ledger->energy_initial_J+ledger->work_source_J+
        ledger->heat_in_J+ledger->energy_in_J-ledger->energy_out_J);
    mass_scale=swa_max4(fabs(ledger->mass_initial_kg),fabs(ledger->mass_final_kg),fabs(ledger->mass_in_kg),fabs(ledger->mass_out_kg));
    charge_scale=swa_max4(fabs(ledger->charge_initial_C),fabs(ledger->charge_final_C),fabs(ledger->charge_in_C),fabs(ledger->charge_out_C));
    momentum_scale=swa_max4(swa_vnorm(ledger->momentum_initial_Ns),swa_vnorm(ledger->momentum_final_Ns),
                            swa_vnorm(ledger->impulse_external_Ns),swa_max2(swa_vnorm(ledger->momentum_in_Ns),swa_vnorm(ledger->momentum_out_Ns)));
    energy_scale=swa_max4(fabs(ledger->energy_initial_J),fabs(ledger->energy_final_J),fabs(ledger->work_source_J),
                          swa_max2(fabs(ledger->energy_in_J),fabs(ledger->energy_out_J)));
    out->mass_relative=fabs(out->mass_residual_kg)/swa_max2(mass_scale,absolute_tolerance);
    out->charge_relative=fabs(out->charge_residual_C)/swa_max2(charge_scale,absolute_tolerance);
    out->momentum_relative=out->momentum_residual_norm_Ns/swa_max2(momentum_scale,absolute_tolerance);
    out->energy_relative=fabs(out->energy_residual_J)/swa_max2(energy_scale,absolute_tolerance);
    out->finite=isfinite(out->mass_relative)&&isfinite(out->charge_relative)&&
                swa_all_finite3(out->momentum_residual_Ns)&&isfinite(out->energy_relative);
    out->passes=out->finite &&
        (fabs(out->mass_residual_kg)<=absolute_tolerance || out->mass_relative<=relative_tolerance) &&
        (fabs(out->charge_residual_C)<=absolute_tolerance || out->charge_relative<=relative_tolerance) &&
        (out->momentum_residual_norm_Ns<=absolute_tolerance || out->momentum_relative<=relative_tolerance) &&
        (fabs(out->energy_residual_J)<=absolute_tolerance || out->energy_relative<=relative_tolerance);
    return 1;
}

int swa_richardson_three_level(double coarse, double medium, double fine,
                               double refinement_ratio, double expected_order,
                               double order_tolerance, swa_convergence_result *out) {
    double dcm, dmf, ratio, denom;
    if (out==NULL || !(refinement_ratio>1.0) || !(order_tolerance>=0.0)) return 0;
    memset(out,0,sizeof(*out));
    dcm=coarse-medium; dmf=medium-fine;
    out->monotone=(dcm*dmf)>0.0;
    if (dmf==0.0 || dcm/dmf<=0.0) return 1;
    ratio=fabs(dcm/dmf);
    out->observed_order=log(ratio)/log(refinement_ratio);
    denom=pow(refinement_ratio,out->observed_order)-1.0;
    if (denom==0.0) return 1;
    out->extrapolated_value=fine+(fine-medium)/denom;
    out->fine_error_estimate=fabs(out->extrapolated_value-fine);
    out->gci_fine=1.25*fabs(fine-medium)/fabs(denom);
    out->asymptotic_ratio=(out->gci_fine>0.0) ?
        (1.25*fabs(medium-coarse)/fabs(denom))/(pow(refinement_ratio,out->observed_order)*out->gci_fine) : 1.0;
    out->finite=isfinite(out->observed_order)&&isfinite(out->extrapolated_value)&&
                isfinite(out->gci_fine)&&isfinite(out->asymptotic_ratio);
    out->passes=out->finite&&out->monotone&&fabs(out->observed_order-expected_order)<=order_tolerance&&
                fabs(out->asymptotic_ratio-1.0)<=0.25;
    return 1;
}

int swa_write_assurance_receipt(FILE *fp, size_t checks, size_t failures,
                                uint64_t deterministic_hash) {
    if (fp==NULL) return 0;
    return fprintf(fp,
        "{\n  \"schema\": \"spacewind.assurance/v1\",\n  \"checks\": %zu,\n"
        "  \"failures\": %zu,\n  \"deterministic_hash\": \"%016llx\",\n"
        "  \"claim\": \"finite deterministic assurance checks only\",\n"
        "  \"nonclaim\": \"passing checks do not establish physical propulsion or net energy extraction\"\n}\n",
        checks,failures,(unsigned long long)deterministic_hash)>0;
}
