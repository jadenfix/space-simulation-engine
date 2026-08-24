#ifndef SPACEWIND_MATH3_H
#define SPACEWIND_MATH3_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double x;
    double y;
    double z;
} sw_vec3;

typedef struct {
    double w;
    double x;
    double y;
    double z;
} sw_quat;

sw_vec3 sw_v3(double x, double y, double z);
sw_vec3 sw_v3_add(sw_vec3 a, sw_vec3 b);
sw_vec3 sw_v3_sub(sw_vec3 a, sw_vec3 b);
sw_vec3 sw_v3_scale(sw_vec3 a, double s);
sw_vec3 sw_v3_mul_add(sw_vec3 a, double s, sw_vec3 b);
double sw_v3_dot(sw_vec3 a, sw_vec3 b);
sw_vec3 sw_v3_cross(sw_vec3 a, sw_vec3 b);
double sw_v3_norm2(sw_vec3 a);
double sw_v3_norm(sw_vec3 a);
sw_vec3 sw_v3_normalize(sw_vec3 a);
sw_vec3 sw_v3_project(sw_vec3 a, sw_vec3 onto);
sw_vec3 sw_v3_reject(sw_vec3 a, sw_vec3 from);
bool sw_v3_isfinite(sw_vec3 a);

sw_quat sw_q(double w, double x, double y, double z);
sw_quat sw_q_identity(void);
sw_quat sw_q_add(sw_quat a, sw_quat b);
sw_quat sw_q_scale(sw_quat q, double s);
sw_quat sw_q_mul(sw_quat a, sw_quat b);
sw_quat sw_q_conj(sw_quat q);
double sw_q_norm(sw_quat q);
sw_quat sw_q_normalize(sw_quat q);
sw_vec3 sw_q_rotate(sw_quat q, sw_vec3 v);
sw_quat sw_q_from_axis_angle(sw_vec3 axis, double angle_rad);
sw_quat sw_q_derivative_body_rate(sw_quat q, sw_vec3 omega_body);

bool sw_mat4_inverse(const double *a, double *inv_out);
void sw_mat4_mul(const double *a, const double *b, double *out);

#endif
