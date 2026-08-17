#include "spacewind/math3.h"

#include <float.h>
#include <math.h>
#include <string.h>

sw_vec3 sw_v3(double x, double y, double z) {
    sw_vec3 out = {x, y, z};
    return out;
}

sw_vec3 sw_v3_add(sw_vec3 a, sw_vec3 b) {
    return sw_v3(a.x + b.x, a.y + b.y, a.z + b.z);
}

sw_vec3 sw_v3_sub(sw_vec3 a, sw_vec3 b) {
    return sw_v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

sw_vec3 sw_v3_scale(sw_vec3 a, double s) {
    return sw_v3(a.x * s, a.y * s, a.z * s);
}

sw_vec3 sw_v3_mul_add(sw_vec3 a, double s, sw_vec3 b) {
    return sw_v3(a.x + s * b.x, a.y + s * b.y, a.z + s * b.z);
}

double sw_v3_dot(sw_vec3 a, sw_vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

sw_vec3 sw_v3_cross(sw_vec3 a, sw_vec3 b) {
    return sw_v3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

double sw_v3_norm2(sw_vec3 a) {
    return sw_v3_dot(a, a);
}

double sw_v3_norm(sw_vec3 a) {
    return sqrt(sw_v3_norm2(a));
}

sw_vec3 sw_v3_normalize(sw_vec3 a) {
    const double n = sw_v3_norm(a);
    if (n <= DBL_MIN) {
        return sw_v3(0.0, 0.0, 0.0);
    }
    return sw_v3_scale(a, 1.0 / n);
}

sw_vec3 sw_v3_project(sw_vec3 a, sw_vec3 onto) {
    const double d = sw_v3_norm2(onto);
    if (d <= DBL_MIN) {
        return sw_v3(0.0, 0.0, 0.0);
    }
    return sw_v3_scale(onto, sw_v3_dot(a, onto) / d);
}

sw_vec3 sw_v3_reject(sw_vec3 a, sw_vec3 from) {
    return sw_v3_sub(a, sw_v3_project(a, from));
}

bool sw_v3_isfinite(sw_vec3 a) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
}

sw_quat sw_q(double w, double x, double y, double z) {
    sw_quat out = {w, x, y, z};
    return out;
}

sw_quat sw_q_identity(void) {
    return sw_q(1.0, 0.0, 0.0, 0.0);
}

sw_quat sw_q_add(sw_quat a, sw_quat b) {
    return sw_q(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);
}

sw_quat sw_q_scale(sw_quat q, double s) {
    return sw_q(q.w * s, q.x * s, q.y * s, q.z * s);
}

sw_quat sw_q_mul(sw_quat a, sw_quat b) {
    return sw_q(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    );
}

sw_quat sw_q_conj(sw_quat q) {
    return sw_q(q.w, -q.x, -q.y, -q.z);
}

double sw_q_norm(sw_quat q) {
    return sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
}

sw_quat sw_q_normalize(sw_quat q) {
    const double n = sw_q_norm(q);
    if (n <= DBL_MIN) {
        return sw_q_identity();
    }
    return sw_q_scale(q, 1.0 / n);
}

sw_vec3 sw_q_rotate(sw_quat q, sw_vec3 v) {
    q = sw_q_normalize(q);
    const sw_quat p = sw_q(0.0, v.x, v.y, v.z);
    const sw_quat r = sw_q_mul(sw_q_mul(q, p), sw_q_conj(q));
    return sw_v3(r.x, r.y, r.z);
}

sw_quat sw_q_from_axis_angle(sw_vec3 axis, double angle_rad) {
    const sw_vec3 n = sw_v3_normalize(axis);
    const double h = 0.5 * angle_rad;
    const double s = sin(h);
    return sw_q_normalize(sw_q(cos(h), n.x * s, n.y * s, n.z * s));
}

sw_quat sw_q_derivative_body_rate(sw_quat q, sw_vec3 omega_body) {
    const sw_quat omega = sw_q(0.0, omega_body.x, omega_body.y, omega_body.z);
    return sw_q_scale(sw_q_mul(q, omega), 0.5);
}

bool sw_mat4_inverse(const double *a, double *inv_out) {
    double aug[4][8];
    size_t i;
    size_t j;
    size_t col;

    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            aug[i][j] = a[i * 4U + j];
            aug[i][j + 4U] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (col = 0U; col < 4U; ++col) {
        size_t pivot = col;
        double max_abs = fabs(aug[col][col]);
        for (i = col + 1U; i < 4U; ++i) {
            const double candidate = fabs(aug[i][col]);
            if (candidate > max_abs) {
                max_abs = candidate;
                pivot = i;
            }
        }
        if (max_abs <= DBL_EPSILON) {
            return false;
        }
        if (pivot != col) {
            for (j = 0U; j < 8U; ++j) {
                const double tmp = aug[col][j];
                aug[col][j] = aug[pivot][j];
                aug[pivot][j] = tmp;
            }
        }

        {
            const double inv_pivot = 1.0 / aug[col][col];
            for (j = 0U; j < 8U; ++j) {
                aug[col][j] *= inv_pivot;
            }
        }

        for (i = 0U; i < 4U; ++i) {
            if (i == col) {
                continue;
            }
            {
                const double f = aug[i][col];
                for (j = 0U; j < 8U; ++j) {
                    aug[i][j] -= f * aug[col][j];
                }
            }
        }
    }

    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            inv_out[i * 4U + j] = aug[i][j + 4U];
        }
    }
    return true;
}

void sw_mat4_mul(const double *a, const double *b, double *out) {
    double tmp[4][4];
    size_t i;
    size_t j;
    size_t k;
    memset(tmp, 0, sizeof(tmp));
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            for (k = 0U; k < 4U; ++k) {
                tmp[i][j] += a[i * 4U + k] * b[k * 4U + j];
            }
        }
    }
    memcpy(out, tmp, sizeof(tmp));
}
