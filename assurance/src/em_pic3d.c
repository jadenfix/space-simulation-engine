#include "spacewind/em_pic3d.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define PIC3D_MAX_PATH_SEGMENTS 4096U

#include "em_pic3d_support.inc"
#include "em_pic3d_deposition.inc"
#include "em_pic3d_spectral.inc"

#define pic3d_electric_constraint_metrics \
    pic3d_electric_constraint_metrics_unconditioned
#define swa_em_pic3d_electric_divergence_max_C_m3 \
    swa_em_pic3d_electric_divergence_max_C_m3_unconditioned
#include "em_pic3d_diagnostics.inc"
#undef swa_em_pic3d_electric_divergence_max_C_m3
#undef pic3d_electric_constraint_metrics

#include "em_pic3d_conditioning.inc"
#include "em_pic3d_step.inc"
