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
#include "em_pic3d_diagnostics.inc"
#include "em_pic3d_step.inc"
