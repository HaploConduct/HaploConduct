#include "EditDistance.h"

signed char MyersEditDistanceIncremental::myers_total[256 * 256];
signed char MyersEditDistanceIncremental::myers_optimum[256 * 256];
uchar       MyersEditDistanceIncremental::myers_best_pos[256 * 256];


