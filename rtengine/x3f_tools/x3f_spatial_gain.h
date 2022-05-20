/* X3F_SPATIAL_GAIN.H
 *
 * Library for adjusting for spatial gain in X3F images.
 *
 * Copyright 2015 - Roland and Erik Karlsson
 * BSD-style - see doc/copyright.txt
 *
 */
#ifndef X3F_SPATIAL_GAIN_H
#define X3F_SPATIAL_GAIN_H

#include "x3f_io.h"

namespace x3ftools{

/* Quattro HP: R, G, B0, B1, B2, B3 */
constexpr uint32_t MAXCORR = 6; 


typedef struct{
    double weight, weight_x, weight_y;
    uint32_t *gain;
    double mingain, delta;
}x3f_spatial_gain_corr_merrill_t;


typedef struct{
    double *gain;  /* final gain (interpolated if necessary) */
    int malloc;    /* 1 if gain is allocated with malloc() */
    int rows, cols;
    int rowoff, coloff, rowpitch, colpitch;
    int chan, channels;
    x3f_spatial_gain_corr_merrill_t mgain[4];  /* raw Merrill-type gains */
    int mgain_num;
}x3f_spatial_gain_corr_t;


extern int x3f_get_merrill_type_spatial_gain(x3f_t *x3f, int hp_flag,
					     x3f_spatial_gain_corr_t *corr, bool useLV);

extern int x3f_get_interp_merrill_type_spatial_gain(x3f_t *x3f, int hp_flag,
						    x3f_spatial_gain_corr_t *corr, bool useLV);

extern int x3f_get_classic_spatial_gain(x3f_t *x3f, const char *wb,
					x3f_spatial_gain_corr_t *corr);

extern int x3f_get_spatial_gain(x3f_t *x3f, const char *wb,
				x3f_spatial_gain_corr_t *corr);

extern void x3f_cleanup_spatial_gain(x3f_spatial_gain_corr_t *corr,
				     int corr_num);

extern double x3f_calc_spatial_gain(x3f_spatial_gain_corr_t *corr, int corr_num,
				    int row, int col, int chan, int rows, int cols);


}// namespace x3ftools
#endif
