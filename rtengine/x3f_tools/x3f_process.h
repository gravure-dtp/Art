/* X3F_PROCESS.H
 *
 * Library for processing X3F data.
 *
 * Copyright 2015 - Roland and Erik Karlsson
 * BSD-style - see doc/copyright.txt
 *
 */

#ifndef X3F_PROCESS_H
#define X3F_PROCESS_H

#include "x3f_io.h"

namespace x3ftools{

typedef enum x3f_color_encoding_e{
    NONE=0,		   /* Preprocessed but unconverted data */
    SRGB=1,		   /* Preproccesed and convered to sRGB */
    ARGB=2,		   /* Preproccesed and convered to Adobee RGB */
    PPRGB=3,	   /* Preproccesed and convered to ProPhoto RGB */
    UNPROCESSED=4, /* RAW data without any preprocessing */
    QTOP=5,		   /* Quattro top layer without any preprocessing */
}x3f_color_encoding_t;


typedef struct{
    double black[3];
    uint32_t white[3];
} x3f_image_levels_t;

extern int x3f_get_gain(x3f_t *x3f, const char *wb, double *gain);
extern int x3f_get_bmt_to_xyz(x3f_t *x3f, const char *wb, double *bmt_to_xyz);
extern int x3f_get_raw_to_xyz(x3f_t *x3f, const char *wb, double *raw_to_xyz);

extern int x3f_get_image(x3f_t *x3f,
			 x3f_area16_t *image,
			 x3f_image_levels_t *ilevels,
			 x3f_color_encoding_t encoding,
			 int crop,
			 int fix_bad,
			 int denoise,
			 int apply_sgain,
			 const char *wb);

extern int x3f_get_preview(x3f_t *x3f,
			   x3f_area16_t *image,
			   x3f_image_levels_t *ilevels,
			   x3f_color_encoding_t encoding,
			   int apply_sgain,
			   const char *wb,
			   uint32_t max_width,
			   x3f_area8_t *preview);

extern int get_intermediate_bias(x3f_t *x3f, const char *wb, double *black_level,
                                 double *black_dev, double *intermediate_bias, int intermediate_depth);
extern int get_max_intermediate(x3f_t *x3f, const char *wb,
        double intermediate_bias, uint32_t *max_intermediate, int intermediate_depth);
extern int get_black_level(x3f_t *x3f, x3f_area16_t *image, int rescale, 
                    int colors, double *black_level, double *black_dev);
extern void get_raw_neutral(double *raw_to_xyz, double *raw_neutral);

}// namespace x3ftools
#endif

