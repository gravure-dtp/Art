/* X3F_META.H
 *
 * Library for accessing the CAMF and PROP meta data in the X3F file.
 *
 * Copyright 2015 - Roland and Erik Karlsson
 * BSD-style - see doc/copyright.txt
 *
 */
#ifndef X3F_META_H
#define X3F_META_H


#include "x3f_io.h"

namespace x3ftools{
    

extern int x3f_get_camf_text(x3f_t *x3f, const char *name, char **text);
extern int x3f_get_camf_matrix_var(x3f_t *x3f, const char *name,
				   int *dim0, int *dim1, int *dim2, 
                   matrix_type_t mtype, void **matrix);
extern int x3f_get_camf_matrix(x3f_t *x3f, const char *name,
			       int dim0, int dim1, int dim2, matrix_type_t mtype,
			       void *matrix);
extern int x3f_get_camf_float(x3f_t *x3f, const char *name,  double *val);
extern int x3f_get_camf_float_vector(x3f_t *x3f, const char *name,  double *val);
extern int x3f_get_camf_unsigned(x3f_t *x3f, const char *name,  uint32_t *val);
extern int x3f_get_camf_signed(x3f_t *x3f, const char *name,  int32_t *val);
extern int x3f_get_camf_signed_vector(x3f_t *x3f, const char *name,  int32_t *val);
extern int x3f_get_camf_property_list(x3f_t *x3f, const char *list,
				      char ***names, char ***values,
				      uint32_t *num);
extern int x3f_get_camf_property(x3f_t *x3f, const char *list,
                      const char *name, char **value);
extern int x3f_get_prop_entry(x3f_t *x3f, const char *name, char **value);
extern const char *x3f_get_wb(x3f_t *x3f);
extern int x3f_get_camf_matrix_for_wb(x3f_t *x3f, const char *list, const char *wb, 
                      int dim0, int dim1, double *matrix);
extern int x3f_get_max_raw(x3f_t *x3f, uint32_t *max_raw);

}// namespace x3ftools
#endif
