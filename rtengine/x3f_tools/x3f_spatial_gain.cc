/* X3F_SPATIAL_GAIN.C
 *
 * Library for adjusting for spatial gain in X3F images.
 *
 * Copyright 2015 - Roland and Erik Karlsson
 * BSD-style - see doc/copyright.txt
 *
 */
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wnarrowing"
#endif

#include "x3f_spatial_gain.h"
#include "x3f_meta.h"
#include "x3f_printf.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace x3ftools{

static double lens_position(double focal_length, double object_distance){
    double lp;
    lp = 1.0 / (1.0 / focal_length - 1.0 / object_distance);
    return lp;
}


static double get_focal_length(x3f_t *x3f){
    char *flength;
    double focal_length;

    if (x3f_get_prop_entry(x3f, "FLENGTH", &flength))
        focal_length = atof(flength);
    else {
        focal_length = 30.0;
        x3f_printf(WARN, "Could not get focal length, assuming %g mm\n", focal_length);
    }
    return focal_length;
}


static double get_object_distance(x3f_t *x3f){
    double object_distance;

    if (x3f_get_camf_float(x3f, "ObjectDistance", &object_distance))
        object_distance *= 10.0;	/* Convert cm to mm */
    else {
        object_distance = INFINITY;
        x3f_printf(WARN, "Could not get object distance, assuming %g mm\n", object_distance);
    }
    return object_distance;
}


/* Get the minimum object distance (MOD) */
static double get_MOD(x3f_t *x3f){
    int32_t lens_information;
    double mod;

    if (!x3f_get_camf_signed(x3f, "LensInformation", &lens_information))
        lens_information = -1;

    /* TODO: is there any better way of obtaining MOD? */
    switch (lens_information){
        case 1003:			/* DP1 Merrill */
            mod = 200.0;
            break;
        case 1004:			/* DP2 Merrill */
            mod = 280.0;
            break;
        case 1005:			/* DP3 Merrill */
            mod = 226.0;
            break;
        default:
            mod = 280.0;
            x3f_printf(WARN, "Could not get MOD, assuming %g mm\n", mod);
    }
    return mod;
}


static int get_merrill_type_gains_table(x3f_t *x3f, const char *name, const char *chan,
					uint32_t **mgain, int *rows, int *cols, double *mingain, double *delta, bool useLV){
    char table[32];
    int rows_tmp, cols_tmp;
    
    if(!useLV){
    	char *val;
		sprintf(table, "GainsTable%s", chan);

		if (!x3f_get_camf_property(x3f, name, table, &val) ||
			!x3f_get_camf_matrix_var(x3f, val, &rows_tmp, &cols_tmp, nullptr, M_UINT, (void **)mgain) ||
			(*rows != -1 && *rows != rows_tmp) ||
			(*cols != -1 && *cols != cols_tmp)){
			return 0;
		}

		*rows = rows_tmp;
		*cols = cols_tmp;

		sprintf(table, "MinGains%s", chan);
		if (!x3f_get_camf_property(x3f, name, table, &val))
			return 0;
		*mingain = atof(val);

		sprintf(table, "Delta%s", chan);
		if (!x3f_get_camf_property(x3f, name, table, &val))
			return 0;
		*delta = atof(val);
    } else { // use SpatialGain#_LV table
    	char val[50] = {0};
    	strcpy(val, name);
    	val[11] = chan[0];
    	if (!x3f_get_camf_matrix_var(x3f, val, &rows_tmp, &cols_tmp, nullptr, M_UINT, (void **)mgain) ||
    		(*rows != -1 && *rows != rows_tmp) ||
    		(*cols != -1 && *cols != cols_tmp)){
			return 0;
		}
		*rows = rows_tmp;
		*cols = cols_tmp;
		*mingain = 0;
    	*delta = 1.0;
    }

    return 1;
}


typedef struct merrill_spatial_gain_s{
    char *name;
    double x,y;
    struct merrill_spatial_gain_s *next;
}merrill_spatial_gain_t;


int x3f_get_merrill_type_spatial_gain(x3f_t *x3f, int hp_flag, x3f_spatial_gain_corr_t *corr, bool useLV){
    char **include_blocks, **include_blocks_val;
    uint32_t include_blocks_num;
    double capture_aperture;
    double *spatial_gain_fstop;
    int num_fstop;
    int corr_num = 3;

    merrill_spatial_gain_t *gains = nullptr, *g;
    merrill_spatial_gain_t *q_closest[4] = {nullptr, nullptr, nullptr, nullptr};
    double q_closest_dx[4] = {INFINITY, -INFINITY, -INFINITY, INFINITY};
    double q_closest_dy[4] = {INFINITY, INFINITY, -INFINITY, -INFINITY};
    double q_closest_d2[4] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double q_weight_x[4], q_weight_y[4], q_weight[4];
 
    double x,y;
    int i;

    if (!x3f_get_camf_float(x3f, "CaptureAperture", &capture_aperture)) 
        return 0;

    if (!x3f_get_camf_property_list(x3f, "IncludeBlocks",
				  &include_blocks, &include_blocks_val,
				  &include_blocks_num))
        return 0;

    if (hp_flag){  /* Quattro HP */
        if (!x3f_get_camf_matrix_var(x3f, "SpatialGainHP_Fstop",
				 &num_fstop, nullptr, nullptr, M_FLOAT, (void **)&spatial_gain_fstop))
            return 0;
        corr_num = 6;  /* R, G, B0, B1, B2, B3 */

        for (i=0; i<include_blocks_num; i++){
            char **names, **values;
            uint32_t num, aperture_index;
            char dummy;

            if (sscanf(include_blocks[i], "SpatialGainHPProps_%d%c",
                &aperture_index, &dummy)==1 &&
                x3f_get_camf_property_list(x3f, include_blocks[i], &names, &values, &num) &&
                aperture_index>=0 && aperture_index<num_fstop){
                    g = static_cast<merrill_spatial_gain_t*>(alloca(sizeof(merrill_spatial_gain_t)));
                    g->name = include_blocks[i];
                    //g->x = 1.0 / spatial_gain_fstop[aperture_index];
                    g->x = get_focal_length(x3f) / spatial_gain_fstop[aperture_index];
                    g->y = 0.0;  /* unused */
                    g->next = gains;
                    gains = g;
            }
        }
        //x = 1.0 / capture_aperture;
        x = get_focal_length(x3f) / capture_aperture;
        y = 0.0;  /* unused */
    } else {  /* Merrill and Quattro */
        if (x3f_get_camf_matrix_var(x3f, "SpatialGain_Fstop",
            &num_fstop, nullptr, nullptr, M_FLOAT, (void **)&spatial_gain_fstop)){

        	if (useLV){ // use SpatialGain#_LV table
        		for (i=0; i<include_blocks_num; i++){
					char **names, **values;
					uint32_t num, aperture, focus_dist_index;
					char dummy;

					if (sscanf(include_blocks[i], "SpatialGainR_LV_%2d_%2d%c",
							&aperture, &focus_dist_index, &dummy)==2){
						double fstop, lenspos;
						fstop = double(aperture) / 10.0;

						if (focus_dist_index==7)
							lenspos = lens_position(get_focal_length(x3f), INFINITY);
						else if (focus_dist_index==0)
							lenspos = lens_position(get_focal_length(x3f), get_MOD(x3f));
						else continue;

						g = static_cast<merrill_spatial_gain_t*>(alloca(sizeof(merrill_spatial_gain_t)));
						g->name = include_blocks[i];
						g->x = get_focal_length(x3f) / fstop;
						g->y = lenspos;
						g->next = gains;
						gains = g;
						x3f_printf(DEBUG, "%s -> x(f/%f) = %fmm, y(lens_position) = %fmm\n",
								include_blocks[i], fstop, g->x, lenspos);
					}
        		}
        	} else {
				for (i=0; i<include_blocks_num; i++){
					char **names, **values;
					uint32_t num, aperture_index;
					char focus_distance[4], dummy;

					if (sscanf(include_blocks[i], "SpatialGainsProps_%d_%3s%c",
							&aperture_index, focus_distance, &dummy)==2 &&
							x3f_get_camf_property_list(x3f, include_blocks[i], &names, &values, &num) &&
							aperture_index>=0 && aperture_index<num_fstop){
						double lenspos;

						if (!strcmp(focus_distance, "INF"))
							lenspos = lens_position(get_focal_length(x3f), INFINITY);
						else if (!strcmp(focus_distance, "MOD"))
							lenspos = lens_position(get_focal_length(x3f), get_MOD(x3f));
						else continue;

						g = static_cast<merrill_spatial_gain_t*>(alloca(sizeof(merrill_spatial_gain_t)));
						g->name = include_blocks[i];
						//g->x = 1.0 / spatial_gain_fstop[aperture_index];
						g->x = get_focal_length(x3f) / spatial_gain_fstop[aperture_index];
						g->y = lenspos;
						g->next = gains;
						gains = g;
					}
				}
        	}
        } else {
            for (i=0; i<include_blocks_num; i++){
                char **names, **values;
                uint32_t num;
                double aperture, lenspos;
                char dummy;

                if (sscanf(include_blocks[i], "SpatialGainsProps_%lf_%lf%c",
                        &aperture, &lenspos, &dummy)==2 &&
                        x3f_get_camf_property_list(x3f, include_blocks[i], 
                        &names, &values, &num)){
                    g = static_cast<merrill_spatial_gain_t*>(alloca(sizeof(merrill_spatial_gain_t)));
                    g->name = include_blocks[i];
                    //g->x = 1.0 / aperture;
                    g->x = get_focal_length(x3f) / aperture;
                    g->y = lenspos;
                    g->next = gains;    
                    gains = g;
                }
            }
        }
        //x = 1.0 / capture_aperture;
        x = get_focal_length(x3f) / capture_aperture;
        y = lens_position(get_focal_length(x3f), get_object_distance(x3f));
    } /* End Merrill and Quattro */

    /* TODO: doing bilinear interpolation with respect 
     * to 1/aperture and lens position. Is that correct? */
    for (g=gains; g; g=g->next){
        double dx = g->x - x;
        double dy = g->y - y;
        double d2 = dx * dx + dy * dy;
        int q;

        if     (dx>0.0 && dy>0.0) q = 0;
        else if(dx>0.0)           q = 3;
        else if(dy>0.0)           q = 1;
        else                      q = 2;

        if (d2<q_closest_d2[q]){
            q_closest[q] = g;
            q_closest_dx[q] = dx;
            q_closest_dy[q] = dy;
            q_closest_d2[q] = d2;
        }
    }

    /* TODO: bilinear interpolation requires that the points be laid out
       in a more or less rectilinear grid. Is that assumption good
       enough? It appears to be valid for the current cameras. */
    q_weight_x[0] = q_closest_dx[1] / (q_closest_dx[1] - q_closest_dx[0]);
    q_weight_x[1] = q_closest_dx[0] / (q_closest_dx[0] - q_closest_dx[1]);
    q_weight_x[2] = q_closest_dx[3] / (q_closest_dx[3] - q_closest_dx[2]);
    q_weight_x[3] = q_closest_dx[2] / (q_closest_dx[2] - q_closest_dx[3]);

    q_weight_y[0] = q_closest_dy[3] / (q_closest_dy[3] - q_closest_dy[0]);
    q_weight_y[1] = q_closest_dy[2] / (q_closest_dy[2] - q_closest_dy[1]);
    q_weight_y[2] = q_closest_dy[1] / (q_closest_dy[1] - q_closest_dy[2]);
    q_weight_y[3] = q_closest_dy[0] / (q_closest_dy[0] - q_closest_dy[3]);

    /* DEBUG */
    x3f_printf(DEBUG, "\n     Capture shot -> x(f/%f) = %fmm, y(lens_position) = %fmm\n", 
                capture_aperture,x, y);
    for (i=0; i<4; i++){
        if (q_closest[i]){
            x3f_printf(DEBUG, "q[%d] %s -> x = %fmm, y = %fmm\n",
                i, q_closest[i]->name, q_closest[i]->x, q_closest[i]->y);
        }else {
            x3f_printf(DEBUG, "q[%d] = nullptr\n", i);
        }
        x3f_printf(DEBUG, "q[%d] dx = %f, dy = %f, d2 = %f, wx = %f, wy = %f\n",
                i, q_closest_dx[i], q_closest_dy[i], q_closest_d2[i],
                q_weight_x[i], q_weight_y[i]);
    }

    for (i=0; i<4; i++){
        if (std::isnan(q_weight_x[i])) 
            q_weight_x[i] = 1.0;
        if (std::isnan(q_weight_y[i])) 
            q_weight_y[i] = 1.0;
        q_weight[i] = q_weight_x[i] * q_weight_y[i];
    }

    /* DEBUG */
    double v = 0;
    for (i=0; i<4; i++){
        v += q_weight[i];
        x3f_printf(DEBUG, "q[%d] %s -> w = %f\n",
	       i, q_closest[i] ? q_closest[i]->name : "nullptr", q_weight[i]);
    }
    x3f_printf(DEBUG, "     Weight sum verification = %f\n\n", v); 
    
    for (i=0; i<corr_num; i++){
        x3f_spatial_gain_corr_t *c = &corr[i];
        c->gain = nullptr;	 /* No interploated gains present */
        c->malloc = 0;
        c->rows = c->cols = -1;
        c->rowoff = c->coloff = 0;
        c->rowpitch = c->colpitch = 1;
        c->chan = i;
        c->channels = 1;
        c->mgain_num = 0;
    }

    if (hp_flag){
        corr[2].rowoff = 0;
        corr[2].coloff = 0;
        corr[2].rowpitch = corr[2].colpitch = 2;
        corr[2].chan = 2;

        corr[3].rowoff = 0;
        corr[3].coloff = 1;
        corr[3].rowpitch = corr[3].colpitch = 2;
        corr[3].chan = 2;

        corr[4].rowoff = 1;
        corr[4].coloff = 0;
        corr[4].rowpitch = corr[4].colpitch = 2;
        corr[4].chan = 2;

        corr[5].rowoff = 1;
        corr[5].coloff = 1;
        corr[5].rowpitch = corr[5].colpitch = 2;
        corr[5].chan = 2;
    }

    for (i=0; i<4; i++){
        char *name;
        const char *channels_normal[3] = {"R", "G", "B"};
        const char *channels_hp[6] = {"R", "G", "B0", "B1", "B2", "B3"};
        const char **channels = hp_flag ? channels_hp : channels_normal;
        int j;

        if (!q_closest[i]) 
            continue;
        name = q_closest[i]->name;

        for (j=0; j<corr_num; j++){
            x3f_spatial_gain_corr_t *c = &corr[j];
            x3f_spatial_gain_corr_merrill_t *m = &c->mgain[c->mgain_num++];
            m->weight = q_weight[i];
            m->weight_x = q_weight_x[i];
            m->weight_y = q_weight_y[i];

            if (!get_merrill_type_gains_table(x3f, name, channels[j],
					&m->gain, &c->rows, &c->cols, &m->mingain, &m->delta, useLV))
                return 0;
        }
    }

    for (i=0; i<corr_num; i++)
        if (corr[i].mgain_num == 0) 
            return 0;

    return corr_num;
} // x3f_get_merrill_type_spatial_gain


int x3f_get_interp_merrill_type_spatial_gain(x3f_t *x3f, int hp_flag,
					     x3f_spatial_gain_corr_t *corr, bool useLV){
    int corr_num = x3f_get_merrill_type_spatial_gain(x3f, hp_flag, corr, useLV);
    int i;
    
    for (i=0; i<corr_num; i++){
        x3f_spatial_gain_corr_t *c = &corr[i];
        int j, num = c->rows * c->cols * c->channels;
        c->gain = static_cast<double*>(malloc(num * sizeof(double)));
        c->malloc = 1;

        for (j=0; j<num; j++){
        	c->gain[j] = 0.0;
        	if (!useLV){
				for (int g=0; g<c->mgain_num; g++){
					x3f_spatial_gain_corr_merrill_t *m = &c->mgain[g];
					c->gain[j] += m->weight * (m->mingain + m->delta * m->gain[j]);
					//c->gain[j] += m->weight * (m->mingain + (m->gain[j] / 4096.0 / m->delta / 100));
				}
				//x3f_printf(DEBUG, "%f \n", c->gain[j]);
        	} else {
        		double dg;
        		for (int g=0; g<c->mgain_num; g++){
					x3f_spatial_gain_corr_merrill_t *m = &c->mgain[g];
					dg = 1.0 / (double(m->gain[j]) / 1024);
					c->gain[j] += m->weight * dg;
        		}
        		//x3f_printf(DEBUG, "%f \n", c->gain[j]);
        	}

        	/*
            x3f_spatial_gain_corr_merrill_t *m0 = &c->mgain[0];
            x3f_spatial_gain_corr_merrill_t *m1 = &c->mgain[1];
            x3f_spatial_gain_corr_merrill_t *m2 = &c->mgain[2];
            x3f_spatial_gain_corr_merrill_t *m3 = &c->mgain[3];

            c->gain[j] = ((m0->weight_x * (m0->mingain + m0->delta * m0->gain[j]) + 
                           m1->weight_x * (m1->mingain + m1->delta * m1->gain[j])) * m0->weight_y) + 
                         ((m2->weight_x * (m2->mingain + m2->delta * m2->gain[j]) + 
                           m3->weight_x * (m3->mingain + m3->delta * m3->gain[j])) * m2->weight_y); 
            */
        }
    }
    
    /*double CCGR[3] = {0.97286, 0.9546, 0.9758};
    x3f_spatial_gain_corr_t *cr = &corr[0];
    x3f_spatial_gain_corr_t *cg = &corr[1];
    x3f_spatial_gain_corr_t *cb = &corr[2];
    int j, num = cr->rows * cr->cols * cr->channels;
    for (j=0; j<num; j++){
        cg->gain[j] = cr->gain[j] * CCGR[0] + cg->gain[j] * CCGR[1] + cb->gain[j] * CCGR[2]; 
    }*/
    return corr_num;
}


int x3f_get_classic_spatial_gain(x3f_t *x3f, const char *wb, x3f_spatial_gain_corr_t *corr){
    char *gain_name;
    if ((x3f_get_camf_property(x3f, "SpatialGainTables", wb, &gain_name) &&
         x3f_get_camf_matrix_var(x3f, gain_name, &corr->rows, &corr->cols, &corr->channels,
			       M_FLOAT, (void **)&corr->gain)) ||
         x3f_get_camf_matrix_var(x3f, "SpatialGain", &corr->rows, &corr->cols, &corr->channels,
			      M_FLOAT, (void **)&corr->gain)){
        corr->malloc = 0;
        corr->rowoff = corr->coloff = 0;
        corr->rowpitch = corr->colpitch = 1;
        corr->chan = 0;
        corr->mgain_num = 0;
        return 1;
    }
    return 0;
}


int x3f_get_spatial_gain(x3f_t *x3f, const char *wb, x3f_spatial_gain_corr_t *corr){
    int corr_num = 0;
    bool useLv = false; //FIXME: testing LV table for merrill
    corr_num += x3f_get_interp_merrill_type_spatial_gain(x3f, 0, &corr[corr_num], useLv);
    if (corr_num == 0)
        corr_num += x3f_get_classic_spatial_gain(x3f, wb, &corr[corr_num]);
    return corr_num;
}


void x3f_cleanup_spatial_gain(x3f_spatial_gain_corr_t *corr, int corr_num){
    int i;
    for (i=0; i<corr_num; i++)
        if (corr[i].malloc) free(corr[i].gain);
}


double x3f_calc_spatial_gain(x3f_spatial_gain_corr_t *corr, int corr_num,
			     int row, int col, int chan, int rows, int cols){
    double gain = 1.0;
    double rrel = (double)row / rows;
    double crel = (double)col / cols;
    int i;

    for (i=0; i<corr_num; i++){
        x3f_spatial_gain_corr_t *c = &corr[i];
        double rc, cc;
        int ri, ci;
        double rf, cf;
        double *r1, *r2;
        int co1, co2;
        double gr1, gr2;
        int ch = chan - c->chan;

        if (ch<0 || ch >= c->channels) continue;
        if (row%c->rowpitch!=c->rowoff) continue;
        if (col%c->colpitch!=c->coloff) continue;

        rc = rrel * (c->rows - 1);
        ri = (int)floor(rc);
        rf = rc - ri;

        cc = crel * (c->cols - 1);
        ci = (int)floor(cc);
        cf = cc - ci;

        if (ri<0)
            r1 = r2 = &c->gain[0];
        else if (ri>=c->rows - 1)
            r1 = r2 = &c->gain[(c->rows -1 ) * c->cols * c->channels];
        else {
            r1 = &c->gain[ri * c->cols * c->channels];
            r2 = &c->gain[(ri + 1) * c->cols * c->channels];
        }

        if (ci<0)
            co1 = co2 = ch;
        if (ci>=c->cols - 1)
            co1 = co2 = (c->cols - 1) * c->channels + ch;
        else {
            co1 = ci * c->channels + ch;
            co2 = (ci + 1) * c->channels + ch;
        }

        gr1 = r1[co1] + cf * (r1[co2] - r1[co1]);
        gr2 = r2[co1] + cf * (r2[co2] - r2[co1]);

        gain *= gr1 + rf * (gr2 - gr1);
    }
    return gain;
}


}// namespace x3ftools
