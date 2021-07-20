/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2021 Gilles Coissac <info@gillescoissac.fr>
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>

#include "settings.h"
#include "camconst.h"
#include "rawimage.h"
#include "foveon.h"

#include "x3f_tools/x3f_io.h"
#include "x3f_tools/x3f_meta.h"
#include "x3f_tools/x3f_image.h"
#include "x3f_tools/x3f_matrix.h"
#include "x3f_tools/x3f_process.h"
#include "x3f_tools/x3f_spatial_gain.h"


namespace X3F=x3ftools;
namespace rtengine{

extern const Settings *settings;


FoveonHelper::FoveonHelper(RawImage const *ri):raw_image(ri),
											   x3f(X3F::x3f_new_from_file(ri->get_file())),
											   data_loaded(false),
											   black({-1,-1,-1}),
											   white({-1,-1,-1}){
    debug({"FoveonHelper for raw ", raw_image->get_filename()});
}


FoveonHelper::~FoveonHelper(){
    if(x3f){
        X3F::x3f_delete(x3f);
        x3f = nullptr;
    }
}


void FoveonHelper::debug(std::vector<std::string> mess) const{
    if (settings->verbose){
        for (size_t i=0; i<mess.size();i++)
            std::cout << mess[i];
        std::cout << std::endl;
    }
}


bool FoveonHelper::load_data(){
	if (!x3f){
		debug({"FoveonHelper: x3f_tools could not read raw ", raw_image->get_filename()});
		return data_loaded;
	}

	if (data_loaded){
		return data_loaded;
	} else {
		debug({"FoveonHelper: x3f_tools loading CAMF Entry..."});
		if (X3F::X3F_OK!=(x3f_load_data(x3f, X3F::x3f_get_camf(x3f)))){
			debug({"x3f_tools could not load CAMF."});
			X3F::x3f_delete(x3f);
			x3f = nullptr;
			return data_loaded;
		} else {
			debug({"FoveonHelper: x3f_tools loading PROP meta data..."});
			X3F::x3f_directory_entry_t *DE_Meta = X3F::x3f_get_prop(x3f);
			if (DE_Meta!=nullptr){ // Not for Quattro
				if (X3F::X3F_OK!=(X3F::x3f_load_data(x3f, DE_Meta))){
					debug({"FoveonHelper: x3f_tools could not load PROP."});
					X3F::x3f_delete(x3f);
					x3f = nullptr;
					return data_loaded;
				}
			}

			debug({"FoveonHelper: x3f_tools loading RAW data..."});
			X3F::x3f_directory_entry_t *DE_Raw = X3F::x3f_get_raw(x3f);
			if (X3F::X3F_OK!=(X3F::x3f_load_data(x3f, DE_Raw))){
				debug({"FoveonHelper: x3f_tools could not load RAW."});
				X3F::x3f_delete(x3f);
				x3f = nullptr;
				return data_loaded;
			}
			data_loaded = true;
		}
	}
	return data_loaded;
}


bool FoveonHelper::get_BlackLevels(int* black_c4){
	if (get_Levels()){
		for (size_t i=0; i<3; i++)
		    black_c4[i] = black[i];
		black_c4[3] = black[1];
		return true;
	}
	return false;
}


bool FoveonHelper::get_WhiteLevels(int* white_c4){
	//return false;
	if (get_Levels()){
		for (size_t i=0; i<3; i++)
			white_c4[i] = white[i];
		white_c4[3] = white[1];
		return true;
	}
	return false;
}


bool FoveonHelper::get_ccMatrix(double *matrix){
    if (!load_data()){
		debug({"FoveonHelper is not in a valid state."});
		return false;
    } else {
        const char *wb = X3F::x3f_get_wb(x3f);
        if (X3F::x3f_get_raw_to_xyz(x3f, wb, matrix)){
        	debug({"FoveonHelper ccMatrix : [", std::to_string(matrix[0]), ", ",
									  		    std::to_string(matrix[1]), ", ",
											    std::to_string(matrix[2]),
				 "\n                         ", std::to_string(matrix[3]), ", ",
											    std::to_string(matrix[4]), ", ",
											    std::to_string(matrix[5]),
				 "\n                         ", std::to_string(matrix[6]), ", ",
											    std::to_string(matrix[7]), ", ",
											    std::to_string(matrix[8]), "]"});
        	return true;
        } else
        	debug({"FoveonHelper can not obtain ccMatrix."});
    }
    return false;
}


bool FoveonHelper::get_cam_xyz(double matrix[4][3]){
    if (!load_data()){
		debug({"FoveonHelper is not in a valid state."});
		return false;
    } else {
        double ccMatrix[9];
        double inverse[9];
        double d65[9];
        double cam_xyz[9];
        if (get_ccMatrix(ccMatrix)){
        	X3F::x3f_3x3_inverse(ccMatrix, inverse);
        	X3F::x3f_Bradford_D50_to_D65(d65);
        	X3F::x3f_3x3_3x3_mul(inverse, d65, cam_xyz);
        	for(int i=0; i<4; i++)
        		for(int j=0; j<3; j++)
        			matrix[i][j] = cam_xyz[(i==3 ? 1 : i) * 3 + j];
        	debug({"FoveonHelper cam_xyz : [", std::to_string(cam_xyz[0]), ", ",
        									   std::to_string(cam_xyz[1]), ", ",
											   std::to_string(cam_xyz[2]),
				 "\n                        ", std::to_string(cam_xyz[3]), ", ",
											   std::to_string(cam_xyz[4]), ", ",
											   std::to_string(cam_xyz[5]),
				 "\n                        ", std::to_string(cam_xyz[6]), ", ",
											   std::to_string(cam_xyz[7]), ", ",
											   std::to_string(cam_xyz[8]), "]"});
        	return true;
        } else
        	debug({"FoveonHelper can not obtain cam_xyz."});
    }
    return false;
}


bool FoveonHelper::get_asShotNeutral(float *gain){
	double total_gain[3];
	double gain_inv[3];
	if (!load_data()){
		debug({"FoveonHelper is not in a valid state."});
	    return false;
	} else {
		const char *wb = X3F::x3f_get_wb(x3f);
		if(X3F::x3f_get_gain(x3f, wb, total_gain)){
			//X3F::x3f_3x1_invert(total_gain, gain_inv);
			gain[0] = float(total_gain[0]);
			gain[1] = float(total_gain[1]);
			gain[2] = float(total_gain[2]);
			debug({"FoveonHelper asShotNeutral: [", std::to_string(gain[0]), ", ",
					std::to_string(gain[1]), ", ", std::to_string(gain[2]), "]"});
			return true;
		} else {
    		debug({"FoveonHelper can not get asShotNeutral multipliers."});
    	}
	}
	return false;
}


bool FoveonHelper::get_wbGain(float *gain){
    double cam_to_xyz[9], wb_correction[9], _gain[3];
    if (!load_data()){
		debug({"FoveonHelper is not in a valid state."});
		return false;
    } else {
    	const char *wb = X3F::x3f_get_wb(x3f);
		if (X3F::x3f_get_camf_matrix_for_wb(x3f, "WhiteBalanceGains", wb, 3, 0, _gain) ||
			X3F::x3f_get_camf_matrix_for_wb(x3f, "DP1_WhiteBalanceGains", wb, 3, 0, _gain)){
			gain[0] = float(_gain[0]);
			gain[1] = float(_gain[1]);
			gain[2] = float(_gain[2]);
			debug({"FoveonHelper wb gain : [", std::to_string(gain[0]), ", ",
							                   std::to_string(gain[1]), ", ",
											   std::to_string(gain[2]), "]"});
			return true;
		} else if (X3F::x3f_get_camf_matrix_for_wb(x3f, "WhiteBalanceIlluminants", wb, 3, 3, cam_to_xyz) &&
				   X3F::x3f_get_camf_matrix_for_wb(x3f, "WhiteBalanceCorrections", wb, 3, 3, wb_correction)){
			double raw_to_xyz[9], raw_neutral[3];
			X3F::x3f_3x3_3x3_mul(wb_correction, cam_to_xyz, raw_to_xyz);
			X3F::get_raw_neutral(raw_to_xyz, raw_neutral);
			X3F::x3f_3x1_invert(raw_neutral, _gain);
			gain[0] = float(_gain[0]);
			gain[1] = float(_gain[1]);
			gain[2] = float(_gain[2]);
			debug({"FoveonHelper wb gain : [", std::to_string(gain[0]), ", ",
			 	                               std::to_string(gain[1]), ", ",
											   std::to_string(gain[2]), "]"});
			return true;
		} else {
    		debug({"FoveonHelper can not get WhiteBalanceGains."});
    	}
    }
    return false;
}


bool FoveonHelper::get_Levels(){
    X3F::x3f_area16_t image, original_image, expanded, crop, image2;
    double black_level[3], black_dev[3];
    //int quattro = x3f_image_area_qtop(x3f, &qtop);
    int colors_in = 3; //quattro ? 2 : 3;
    uint32_t max_raw[3];
    uint32_t white_levels[3];
    double intermediate_bias;
    uint32_t depth;
    const char *wb;

    debug({"FoveonHelper->get_Levels()"});

    if (black[0]!=-1)
    	// we have already computed levels
    	return true;

    if (!load_data()){
        debug({"FoveonHelper is not in a valid state."});
        return false;
    }

    wb = X3F::x3f_get_wb(x3f);
    debug({"FoveonHelper: White Balance camera setup : ", wb});

    // Image depth
    if(!X3F::x3f_get_camf_unsigned(x3f, "ImageDepth", &depth)){
    	debug({"FoveonHelper: Could not get image sensor depth."});
    	return false;
    }

    // preprocess
    if (!X3F::x3f_image_area(x3f, &image) || image.channels<3){
        debug({"FoveonHelper: Could not get image area."});
        return false;
    }

    // TODO: reactivate quattro case
    /*if (quattro && (qtop.channels<1 ||
        qtop.rows<2*image_crop.rows || qtop.columns<2*image_crop.columns))
        return false;
    */

    if (!X3F::get_black_level(x3f, &image, 1, colors_in, black_level, black_dev)){
        //|| (quattro && !get_black_level(x3f, &qtop, 0, 1, &black_level[2], &black_dev[2]))){
        debug({"FoveonHelper: x3f_tools could not get black level."});
        return false;
    }

    if (!X3F::x3f_get_max_raw(x3f, max_raw) ||
        !X3F::get_intermediate_bias(x3f, wb, black_level, black_dev, &intermediate_bias, depth) ||
        !X3F::get_max_intermediate(x3f, wb, intermediate_bias, white_levels, depth)){
    	debug({"FoveonHelper: x3f_tools could not get camera levels."});
    	return false;
    }
    for (size_t i=0; i<3; i++)
        	black[i] = int(intermediate_bias);
        	//black[i] = int(black_level[i]);
    for (size_t i=0; i<3; i++)
        	white[i] = white_levels[i];

    debug({"FoveonHelper: x3f_tools max raw[", std::to_string(max_raw[0]), ", ",
                							   std::to_string(max_raw[1]), ", ",
											   std::to_string(max_raw[2]), "]"});
    debug({"FoveonHelper: x3f_tools black levels[", std::to_string(black_level[0]), ", ",
            										std::to_string(black_level[1]), ", ",
													std::to_string(black_level[2]), "]"});
    debug({"FoveonHelper: x3f_tools white levels[", std::to_string(white[0]), ", ",
                									std::to_string(white[1]), ", ",
    												std::to_string(white[2]), "]"});
    return true;
} // get_BlackLevels


bool FoveonHelper::get_camf_rect(uint32_t *crect) const{
    CameraConstantsStore* ccs = CameraConstantsStore::getInstance();
    CameraConst *cc = ccs->get(raw_image->get_maker().c_str(), raw_image->get_model().c_str());
    
    if(cc->has_rawCrop(raw_image->get_width(), raw_image->get_height())){
        int lm, tm, w, h;
        cc->get_rawCrop(raw_image->get_width(), raw_image->get_height(), lm, tm, w, h);
        //debug({"cc_crop[", std::to_string(lm), ", ", std::to_string(tm),
        //       ", ", std::to_string(w), ", ", std::to_string(h), "]"});
        crect[0] = tm;
        crect[1] = lm;
        crect[2] = raw_image->get_height() + h;
        crect[3] = raw_image->get_width() + w ; 
    } else {
    	debug({"FoveonHelper: could not get camf_rect."});
    	return false;
    }
    return true;
}
 
 
std::vector<GainMap> FoveonHelper::read_foveon_spatial_gain(){
    debug({"FoveonHelper: reading foveon spatial gain correction."});
    std::vector<GainMap> gain_maps_vector;

    if (!load_data()){
        debug({"FoveonHelper is not in a valid state."});
    } else {
        /* Quattro raws are already corrected for spatial gain. 
         * So it is disabled by default. */
        if (x3f->header.version < X3F::X3F_VERSION_4_0){
            X3F::x3f_spatial_gain_corr_t corr[X3F::MAXCORR];
            int corr_num;
            uint32_t active_area[4];
            double originv, originh, scalev, scaleh;
            
            corr_num = X3F::x3f_get_spatial_gain(x3f, X3F::x3f_get_wb(x3f), corr);
            if (corr_num == 0){
                debug({"FoveonHelper: no spatial gain correction found..."});
                X3F::x3f_delete(x3f);
                return gain_maps_vector;
            }
            debug({"FoveonHelper: x3f_tools found ", std::to_string(corr_num), " spatial gain correction(s)..."});
            
            /* Spatial gain in X3F refers to the entire image, 
             * But correction will be applied after cropping to ActiveArea. */
            if (!get_camf_rect(active_area)){
                X3F::x3f_delete(x3f);
                return gain_maps_vector;
            } 
            debug({"FoveonHelper: camf rect ActiveArea[", std::to_string(active_area[0]), ", ",
                    									  std::to_string(active_area[1]), ", ",
														  std::to_string(active_area[2]), ", ",
														  std::to_string(active_area[3]), "]"});
            
            originv = -(double)active_area[0] / (active_area[2] - active_area[0]); // 0
            originh = -(double)active_area[1] / (active_area[3] - active_area[1]); // 0
            scalev = (double)raw_image->get_height() / (active_area[2] - active_area[0]); // 1.0
            scaleh = (double)raw_image->get_width() / (active_area[3] - active_area[1]); // 1.0
            
            /* Compute GainMaps */            
            for (size_t i=0; i<(size_t)corr_num; i++){
                GainMap gain_map;
                X3F::x3f_spatial_gain_corr_t *c = &corr[i];
                gain_map.top = (uint32_t) c->rowoff;
                gain_map.left = (uint32_t) c->coloff;
                gain_map.bottom = (uint32_t) active_area[2] - active_area[0];// this->get_height();
                gain_map.right = (uint32_t) active_area[3] - active_area[1];// this->get_width();
                gain_map.plane = (uint32_t) c->chan;
                gain_map.planes = (uint32_t) c->channels;
                gain_map.row_pitch = (uint32_t) c->rowpitch;
                gain_map.col_pitch = (uint32_t) c->colpitch;
                gain_map.map_points_v = (uint32_t) c->rows;
                gain_map.map_points_h = (uint32_t) c->cols;
                gain_map.map_spacing_v = (double) (scalev/(c->rows-1));
                gain_map.map_spacing_h = (double) (scaleh/(c->cols-1));
                gain_map.map_origin_v = originv;
                gain_map.map_origin_h = originh;
                gain_map.map_planes = (uint32_t) c->channels;
                
                size_t n = gain_map.map_points_v * gain_map.map_points_h * gain_map.map_planes;
                gain_map.map_gain.reserve(n);
                for (size_t j = 0; j < n; ++j){
                    gain_map.map_gain.push_back((float)c->gain[j]);
                }
                gain_maps_vector.push_back(gain_map);
            }
            X3F::x3f_cleanup_spatial_gain(corr, corr_num);
        } else {  // Quattro files
            debug({"FoveonHelper: Quattro files ares already corrected for spatial gain."});
        }
    }
    return gain_maps_vector;
} // read_foveon_spatial_gain


} // rtengine
