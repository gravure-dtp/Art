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
											   white({-1,-1,-1}),
											   scale({1.0, 1.0, 1.0}),
											   cam_wb(nullptr),
											   is_quattro(false),
											   crop_area({0,0,0,0}){
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
	if (get_raw_property()){
		for (size_t i=0; i<3; i++)
		    black_c4[i] = int(black[i]);
		black_c4[3] = int(black[1]);
		return true;
	}
	return false;
}


bool FoveonHelper::get_BlackDev(int* black){
	if (get_raw_property()){
		for (size_t i=0; i<3; i++)
		    black[i] = int(black_dev[i]);
		black[3] = int(black_dev[1]);
		return true;
	}
	return false;
}


bool FoveonHelper::get_BlackBias(int* black){
	if (get_raw_property()){
		for (size_t i=0; i<3; i++)
		    black[i] = int(black_bias[i]);
		black[3] = int(black_bias[1]);
		return true;
	}
	return false;
}


bool FoveonHelper::get_WhiteLevels(int* white_c4){
	if (get_raw_property()){
		for (size_t i=0; i<3; i++)
			white_c4[i] = white[i];
		white_c4[3] = white[1];
		return true;
	}
	return false;
}


bool FoveonHelper::get_raw_to_xyz(const char *wb, double *matrix){
	double bmt_to_xyz[9], gain[9], gain_mat[9];

	if (get_gain(wb, gain) && X3F::x3f_get_bmt_to_xyz(x3f, wb, bmt_to_xyz)){
		X3F::x3f_3x3_diag(gain, gain_mat);
		X3F::x3f_3x3_3x3_mul(bmt_to_xyz, gain_mat, matrix);
	} else
    	debug({"FoveonHelper can not obtain raw_to_xyz."});
	return false;
}


bool FoveonHelper::get_ccMatrix(double *matrix){
    if (!get_raw_property()){
		debug({"FoveonHelper is not in a valid state."});
		return false;
    } else {
        /*double bmt_xyz[9];
        double neutral_mat[9];
        double gain[3];
        double neutral[3];
        if (X3F::x3f_get_bmt_to_xyz(x3f, cam_wb, bmt_xyz) && get_asShotNeutral(gain)){
        	X3F::x3f_3x1_invert(gain, neutral);
        	X3F::x3f_3x3_diag(neutral, neutral_mat);
        	X3F::x3f_3x3_3x3_mul(bmt_xyz, neutral_mat, matrix);
        */
        if (X3F::x3f_get_raw_to_xyz(x3f, cam_wb, matrix)){
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
} // get_ccMatrix


bool FoveonHelper::get_cam_xyz(double matrix[4][3]){
    if (!get_raw_property()){
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
} // get_cam_xyz


bool FoveonHelper::get_gain(const char *wb, double *gain){
	double gain_fact[3];

	if (get_wbgain(wb, gain)){
		get_sensor_gain(gain_fact);
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

		debug({"FoveonHelper color gain corrections: [", std::to_string(gain[0]), ", ",
							           	   	   	   	   	 std::to_string(gain[1]), ", ",
														 std::to_string(gain[2]), "]"});
		return true;
	} else {
	    debug({"FoveonHelper can not get color gain corrections."});
	}
	return false;
} // get_gain


bool FoveonHelper::get_gain(double *gain){
	return (get_gain(cam_wb, gain));
}


bool FoveonHelper::get_wbgain(const char *wb, double *gain){
    double cam_to_xyz[9], wb_correction[9], _gain[3];
    if (!get_raw_property()){
		debug({"FoveonHelper is not in a valid state."});
		return false;
    } else {
		if (X3F::x3f_get_camf_matrix_for_wb(x3f, "WhiteBalanceGains", wb, 3, 0, _gain) ||
			X3F::x3f_get_camf_matrix_for_wb(x3f, "DP1_WhiteBalanceGains", wb, 3, 0, _gain)){
			gain[0] = (_gain[0]);
			gain[1] = (_gain[1]);
			gain[2] = (_gain[2]);
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
			gain[0] = (_gain[0]);
			gain[1] = (_gain[1]);
			gain[2] = (_gain[2]);
			debug({"FoveonHelper wb gain : [", std::to_string(gain[0]), ", ",
			 	                               std::to_string(gain[1]), ", ",
											   std::to_string(gain[2]), "]"});
			return true;
		} else {
    		debug({"FoveonHelper can not get WhiteBalanceGains."});
    	}
    }
    return false;
} // get_wbgain


void FoveonHelper::get_sensor_gain(double *gain){
	double gain_fact[3];
	int i_gain_fact[3];
	gain[0] = gain[1] = gain[2] = 1.0;

	if (X3F::x3f_get_camf_float_vector(x3f, "SensorAdjustmentGainFact", gain_fact))
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	if (X3F::x3f_get_camf_float_vector(x3f, "TempGainFact", gain_fact))
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	if (X3F::x3f_get_camf_float_vector(x3f, "FNumberGainFact", gain_fact))
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	/*
	if (X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_BR", gain_fact))
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	if (X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_GR", gain_fact))
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	if (X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_RR", gain_fact))
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	*/
	/*double maxd = 0;
	if (X3F::x3f_get_camf_signed_vector(x3f, "DespAdjust", i_gain_fact)){
		for (int i=0; i<3; i++){
			gain_fact[i] = double(i_gain_fact[i]);
			if (maxd < gain_fact[i])
				maxd = gain_fact[i];
		}
		for (int i=0; i<3; i++)
			gain_fact[i] = 1.0 / (gain_fact[i] / maxd);
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);
	}*/

	debug({"FoveonHelper: sensor gain [", std::to_string(gain[0]), ", ",
		                				  std::to_string(gain[1]), ", ",
										  std::to_string(gain[2]), "]"});
} // get_sensor_gain


void FoveonHelper::get_correct_color_gain(std::array<std::array<float, 3>, 3> &gain_adj){
	double gain_fact[9];
	double gain[9];
	gain[0] = gain[1] = gain[2] = 1.0;

	(X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_BR", gain_fact));
		//X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	(X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_GR", &gain_fact[3]));
		//X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	(X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_RR", &gain_fact[6]));
		//X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);

	for (int c=0; c<3; c++)
		for(int r=0; r<3; r++)
			gain_adj[c][r] = float(gain_fact[r + 3 * c]);
	debug({"FoveonHelper: spatial color correct gain \n[", std::to_string(gain_adj[0][0]), ", ",
	                							         std::to_string(gain_adj[0][1]), ", ",
												         std::to_string(gain_adj[0][2]), "\n ",
														 std::to_string(gain_adj[1][0]), ", ",
														 std::to_string(gain_adj[1][1]), ", ",
														 std::to_string(gain_adj[1][2]), "\n ",
														 std::to_string(gain_adj[2][0]), ", ",
														 std::to_string(gain_adj[2][1]), ", ",
														 std::to_string(gain_adj[2][2]), "\n ",
														 "]"});
} // get_correct_color_gain


void FoveonHelper::get_spatial_gain_adjust(float *adjust){
	double gain_fact[3];
	int i_gain_fact[3];
	double gain[3] = {1.0, 1.0, 1.0}, sgain[3];

	get_sensor_gain(sgain);
//	(X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_BR", gain_fact));
//	X3F::x3f_3x1_comp_mul(gain_fact, sgain, sgain);
//	(X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_GR", gain_fact));
//	X3F::x3f_3x1_comp_mul(gain_fact, sgain, sgain);
//	(X3F::x3f_get_camf_float_vector(x3f, "CorrectColorGain_RR", gain_fact));
//	X3F::x3f_3x1_comp_mul(gain_fact, sgain, sgain);

	double maxd = 0;
	if (get_raw_property()){
		for (int c=0; c<3; c++){
			//gain_fact[c] = black_dev[c] / black[c] ;
			gain_fact[c] = (black_dev[c]) / (black[c]) ;
			if (maxd < gain_fact[c])
				maxd = gain_fact[c];
		}
		for (int c=0; c<3; c++)
			gain_fact[c] = maxd / gain_fact[c];
		X3F::x3f_3x1_comp_mul(gain_fact, gain, gain);
	}

	for (int c=0; c<3; c++) adjust[c] = gain[c];
	debug({"FoveonHelper: spatial gain adujstement [", std::to_string(adjust[0]), ", ",
		                							   std::to_string(adjust[1]), ", ",
													   std::to_string(adjust[2]), "]"});
} // get_spatial_gain_adjust


bool FoveonHelper::get_lin_lut(){
	int dim0=0;
	int dim1=0;
	uint32_t lut_depth;
	uint16_t lin_lut[3][4096];
	int32_t *matrix_RR, *matrix_GR, *matrix_BR;
	double maxy, x, dx, tmp;
	int ix, tx;


	if(X3F::x3f_get_camf_unsigned(x3f, "LinLUTBitDepth", &lut_depth) &&
	   X3F::x3f_get_camf_matrix_var(x3f, "LinLUTs_RR", &dim0, &dim1, nullptr, X3F::M_INT, (void **)&matrix_RR) &&
	   X3F::x3f_get_camf_matrix_var(x3f, "LinLUTs_BR", &dim0, &dim1, nullptr, X3F::M_INT, (void **)&matrix_BR) &&
	   X3F::x3f_get_camf_matrix_var(x3f, "LinLUTs_GR", &dim0, &dim1, nullptr, X3F::M_INT, (void **)&matrix_GR)){

		maxy = (1 << lut_depth) - 1;

		for (int d=0; d<4096; d++){
			x = (d / 4095.0) * (dim1 - 1);
			ix = int(x);
			dx = x - ix;
			//debug({std::to_string(x), ", ", std::to_string(ix), ", ", std::to_string(dx)});
			for (int c=0; c<3; c++){
				tx = ix + (dim1 * c);
				if(x >= dim1 -1)
					tmp = double(matrix_GR[tx]);
				else
					tmp = double(matrix_GR[tx]) * (1.0 - dx) + double(matrix_GR[tx + 1]) * dx;
				lin_lut[c][d] = static_cast<uint16_t>(tmp / maxy * 4095);

				/*if(c==0)
					debug({std::to_string(tmp), "->", std::to_string(lin_lut[c][d]), "  ",
					       std::to_string(matrix_RR[tx]), ", ", std::to_string(matrix_RR[tx+1])});
				*/
			}
		}
		/*
		for (int i=0; i<4096; i++)
			debug({"(", std::to_string(lin_lut[0][i]), ", ",
			            std::to_string(lin_lut[1][i]), ", ",
						std::to_string(lin_lut[2][i]), ")" });
		*/
		return true;
	} else {
    	debug({"FoveonHelper: could not get LinLut Matrix."});
    }
	return false;
} // get_lin_lut


bool FoveonHelper::get_raw_property(){
    X3F::x3f_area16_t full_image, expanded, crop_image, qtop;
    double black_level[3], _black_dev[3];
    uint32_t max_raw[3], white_levels[3];
    double intermediate_bias;
    uint32_t depth;
    int channels;
    bool crop = false;
    uint32_t area[4];
    const char *area_string[2] = {"ActiveImageArea", "KeepImageArea"};
    int area_index = 1;

    if (black[0]!=-1) // we have already computed levels
    	return true;

    if (!load_data()){
        debug({"FoveonHelper is not in a valid state."});
        return false;
    }

    // camera wb setup
    cam_wb = X3F::x3f_get_wb(x3f);
    debug({"FoveonHelper: White Balance camera setup : ", cam_wb});

    // Image depth
    if(!X3F::x3f_get_camf_unsigned(x3f, "ImageDepth", &depth)){
    	debug({"FoveonHelper: Could not get image sensor depth."});
    	return false;
    }

    // getting image_areas
    is_quattro = bool(x3f_image_area_qtop(x3f, &qtop));
    channels = is_quattro ? 2 : 3;
    if (is_quattro){ // TODO: support quattro files
    	debug({"FoveonHelper: Doesn't support quattro camera."});
    	return false;
    }

    if (!X3F::x3f_image_area(x3f, &full_image) || full_image.channels < 3){
        debug({"FoveonHelper: Could not get image area."});
        return false;
    }
//    if (!crop || !X3F::x3f_crop_area_camf(x3f, "ActiveImageArea", &full_image, 1, &crop_image))
//    	crop_image = full_image;

//    if (quattro && (qtop.channels<1 ||
//        qtop.rows<2*image_crop.rows || qtop.columns<2*image_crop.columns))
//        return false;

    if(crop) area_index = 0;
    if (!X3F::x3f_get_camf_rect(x3f, area_string[area_index], &full_image, 1, area)){
    	debug({"FoveonHelper: Could not get camf rect."});
		return false;
    }
    // Translate from Sigma's to Adobe's view on rectangles
	crop_area[0] = area[1];
	crop_area[1] = area[0];
	crop_area[2] = area[3] + 1;
	crop_area[3] = area[2] + 1;

    // Raw Levels
    if (!X3F::get_black_level(x3f, &full_image, 1, channels, black_level, _black_dev)){
        //|| (quattro && !get_black_level(x3f, &qtop, 0, 1, &black_level[2], &black_dev[2]))){
        debug({"FoveonHelper: x3f_tools could not get black level."});
        return false;
    }

    double gain[3], max_intermediate[3], maxgain=0;
    if (X3F::x3f_get_max_raw(x3f, max_raw)){
    	get_sensor_gain(gain);
    	for (int i=0; i<3; i++) if (gain[i] > maxgain) maxgain = gain[i];
    	for (int i=0; i<3; i++)
    		max_intermediate[i] = gain[i] * ((1 << depth) - 1) / maxgain;

    	intermediate_bias = 0.0;
		for (int i=0; i<3; i++){
			double bias = _black_dev[i] * max_intermediate[i] / (max_raw[i] - black_level[i]);
			if (bias > intermediate_bias) intermediate_bias = bias;
		}

		for (int i=0; i<3; i++)
			max_intermediate[i] = gain[i] * (((1 << depth) - 1) - intermediate_bias) / maxgain + intermediate_bias;
    } else{
    	debug({"FoveonHelper: x3f_tools could not get white levels."});
    	return false;
    }

    for (size_t i=0; i<3; i++){
        black[i] = black_level[i];
    	black_bias[i] = intermediate_bias;
    	black_dev[i] = _black_dev[i];
        white[i] = int(max_intermediate[i]);
        scale[i] = (max_intermediate[i] - intermediate_bias) / (max_raw[i] - black_level[i]);
    }

    debug({"FoveonHelper: max raw[", std::to_string(max_raw[0]), ", ",
                					 std::to_string(max_raw[1]), ", ",
									 std::to_string(max_raw[2]), "]"});

    debug({"FoveonHelper: black levels [", std::to_string(black[0]), ", ",
            							  std::to_string(black[1]), ", ",
										  std::to_string(black[2]), "]"});

    debug({"FoveonHelper: black levels deviation [", std::to_string(black_dev[0]), ", ",
                        				   	   	   	 std::to_string(black_dev[1]), ", ",
													 std::to_string(black_dev[2]), "]"});

    debug({"FoveonHelper: black bias [", std::to_string(black_bias[0]), ", ",
                						std::to_string(black_bias[1]), ", ",
    									std::to_string(black_bias[2]), "]"});

    debug({"FoveonHelper: white levels [", std::to_string(white[0]), ", ",
                						  std::to_string(white[1]), ", ",
    									  std::to_string(white[2]), "]"});

    debug({"FoveonHelper: scale [", std::to_string(scale[0]), ", ",
                    				std::to_string(scale[1]), ", ",
        							std::to_string(scale[2]), "]"});

    return true;
} // get_raw_property


bool FoveonHelper::get_crop_area(uint32_t *rect) const{
    /*CameraConstantsStore* ccs = CameraConstantsStore::getInstance();
    CameraConst *cc = ccs->get(raw_image->get_maker().c_str(), raw_image->get_model().c_str());

    if(cc->has_rawCrop(raw_image->get_width(), raw_image->get_height())){
        int lm, tm, w, h;
        cc->get_rawCrop(raw_image->get_width(), raw_image->get_height(), lm, tm, w, h);
        //debug({"cc_crop[", std::to_string(lm), ", ", std::to_string(tm),
        //       ", ", std::to_string(w), ", ", std::to_string(h), "]"});
        rect[0] = tm;
        rect[1] = lm;
        rect[2] = raw_image->get_height() + h;
        rect[3] = raw_image->get_width() + w ;
        debug({"FoveonHelper: Active Image Area II [", std::to_string(rect[0]), ", ",
													  std::to_string(rect[1]), ", ",
													  std::to_string(rect[2]), ", ",
													  std::to_string(rect[3]), "]"});*/
	if (crop_area[2]){
		rect[0] = crop_area[0];
		rect[1] = crop_area[1];
		rect[2] = crop_area[2];
		rect[3] = crop_area[3];
		debug({"FoveonHelper: Active Image Area II [", std::to_string(crop_area[0]), ", ",
													   std::to_string(crop_area[1]), ", ",
													   std::to_string(crop_area[2]), ", ",
													   std::to_string(crop_area[3]), "]"});
    } else {
    	debug({"FoveonHelper: could not get crop area."});
    	return false;
    }
    return true;
} // get_crop_area


std::vector<GainMap> FoveonHelper::read_foveon_spatial_gain(){
    debug({"FoveonHelper: reading foveon spatial gain correction."});
    std::vector<GainMap> gain_maps_vector;

    if (!get_raw_property()){
        debug({"FoveonHelper is not in a valid state."});
    } else {
        /* Quattro raws are already corrected for spatial gain. 
         * So it is disabled by default. */
        if (x3f->header.version < X3F::X3F_VERSION_4_0 || !is_quattro){
            X3F::x3f_spatial_gain_corr_t corr[X3F::MAXCORR];
            int corr_num;
            uint32_t active_area[4];
            double originv, originh, scalev, scaleh;

            corr_num = X3F::x3f_get_spatial_gain(x3f, cam_wb, corr);

            if (corr_num == 0){
                debug({"FoveonHelper: no spatial gain correction found..."});
                X3F::x3f_delete(x3f);
                return gain_maps_vector;
            }
            debug({"FoveonHelper: x3f_tools found ", std::to_string(corr_num), " spatial gain correction(s)..."});
            
            /* Spatial gain in X3F refers to the entire image, 
             * But correction will be applied after cropping to ActiveArea. */
            if (!get_crop_area(active_area)){
                X3F::x3f_delete(x3f);
                return gain_maps_vector;
            } 

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
                    gain_map.map_gain.push_back(float(c->gain[j]));
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
