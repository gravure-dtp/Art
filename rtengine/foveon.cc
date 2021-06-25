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

#include "settings.h"
#include "myfile.h"
#include "camconst.h"
#include "rawimage.h"
#include "foveon.h"
#include "x3f_tools/x3f_io.h"
#include "x3f_tools/x3f_meta.h"
#include "x3f_tools/x3f_spatial_gain.h"


namespace rtengine{
extern const Settings *settings;


FoveonHelper::FoveonHelper(RawImage const *ri)
:raw_image(ri), x3f(nullptr){
    x3f = x3f_new_from_file(raw_image->get_file());
    if (!x3f){
        debug({"x3f_tools could not read raw ", raw_image->get_filename()});
    } else {
        debug({"FoveonHelper for raw ", raw_image->get_filename()});
        
        x3f_directory_entry_t *DE = x3f_get_prop(x3f);
        if (X3F_OK!=(x3f_load_data(x3f, x3f_get_camf(x3f)))){
            debug({"x3f_tools could not load CAMF"});
            x3f_delete(x3f);
            x3f = nullptr;
        } else{
            debug({"x3f_tools accessing metadata..."}); 
            if (DE!=NULL){ /* Not for Quattro */
                if (X3F_OK != (x3f_load_data(x3f, DE))){
                    debug({"x3f_tools could not load PROP"});
                    x3f_delete(x3f);
                    x3f = nullptr;
                }
            }
            
            /*
            if (NULL==(DE=x3f_get_raw(x3f))){
                std::cout << "x3f_tools could not find any matching RAW format" << std::endl;
                x3f_delete(x3f);
                x3f = nullptr;
            }
            
            if (X3F_OK!=(x3f_load_data(x3f, DE))){
                std::cout << "x3f_tools could not load RAW" << std::endl;
                x3f_delete(x3f);
                x3f = nullptr;
            }*/
        }
    }
}

FoveonHelper::~FoveonHelper(){
    if(x3f){
        x3f_delete(x3f);
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


bool FoveonHelper::get_camf_rect(uint32_t *crect) const{
    CameraConstantsStore* ccs = CameraConstantsStore::getInstance();
    CameraConst *cc = ccs->get(raw_image->get_maker().c_str(), raw_image->get_model().c_str());
    
    if(cc->has_rawCrop(raw_image->get_width(), raw_image->get_height())){
        int lm, tm, w, h;
        cc->get_rawCrop(raw_image->get_width(), raw_image->get_height(), lm, tm, w, h);
        debug({"cc_crop[", std::to_string(lm), ", ", std::to_string(tm),
               ", ", std::to_string(w), ", ", std::to_string(h), "]"});
        crect[0] = tm;
        crect[1] = lm;
        crect[2] = raw_image->get_height() + h;
        crect[3] = raw_image->get_width() + w ; 
    } else return false;
    return true;
}
 
 
std::vector<GainMap> FoveonHelper::read_foveon_spatial_gain() const{    
    debug({"read_foveon_spatial_gain()"});
    std::vector<GainMap> gain_maps_vector;
    
    if (!x3f){
        debug({"FoveonHelper is not in a valid state"});
    } else {
        /* Quattro raws are already corrected for spatial gain. 
         * So it is disabled by default. */
        if (x3f->header.version < X3F_VERSION_4_0){
            x3f_spatial_gain_corr_t corr[MAXCORR];
            int corr_num;
            uint32_t active_area[4];
            double originv, originh, scalev, scaleh;
            
            corr_num = x3f_get_spatial_gain(x3f, x3f_get_wb(x3f), corr);
            if (corr_num == 0){
                debug({"x3f_tools : no spatial gain correction found..."});
                x3f_delete(x3f);
                return gain_maps_vector;
            }
            debug({"x3f_tools found ", std::to_string(corr_num), " spatial gain correction(s)..."});
            
            /* Spatial gain in X3F refers to the entire image, 
             * But correction will be applied after cropping to ActiveArea. */
            if (!get_camf_rect(active_area)){
                debug({"couldn't get camf rect..."});
                x3f_delete(x3f);
                return gain_maps_vector;
            } 
            debug({"camf rect ActiveArea[", std::to_string(active_area[0]), ", ", 
                    std::to_string(active_area[1]), ", ", std::to_string(active_area[2]),
                    ", ", std::to_string(active_area[3]), "]"});
            
            originv = -(double)active_area[0] / (active_area[2] - active_area[0]); // 0
            originh = -(double)active_area[1] / (active_area[3] - active_area[1]); // 0
            scalev = (double)raw_image->get_height() / (active_area[2] - active_area[0]); // 1.0
            scaleh = (double)raw_image->get_width() / (active_area[3] - active_area[1]); // 1.0
            
            /* Compute GainMaps */            
            for (size_t i=0; i<(size_t)corr_num; i++){
                GainMap gain_map;
                x3f_spatial_gain_corr_t *c = &corr[i];
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
            x3f_cleanup_spatial_gain(corr, corr_num);            
        } else {  // Quattro files
            debug({"x3f_tools : Quattro files ares already corrected for spatial gain."});
        }
    }
    return gain_maps_vector;
}    


} // rtengine