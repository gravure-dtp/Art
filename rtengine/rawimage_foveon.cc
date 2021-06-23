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

#include "myfile.h"
#include "gainmap.h"
#include "rawimage.h"
#include "rawimagesource.h"
#include "camconst.h"
#include "rescale.h"
#include "x3f_tools/x3f_io.h"
#include "x3f_tools/x3f_meta.h"
#include "x3f_tools/x3f_spatial_gain.h"

namespace rtengine {
    
extern const Settings *settings;

void RawImageSource::apply_foveon_spatial_gain(unsigned short black[4], std::vector<GainMap> maps){
    if (settings->verbose){
        std::cout << "GAIN MAP: found " << maps.size() << " maps..." << std::endl;
    }
    
    float fblack[4];
    for (int i = 0; i < 4; ++i){
        fblack[i] = black[i];
    }
    
    // now we can apply each gain map to raw_data
    array2D<float> mvals;
    
    for (auto &m : maps){
        if (settings->verbose){
            std::cout << "GAIN MAP: RawSource[W:" << W << ", H:" << H << "]" << std::endl;
            //std::cout << "GAIN MAP: " << m.to_str() << std::endl;
        }
        
        mvals(m.map_points_h, m.map_points_v, &(m.map_gain[0]), 0);
        const float col_scale = float(m.map_points_h-1) / float(W);
        const float row_scale = float(m.map_points_v-1) / float(H);

#ifdef _OPENMP
#       pragma omp parallel for
#endif
        for (unsigned y = m.top; y < m.bottom; y += m.row_pitch){
            float ys = y * row_scale;
            for (unsigned x = m.left; x < m.right; x += m.col_pitch){
                float xs = x * col_scale;
                float f = getBilinearValue(mvals, xs, ys);
                int c  = FC(y, x);
                int c4 = ( c == 1 && !(y & 1) ) ? 3 : c;
                float b = fblack[c4];
                int z = 3 * x + m.plane;
                rawData[y][z] = CLIP((rawData[y][z] - b) * f + b);
            }
        }
    }    
}


bool RawImage::get_foveon_camf_rect(uint32_t *crect){
    CameraConstantsStore* ccs = CameraConstantsStore::getInstance();
    CameraConst *cc = ccs->get(make, model);
    if(cc->has_rawCrop(width, height)){
        int lm, tm, w, h;
        cc->get_rawCrop(width, height, lm, tm, w, h);
        if (settings->verbose)
            std::cout << "cc_crop[" << lm
                              << ", " << tm
                              << ", " << w
                              << ", " << h << "]" << std::endl;
        crect[0] = tm;
        crect[1] = lm;
        crect[2] = height + h;
        crect[3] = width + w ; 
    } else return false;
    return true;
}

    
std::vector<GainMap> RawImage::read_foveon_spatial_gain(){    
    if (settings->verbose)
        std::cout << "read_foveon_spatial_gain()" << std::endl;
                
    std::vector<GainMap> gain_maps_vector;
    x3f_t *x3f = NULL;
    x3f = x3f_new_from_file(ifp);
    
    if (x3f==NULL){
        if (settings->verbose)
            std::cout << "x3f_tools could not read image file" << std::endl;
    } else {
        if (settings->verbose)
            std::cout << "x3f_tools accessing metadata..." << std::endl; 

        x3f_directory_entry_t *DE = x3f_get_prop(x3f);
        
        if (X3F_OK!=(x3f_load_data(x3f, x3f_get_camf(x3f)))){
            if (settings->verbose)
                std::cout << "x3f_tools could not load CAMF" << std::endl;
            x3f_delete(x3f);
            return gain_maps_vector;
        }
        
        //TODO: is this step mandatory?
        if (DE!=NULL){
            /* Not for Quattro */
            if (X3F_OK != (x3f_load_data(x3f, DE))){
                if (settings->verbose)
                    std::cout << "x3f_tools could not load PROP" << std::endl;
                x3f_delete(x3f);
                return gain_maps_vector;
            }
        }
        /*
        if (NULL==(DE=x3f_get_raw(x3f))){
            std::cout << "x3f_tools could not find any matching RAW format" << std::endl;
            x3f_delete(x3f);
            return gain_maps_vector;
        }
        
        if (X3F_OK!=(x3f_load_data(x3f, DE))){
            std::cout << "x3f_tools could not load RAW" << std::endl;
            x3f_delete(x3f);
            return gain_maps_vector;
        }*/
        
        /* Quattro raws are already corrected for spatial gain. 
         * So it is disabled by default. */
        if (x3f->header.version < X3F_VERSION_4_0){
            x3f_spatial_gain_corr_t corr[MAXCORR];
            int corr_num;
            uint32_t active_area[4];
            double originv, originh, scalev, scaleh;
            
            corr_num = x3f_get_spatial_gain(x3f, x3f_get_wb(x3f), corr);
            if (corr_num == 0){
                if (settings->verbose)
                    std::cout << "x3f_tools : no spatial gain correction found..." << std::endl;
                x3f_delete(x3f);
                return gain_maps_vector;
            }
            
            if (settings->verbose)
                std::cout << "x3f_tools found " << corr_num << " spatial gain correction(s)..." << std::endl;
            
            /* Spatial gain in X3F refers to the entire image, 
             * But correction will be applied after cropping to ActiveArea. */
            if (!get_foveon_camf_rect(active_area)){
                if (settings->verbose)
                    std::cout << "couldn't get camf rect..." << std::endl;
                x3f_delete(x3f);
                return gain_maps_vector;
            } else if (settings->verbose){
                std::cout << "camf rect ActiveArea[" << active_area[0]
                          << ", " << active_area[1]
                          << ", " << active_area[2]
                          << ", " << active_area[3] << "]" << std::endl;
            }
            
            originv = -(double)active_area[0] / (active_area[2] - active_area[0]); // 0
            originh = -(double)active_area[1] / (active_area[3] - active_area[1]); // 0
            scalev = (double)this->get_height() / (active_area[2] - active_area[0]); // 1.0
            scaleh = (double)this->get_width() / (active_area[3] - active_area[1]); // 1.0
            
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
            std::cout << "x3f_tools : Quattro files ares already corrected for spatial gain." << std::endl;
        }
    }
    // Clean up
    x3f_delete(x3f);
    return gain_maps_vector;
}    
    
} // namespace rtengine