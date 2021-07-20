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
#pragma once

#include <vector>
 
#include "gainmap.h"
#include "rawimage.h"

#include "x3f_tools/x3f_io.h"
#include "x3f_tools/x3f_meta.h"

namespace rtengine {

class RawImage; // forward declaration avoiding circular ref 
    
class FoveonHelper{
public:
    explicit FoveonHelper(RawImage const *ri);
    ~FoveonHelper();
    
protected:
    RawImage const *raw_image;
    x3ftools::x3f_t *x3f;
    bool data_loaded;
    int black[3];
    int white[3];
    
public:
    std::vector<GainMap> read_foveon_spatial_gain();
    bool get_camf_rect(uint32_t *crect) const;
    bool get_BlackLevels(int* black_c4);
    bool get_WhiteLevels(int* white_c4);
    bool get_ccMatrix(double *matrix);
    bool get_cam_xyz(double matrix[4][3]);
    bool get_asShotNeutral(float *gain);
    bool get_wbGain(float *gain);
    
    static bool is_supported(const std::string& camera){
        static const std::vector<std::string> supported_cams = 
            {"DP1 Merrill", "DP2 Merrill", "DP3 Merrill"};
        for(size_t i=0; i< supported_cams.size(); i++){
            if (supported_cams[i] == camera)
                return true;
            else continue;
        }
        return false;
    }
    
protected:
    void debug(std::vector<std::string> mess) const;
    bool load_data();
    bool get_Levels();
};


} // namespace rtengine

