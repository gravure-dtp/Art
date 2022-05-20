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

#include "gainmap.h"
#include "rawimagesource.h"
#include "rescale.h"


namespace rtengine {
extern const Settings *settings;


void RawImageSource::apply_foveon_spatial_gain(unsigned short black[4], std::vector<GainMap> maps){
    if (settings->verbose){
        std::cout << "GAIN MAP: found " << maps.size() << " maps..." << std::endl;
    }
    
    int iblack[4];
    float fblack[4];
	float adjust[3], coef[3];
    FoveonHelper *fv = ri->get_foveon_helper();
    fv->get_BlackDev(iblack);
    for (int i = 0; i < 4; ++i){
        fblack[i] = iblack[i];
    }
    fv->get_spatial_gain_adjust(adjust);
//    adjust[0] = 1.1;
//    adjust[1] = 1.3;
//    adjust[2] = 1.0;
    for (int c=0; c<3; c++)
    	coef[c] = adjust[c] - 1.0;
    
    // now we can apply each gain map to raw_data
    array2D<float> mvals;
    
    for (auto &m : maps){
        if (settings->verbose){
            std::cout << "GAIN MAP: RawSource[W:" << W << ", H:" << H << "] plane=" << std::to_string(m.plane) << std::endl;
            std::cout << "GAIN MAP: " << m.to_str() << std::endl;
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
				int z = 3 * x + m.plane;

				f = (f + coef[m.plane]) / adjust[m.plane]; //f *= adjust[m.plane];
				rawData[y][z] = CLIP((rawData[y][z]) * f);
				// Gain with black removal
				//rawData[y][z] = CLIP((rawData[y][z] - fblack[m.plane]) * f + fblack[m.plane]);
            }
        }
    }    
}
} // namespace rtengine
