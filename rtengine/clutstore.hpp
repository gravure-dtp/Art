/** -*- C++ -*-
 *  
 *  This file is part of ART
 *
 *  Copyright (c) 2022 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "clutstore.h"

namespace rtengine {

inline void HaldCLUTApplication::apply_single(float &r, float &g, float &b)
{
    if (!ok_) {
        return;
    }

#ifdef ART_USE_OCIO
    if (ocio_processor_) {
        Vec3<float> v(r / 65535.f, g / 65535.f, b / 65535.f);
        v = dot_product(conv_, v);

        OCIO::PackedImageDesc pd(v, 1, 1, 3);
        ocio_processor_->apply(pd);
            
        v = dot_product(iconv_, v);
        r = v[0];
        g = v[1];
        b = v[2];
    } else
#endif // ART_USE_OCIO
    {
        float out_rgbx[4] ALIGNED16; // Line buffer for CLUT
        float clutr[1] = {r};
        float clutg[1] = {g};
        float clutb[1] = {b};
    
        if (!clut_and_working_profiles_are_same_) {
            // Convert from working to clut profile
            float x, y, z;
            Color::rgbxyz(r, g, b, x, y, z, wprof_);
            Color::xyz2rgb(x, y, z, clutr[0], clutg[0], clutb[0], xyz2clut_);
        }

        // Apply gamma sRGB (default RT)
        clutr[0] = Color::gamma_srgbclipped(clutr[0]);
        clutg[0] = Color::gamma_srgbclipped(clutg[0]);
        clutb[0] = Color::gamma_srgbclipped(clutb[0]);

        hald_clut_->getRGB(strength_, 1, clutr, clutg, clutb, out_rgbx);

        // Apply inverse gamma sRGB
        clutr[0] = Color::igamma_srgb(out_rgbx[0]);
        clutg[0] = Color::igamma_srgb(out_rgbx[1]);
        clutb[0] = Color::igamma_srgb(out_rgbx[2]);

        if (!clut_and_working_profiles_are_same_) {
            // Convert from clut to working profile
            float x, y, z;
            Color::rgbxyz(clutr[0], clutg[0], clutb[0], x, y, z, clut2xyz_);
            Color::xyz2rgb(x, y, z, clutr[0], clutg[0], clutb[0], wiprof_);
        }

        r = clutr[0];
        g = clutg[0];
        b = clutb[0];
    }
}    

} // namespace rtengine
