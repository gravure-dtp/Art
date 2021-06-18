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

#include "myfile.h"
#include "gainmap.h"
#include "rawimage.h"
#include "x3f_tools/x3f_io.h"

namespace rtengine {
    
extern const Settings *settings;
    
void RawImage::read_foveon_spatial_gain(){    
    if (settings->verbose){
        std::cout << "read_foveon_spatial_gain()" << std::endl;
    }
    
    x3f_t *x3f = NULL;
    x3f = x3f_new_from_file(ifp);
    if (x3f==NULL){
      std::cout << "x3f_tools could not read image file" << std::endl;
    } else {
      std::cout << "x3f_tools accessing metadata..." << std::endl;  
    }
    
    // Clean up
    x3f_delete(x3f);
}    
    
} // namespace rtengine