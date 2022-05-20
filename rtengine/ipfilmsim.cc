/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "improcfun.h"
#include "curves.h"
#include "color.h"
#include "clutstore.h"
#include "../rtgui/multilangmgr.h"

namespace rtengine {

void ImProcFunctions::filmSimulation(Imagefloat *img)
{
    if (!params->filmSimulation.enabled) {
        return;
    }

    img->setMode(Imagefloat::Mode::RGB, multiThread);

    HaldCLUTApplication hald_clut(params->filmSimulation.clutFilename, params->icm.workingProfile);
    constexpr int TS = 112;
    hald_clut.init(float(params->filmSimulation.strength)/100.f, TS);

    if (hald_clut) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < img->getHeight(); ++y) {
            for (int jj = 0; jj < img->getWidth(); jj += TS) {
                int jstart = jj;
                float *r = img->r(y)+jstart;
                float *g = img->g(y)+jstart;
                float *b = img->b(y)+jstart;
                int tW = min(jj + TS, img->getWidth());
                hald_clut(r, g, b, 0, jstart, tW, 1);
            }
        }
    } else if (plistener) {
        plistener->error(Glib::ustring::compose(M("TP_FILMSIMULATION_LABEL") + " - " + M("ERROR_MSG_FILE_READ"), params->filmSimulation.clutFilename.empty() ? "(" + M("GENERAL_NONE") + ")" : params->filmSimulation.clutFilename));
    }
}

} // namespace rtengine
