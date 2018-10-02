/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "improcfun.h"
#include "gauss.h"
#include "sleef.c"
#include "opthelper.h"
#include "guidedfilter.h"

namespace rtengine {

void ImProcFunctions::shadowsHighlights(LabImage *lab)
{
    if (!params->sh.enabled || (!params->sh.highlights && !params->sh.shadows)){
        return;
    }

    const int width = lab->W;
    const int height = lab->H;

    array2D<float> mask(width, height);
    array2D<float> L(width, height);
    const float sigma = params->sh.radius * 5.f / scale;
    LUTf f(65536);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);

    const auto rgb2lab =
        [&](float R, float G, float B, float &l, float &a, float &b) -> void
        {
            float x, y, z;
            Color::rgbxyz(R, G, B, x, y, z, ws);
            Color::XYZ2Lab(x, y, z, l, a, b);
        };

    const auto lab2rgb =
        [&](float l, float a, float b, float &R, float &G, float &B) -> void
        {
            float x, y, z;
            Color::Lab2XYZ(l, a, b, x, y, z);
            Color::xyz2rgb(x, y, z, R, G, B, iws);
        };
    
    const auto apply =
        [&](int amount, int tonalwidth, bool hl) -> void
        {
            const float thresh = tonalwidth * 327.68f;
            const float scale = hl ? (thresh > 0.f ? 0.9f / thresh : 1.f) : thresh * 0.9f;

#ifdef _OPENMP
            #pragma omp parallel if (multiThread)
#endif
            {

#ifdef _OPENMP
                #pragma omp for
#endif
                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        float l = lab->L[y][x];
                        if (hl) {
                            mask[y][x] = (l > thresh) ? 1.f : pow4(l * scale);
                            L[y][x] = 1.f - (l / 32768.f);
                        } else {
                            mask[y][x] = l <= thresh ? 1.f : pow4(scale / l);
                            L[y][x] = l / 32768.f;
                        }
                        // if (L[y][x] < 0.f || L[y][x] > 1.f) {
                        //     std::cerr << "BAD GUIDE!" << std::endl;
                        //     abort();
                        // }
                    }
                }

                //gaussianBlur(mask, mask, width, height, sigma);
            }

            guidedFilter(L, mask, mask, sigma, 0.01, multiThread, 1);

            const float base = std::pow(4.f, float(amount)/100.f);
            const float gamma = hl ? base : 1.f / base;

            const float contrast = std::pow(2.f, float(amount)/100.f);
            DiagonalCurve sh_contrast({
                    DCT_NURBS,
                    0, 0,
                    0.125, std::pow(0.125 / 0.25, contrast) * 0.25, 
                    0.25, 0.25,
                    0.375, std::pow(0.375 / 0.25, contrast) * 0.25,
                    1, 1
                });

            if(!hl) {
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int c = 0; c < 65536; ++c) {
                    float l, a, b;
                    float R = c, G = c, B = c;
                    rgb2lab(R, G, B, l, a, b);
                    auto base = pow_F(l / 32768.f, gamma);
                    // get a bit more contrast in the shadows
                    base = sh_contrast.getVal(base);
                    l = base * 32768.f;
                    lab2rgb(l, a, b, R, G, B);
                    f[c] = G;
                }
            } else {
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int c = 0; c < 65536; ++c) {
                    float l, a, b;
                    float R = c, G = c, B = c;
                    rgb2lab(R, G, B, l, a, b);
                    auto base = pow_F(l / 32768.f, gamma);
                    l = base * 32768.f;
                    lab2rgb(l, a, b, R, G, B);
                    f[c] = G;
                }
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    // float l = lab->L[y][x];
                    float blend = LIM01(mask[y][x]);
                    if (blend < 0.f || blend > 1.f) {
                        abort();
                    }
                    float orig = 1.f - blend;
                    if (lab->L[y][x] >= 0.f && lab->L[y][x] < 32768.f) {
                        float rgb[3];
                        lab2rgb(lab->L[y][x], lab->a[y][x], lab->b[y][x], rgb[0], rgb[1], rgb[2]);
                        for (int i = 0; i < 3; ++i) {
                            rgb[i] = f[rgb[i]] * blend + rgb[i] * orig;
                        }
                        rgb2lab(rgb[0], rgb[1], rgb[2], lab->L[y][x], lab->a[y][x], lab->b[y][x]);
                    }
                }
            }
        };

    if (params->sh.highlights > 0) {
        apply(params->sh.highlights, params->sh.htonalwidth, true);
    }

    if (params->sh.shadows > 0) {
        apply(params->sh.shadows, params->sh.stonalwidth, false);
    }
}

} // namespace rtengine
