/*
 *  This file is part of RawTherapee.
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but widthITheightOUT ANY widthARRANTY; without even the implied warranty of
 *  MERCheightANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  2010 Emil Martinec <ejmartin@uchicago.edu>
 *
 */
#include <cstddef>
#include "rt_math.h"
#include "labimage.h"
#include "improcfun.h"
#include "cieimage.h"
#include "sleef.c"
#include "opthelper.h"
#include "gauss.h"
#include "rt_algo.h"

using namespace std;

namespace rtengine {

void ImProcFunctions::impulse_nr (LabImage* lab, double thresh)
{
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // impulse noise removal
    // local variables

    int width = lab->W;
    int height = lab->H;

    // buffer for the lowpass image
    // float * lpf[height] ALIGNED16;
    // lpf[0] = new float [width * height];
    // buffer for the highpass image
    char * impish[height] ALIGNED16;
    impish[0] = new char [width * height];

    for (int i = 1; i < height; i++) {
        // lpf[i] = lpf[i - 1] + width;
        impish[i] = impish[i - 1] + width;
    }

    markImpulse(width, height, lab->L, impish, thresh);

//     //The cleaning algorithm starts here

//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     // modified bilateral filter for lowpass image, omitting input pixel; or Gaussian blur

    const float eps = 1.0;

// #ifdef _OPENMP
//     #pragma omp parallel
// #endif
//     {
//         gaussianBlur (lab->L, lpf, width, height, max(2.0, thresh - 1.0));
//     }

//     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//     float impthr = max(1.0, 5.5 - thresh);
//     float impthrDiv24 = impthr / 24.0f;         //Issue 1671: moved the Division outside the loop, impthr can be optimized out too, but I let in the code at the moment


// #ifdef _OPENMP
//     #pragma omp parallel
// #endif
//     {
//         int i1, j1, j;
//         float hpfabs, hfnbrave;
// #ifdef __SSE2__
//         vfloat hfnbravev, hpfabsv;
//         vfloat impthrDiv24v = F2V( impthrDiv24 );
// #endif
// #ifdef _OPENMP
//         #pragma omp for
// #endif

//         for (int i = 0; i < height; i++) {
//             for (j = 0; j < 2; j++) {
//                 hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

//                 //block average of high pass data
//                 for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
//                     for (j1 = 0; j1 <= j + 2; j1++) {
//                         hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
//                     }

//                 impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
//             }

// #ifdef __SSE2__

//             for (; j < width - 5; j += 4) {
//                 hfnbravev = ZEROV;
//                 hpfabsv = vabsf(LVFU(lab->L[i][j]) - LVFU(lpf[i][j]));

//                 //block average of high pass data
//                 for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ ) {
//                     for (j1 = j - 2; j1 <= j + 2; j1++) {
//                         hfnbravev += vabsf(LVFU(lab->L[i1][j1]) - LVFU(lpf[i1][j1]));
//                     }
//                 }

//                 int mask = _mm_movemask_ps((hfnbravev - hpfabsv) * impthrDiv24v - hpfabsv);
//                 impish[i][j] = (mask & 1);
//                 impish[i][j + 1] = ((mask & 2) >> 1);
//                 impish[i][j + 2] = ((mask & 4) >> 2);
//                 impish[i][j + 3] = ((mask & 8) >> 3);
//             }

// #endif

//             for (; j < width - 2; j++) {
//                 hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

//                 //block average of high pass data
//                 for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
//                     for (j1 = j - 2; j1 <= j + 2; j1++) {
//                         hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
//                     }

//                 impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
//             }

//             for (; j < width; j++) {
//                 hpfabs = fabs(lab->L[i][j] - lpf[i][j]);

//                 //block average of high pass data
//                 for (i1 = max(0, i - 2), hfnbrave = 0; i1 <= min(i + 2, height - 1); i1++ )
//                     for (j1 = j - 2; j1 < width; j1++) {
//                         hfnbrave += fabs(lab->L[i1][j1] - lpf[i1][j1]);
//                     }

//                 impish[i][j] = (hpfabs > ((hfnbrave - hpfabs) * impthrDiv24));
//             }
//         }
//     }

//now impulsive values have been identified

// Issue 1671:
// often, noise isn't evenly distributed, e.g. only a few noisy pixels in the bright sky, but many in the dark foreground,
// so it's better to schedule dynamic and let every thread only process 16 rows, to avoid running big threads out of work
// Measured it and in fact gives better performance than without schedule(dynamic,16). Of course, there could be a better
// choice for the chunk_size than 16
// race conditions are avoided by the array impish
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int i1, j1, j;
        float wtdsum[3], dirwt, norm;
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int i = 0; i < height; i++) {
            for (j = 0; j < 2; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = 0; j1 <= j + 2; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1 / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * lab->L[i1][j1];
                        wtdsum[1] += dirwt * lab->a[i1][j1];
                        wtdsum[2] += dirwt * lab->b[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    lab->L[i][j] = wtdsum[0] / norm; //low pass filter
                    lab->a[i][j] = wtdsum[1] / norm; //low pass filter
                    lab->b[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }

            for (; j < width - 2; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 <= j + 2; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1 / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * lab->L[i1][j1];
                        wtdsum[1] += dirwt * lab->a[i1][j1];
                        wtdsum[2] += dirwt * lab->b[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    lab->L[i][j] = wtdsum[0] / norm; //low pass filter
                    lab->a[i][j] = wtdsum[1] / norm; //low pass filter
                    lab->b[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }

            for (; j < width; j++) {
                if (!impish[i][j]) {
                    continue;
                }

                norm = 0.0;
                wtdsum[0] = wtdsum[1] = wtdsum[2] = 0.0;

                for (i1 = max(0, i - 2); i1 <= min(i + 2, height - 1); i1++ )
                    for (j1 = j - 2; j1 < width; j1++ ) {
                        if (impish[i1][j1]) {
                            continue;
                        }

                        dirwt = 1 / (SQR(lab->L[i1][j1] - lab->L[i][j]) + eps); //use more sophisticated rangefn???
                        wtdsum[0] += dirwt * lab->L[i1][j1];
                        wtdsum[1] += dirwt * lab->a[i1][j1];
                        wtdsum[2] += dirwt * lab->b[i1][j1];
                        norm += dirwt;
                    }

                if (norm) {
                    lab->L[i][j] = wtdsum[0] / norm; //low pass filter
                    lab->a[i][j] = wtdsum[1] / norm; //low pass filter
                    lab->b[i][j] = wtdsum[2] / norm; //low pass filter
                }
            }
        }
    }
//now impulsive values have been corrected

    // delete [] lpf[0];
    delete [] impish[0];

}

} // namespace rtengine
