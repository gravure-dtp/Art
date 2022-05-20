/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include <cmath>
#include <glib.h>
#include <glibmm.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "alignedbuffer.h"
#include "rtengine.h"
#include "improcfun.h"
#include "curves.h"
#include "mytime.h"
#include "iccstore.h"
#include "imagesource.h"
#include "rtthumbnail.h"
#include "utils.h"
#include "iccmatrices.h"
#include "color.h"
#include "calc_distort.h"
#include "rt_math.h"
#include "improccoordinator.h"
#include "clutstore.h"
#include "StopWatch.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/guiutils.h"
#include "refreshmap.h"

namespace rtengine {

using namespace procparams;

extern const Settings* settings;

ImProcFunctions::ImProcFunctions(const ProcParams* iparams, bool imultiThread):
    monitor(nullptr),
    monitorTransform(nullptr),
    params(iparams),
    scale(1),
    multiThread(imultiThread),
    cur_pipeline(Pipeline::OUTPUT),
    dcpProf(nullptr),
    dcpApplyState(nullptr),
    pipetteBuffer(nullptr),
    lumimul{},
    offset_x(0),
    offset_y(0),
    full_width(-1),
    full_height(-1),
    histToneCurve(nullptr),
    histCCurve(nullptr),
    histLCurve(nullptr),
    show_sharpening_mask(false),
    plistener(nullptr),
    progress_step(0),
    progress_end(1)
{
}


ImProcFunctions::~ImProcFunctions ()
{
    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }
}

void ImProcFunctions::setScale (double iscale)
{
    scale = iscale;
}


void ImProcFunctions::updateColorProfiles (const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, GamutCheck gamutCheck)
{
    // set up monitor transform
    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }
    gamutWarning.reset(nullptr);

    monitorTransform = nullptr;
    monitor = nullptr;

    if (!monitorProfile.empty()) {
#if !defined(__APPLE__) // No support for monitor profiles on OS X, all data is sRGB
        monitor = ICCStore::getInstance()->getProfile (monitorProfile);
#else
        monitor = ICCStore::getInstance()->getProfile (options.rtSettings.srgb);
#endif
    }

    if (monitor) {
        MyMutex::MyLock lcmsLock (*lcmsMutex);

        cmsUInt32Number flags;
        //cmsHPROFILE iprof  = cmsCreateLab4Profile (nullptr);
        cmsHPROFILE iprof = ICCStore::getInstance()->getProfile(params->icm.outputProfile);
        if (!iprof) {
            iprof = ICCStore::getInstance()->getsRGBProfile();
        }
        
        cmsHPROFILE gamutprof = nullptr;
        cmsUInt32Number gamutbpc = 0;
        RenderingIntent gamutintent = RI_RELATIVE;

        bool softProofCreated = false;
        cmsHPROFILE oprof = nullptr;
        RenderingIntent outIntent;

        if (softProof) {
            flags = cmsFLAGS_SOFTPROOFING | cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (!settings->printerProfile.empty()) {
                oprof = ICCStore::getInstance()->getProfile (settings->printerProfile);
                if (settings->printerBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }
                outIntent = settings->printerIntent;
            // } else {
            //     oprof = ICCStore::getInstance()->getProfile(params->icm.outputProfile);
            //     if (params->icm.outputBPC) {
            //         flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            //     }
            //     outIntent = params->icm.outputIntent;
            }

            if (oprof) {
                // NOCACHE is for thread safety, NOOPTIMIZE for precision

                // if (gamutCheck) {
                //     flags |= cmsFLAGS_GAMUTCHECK;
                // }

                const auto make_gamma_table =
                    [](cmsHPROFILE prof, cmsTagSignature tag) -> void
                    {
                        cmsToneCurve *tc = static_cast<cmsToneCurve *>(cmsReadTag(prof, tag));
                        if (tc) {
                            const cmsUInt16Number *table = cmsGetToneCurveEstimatedTable(tc);
                            cmsToneCurve *tc16 = cmsBuildTabulatedToneCurve16(nullptr, cmsGetToneCurveEstimatedTableEntries(tc), table);
                            if (tc16) {
                                cmsWriteTag(prof, tag, tc16);
                                cmsFreeToneCurve(tc16);
                            }
                        }
                    };

                cmsHPROFILE softproof = ProfileContent(oprof).toProfile();
                if (softproof) {
                    make_gamma_table(softproof, cmsSigRedTRCTag);
                    make_gamma_table(softproof, cmsSigGreenTRCTag);
                    make_gamma_table(softproof, cmsSigBlueTRCTag);
                }

                monitorTransform = cmsCreateProofingTransform (
                    iprof, TYPE_RGB_FLT, //TYPE_Lab_FLT,
                                       monitor, TYPE_RGB_FLT,
                                       softproof,
                                       monitorIntent, outIntent,
                                       flags
                                   );

                if (softproof) {
                    cmsCloseProfile(softproof);
                }

                if (monitorTransform) {
                    softProofCreated = true;
                }

                // if (gamutCheck == GAMUT_CHECK_OUTPUT) {
                //     gamutprof = oprof;
                //     if (params->icm.outputBPC) {
                //         gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
                //     }
                //     gamutintent = outIntent;
                // }
            }
        }

        if (gamutCheck == GAMUT_CHECK_MONITOR) {
            gamutprof = monitor;
            if (settings->monitorBPC) {
                gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
            }
            gamutintent = monitorIntent;
        } else if (gamutCheck == GAMUT_CHECK_OUTPUT) {
            if (!oprof) {
                oprof = ICCStore::getInstance()->getProfile(params->icm.outputProfile);
                if (params->icm.outputBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }
                outIntent = params->icm.outputIntent;
            }
            gamutprof = oprof;
            if (params->icm.outputBPC) {
                gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
            }
            gamutintent = outIntent;
        }

        if (!softProofCreated) {
            flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (settings->monitorBPC) {
                flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            }

            monitorTransform = cmsCreateTransform (iprof, TYPE_RGB_FLT, 
                                                   monitor, TYPE_RGB_FLT, monitorIntent, flags);
        }

        if (gamutCheck && gamutprof) {
            gamutWarning.reset(new GamutWarning(gamutprof, gamutintent, gamutbpc));
        }

//        cmsCloseProfile (iprof);
    }
}

void ImProcFunctions::firstAnalysis (const Imagefloat* const original, const ProcParams &params, LUTu & histogram)
{

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.workingProfile);

    lumimul[0] = wprof[1][0];
    lumimul[1] = wprof[1][1];
    lumimul[2] = wprof[1][2];
    int W = original->getWidth();
    int H = original->getHeight();

    float lumimulf[3] = {static_cast<float> (lumimul[0]), static_cast<float> (lumimul[1]), static_cast<float> (lumimul[2])};

    // calculate histogram of the y channel needed for contrast curve calculation in exposure adjustments
    histogram.clear();

    if (multiThread) {

#ifdef _OPENMP
        const int numThreads = min (max (W * H / (int)histogram.getSize(), 1), omp_get_max_threads());
        #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
        {
            LUTu hist (histogram.getSize());
            hist.clear();
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {

                    float r = original->r (i, j);
                    float g = original->g (i, j);
                    float b = original->b (i, j);

                    int y = (lumimulf[0] * r + lumimulf[1] * g + lumimulf[2] * b);
                    hist[y]++;
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            histogram += hist;

        }
#ifdef _OPENMP
        static_cast<void> (numThreads); // to silence cppcheck warning
#endif
    } else {
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {

                float r = original->r (i, j);
                float g = original->g (i, j);
                float b = original->b (i, j);

                int y = (lumimulf[0] * r + lumimulf[1] * g + lumimulf[2] * b);
                histogram[y]++;
            }
        }
    }
}


namespace {

void proPhotoBlue(Imagefloat *rgb, bool multiThread)
{
    // TODO
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            vfloat rv = LVF(rgb->r(y, x));
            vfloat gv = LVF(rgb->g(y, x));
            vmask zeromask = vorm(vmaskf_eq(rv, ZEROV), vmaskf_eq(gv, ZEROV));
            if (_mm_movemask_ps((vfloat)zeromask)) {
                for (int k = 0; k < 4; ++k) {
                    float r = rgb->r(y, x+k);
                    float g = rgb->g(y, x+k);
                    float b = rgb->b(y, x+k);
                    
                    if ((r == 0.0f || g == 0.0f) && rtengine::min(r, g, b) >= 0.f) {
                        float h, s, v;
                        Color::rgb2hsv(r, g, b, h, s, v);
                        s *= 0.99f;
                        Color::hsv2rgb(h, s, v, rgb->r(y, x+k), rgb->g(y, x+k), rgb->b(y, x+k));
                    }
                }
            }
        }
#endif
        for (; x < W; ++x) {
            float r = rgb->r(y, x);
            float g = rgb->g(y, x);
            float b = rgb->b(y, x);

            if ((r == 0.0f || g == 0.0f) && rtengine::min(r, g, b) >= 0.f) {
                float h, s, v;
                Color::rgb2hsv(r, g, b, h, s, v);
                s *= 0.99f;
                Color::hsv2rgb(h, s, v, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x));
            }
        }
    }
}


void dcpProfile(Imagefloat *img, DCPProfile *dcp, const DCPProfile::ApplyState *as, bool multithread)
{
    if (dcp && as) {
        img->setMode(Imagefloat::Mode::RGB, multithread);
        
        const int H = img->getHeight();
        const int W = img->getWidth();
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            float *r = img->r(y);
            float *g = img->g(y);
            float *b = img->b(y);
            dcp->step2ApplyTile(r, g, b, W, 1, 1, *as);
        }
    }
}

} // namespace


void ImProcFunctions::getAutoExp  (const LUTu &histogram, int histcompr, double clip,
                                   double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh)
{

    float scale = 65536.0f;
    float midgray = 0.1842f; //middle gray in linear gamma =1 50% luminance

    int imax = 65536 >> histcompr;
    int overex = 0;
    float sum = 0.f, hisum = 0.f, losum = 0.f;
    float ave = 0.f, hidev = 0.f, lodev = 0.f;

    //find average luminance
    histogram.getSumAndAverage (sum, ave);

    //find median of luminance
    size_t median = 0, count = histogram[0];

    while (count < sum / 2) {
        median++;
        count += histogram[median];
    }

    if (median == 0 || ave < 1.f) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }

    // compute std dev on the high and low side of median
    // and octiles of histogram
    float octile[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}, ospread = 0.f;
    count = 0;

    int i = 0;

    for (; i < min ((int)ave, imax); i++) {
        if (count < 8) {
            octile[count] += histogram[i];

            if (octile[count] > sum / 8.f || (count == 7 && octile[count] > sum / 16.f)) {
                octile[count] = xlog (1. + (float)i) / log (2.f);
                count++;// = min(count+1,7);
            }
        }

        //lodev += SQR(ave-i)*histogram[i];
        lodev += (xlog (ave + 1.f) - xlog ((float)i + 1.)) * histogram[i];
        losum += histogram[i];
    }

    for (; i < imax; i++) {
        if (count < 8) {
            octile[count] += histogram[i];

            if (octile[count] > sum / 8.f || (count == 7 && octile[count] > sum / 16.f)) {
                octile[count] = xlog (1. + (float)i) / log (2.f);
                count++;// = min(count+1,7);
            }
        }

        //hidev += SQR(i-ave)*histogram[i];
        hidev += (xlog ((float)i + 1.) - xlog (ave + 1.f)) * histogram[i];
        hisum += histogram[i];

    }

    if (losum == 0 || hisum == 0) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }

//    lodev = (lodev / (log(2.f) * losum));
//    hidev = (hidev / (log(2.f) * hisum));

    if (octile[6] > log1p ((float)imax) / log2 (2.f)) { //if very overxposed image
        octile[6] = 1.5f * octile[5] - 0.5f * octile[4];
        overex = 2;
    }

    if (octile[7] > log1p ((float)imax) / log2 (2.f)) { //if overexposed
        octile[7] = 1.5f * octile[6] - 0.5f * octile[5];
        overex = 1;
    }

    // store values of octile[6] and octile[7] for calculation of exposure compensation
    // if we don't do this and the pixture is underexposed, calculation of exposure compensation assumes
    // that it's overexposed and calculates the wrong direction
    float oct6, oct7;
    oct6 = octile[6];
    oct7 = octile[7];


    for (int i = 1; i < 8; i++) {
        if (octile[i] == 0.0f) {
            octile[i] = octile[i - 1];
        }
    }

    // compute weighted average separation of octiles
    // for future use in contrast setting
    for (int i = 1; i < 6; i++) {
        ospread += (octile[i + 1] - octile[i]) / max (0.5f, (i > 2 ? (octile[i + 1] - octile[3]) : (octile[3] - octile[i])));
    }

    ospread /= 5.f;

    if (ospread <= 0.f) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }


    // compute clipping points based on the original histograms (linear, without exp comp.)
    unsigned int clipped = 0;
    int rawmax = (imax) - 1;

    while (histogram[rawmax] + clipped <= 0 && rawmax > 1) {
        clipped += histogram[rawmax];
        rawmax--;
    }

    //compute clipped white point
    unsigned int clippable = (int) (sum * clip / 100.f );
    clipped = 0;
    int whiteclip = (imax) - 1;

    while (whiteclip > 1 && (histogram[whiteclip] + clipped) <= clippable) {
        clipped += histogram[whiteclip];
        whiteclip--;
    }

    //compute clipped black point
    clipped = 0;
    int shc = 0;

    while (shc < whiteclip - 1 && histogram[shc] + clipped <= clippable) {
        clipped += histogram[shc];
        shc++;
    }

    //rescale to 65535 max
    rawmax <<= histcompr;
    whiteclip <<= histcompr;
    ave = ave * (1 << histcompr);
    median <<= histcompr;
    shc <<= histcompr;

//    //prevent division by 0
//    if (lodev == 0.f) {
//        lodev = 1.f;
//    }

    //compute exposure compensation as geometric mean of the amount that
    //sets the mean or median at middle gray, and the amount that sets the estimated top
    //of the histogram at or near clipping.
    //float expcomp1 = (log(/*(median/ave)*//*(hidev/lodev)*/midgray*scale/(ave-shc+midgray*shc))+log((hidev/lodev)))/log(2.f);
    float expcomp1 = (log (/*(median/ave)*//*(hidev/lodev)*/midgray * scale / (ave - shc + midgray * shc))) / log (2.f);
    float expcomp2;

    if (overex == 0) { // image is not overexposed
        expcomp2 = 0.5f * ( (15.5f - histcompr - (2.f * oct7 - oct6)) + log (scale / rawmax) / log (2.f) );
    } else {
        expcomp2 = 0.5f * ( (15.5f - histcompr - (2.f * octile[7] - octile[6])) + log (scale / rawmax) / log (2.f) );
    }

    if (fabs (expcomp1) - fabs (expcomp2) > 1.f) { //for great expcomp
        expcomp = (expcomp1 * fabs (expcomp2) + expcomp2 * fabs (expcomp1)) / (fabs (expcomp1) + fabs (expcomp2));
    } else {
        expcomp = 0.5 * (double)expcomp1 + 0.5 * (double) expcomp2; //for small expcomp
    }

    float gain = exp ((float)expcomp * log (2.f));

    float corr = sqrt (gain * scale / rawmax);
    black = (int) shc * corr;


    //now tune hlcompr to bring back rawmax to 65535
    hlcomprthresh = 0;
    //this is a series approximation of the actual formula for comp,
    //which is a transcendental equation
    float comp = (gain * ((float)whiteclip) / scale - 1.f) * 2.3f; // 2.3 instead of 2 to increase slightly comp
    hlcompr = (int) (100.*comp / (max (0.0, expcomp) + 1.0));
    hlcompr = max (0, min (100, hlcompr));

    //now find brightness if gain didn't bring ave to midgray using
    //the envelope of the actual 'control cage' brightness curve for simplicity
    float midtmp = gain * sqrt (median * ave) / scale;

    if (midtmp < 0.1f) {
        bright = (midgray - midtmp) * 15.0 / (midtmp);
    } else {
        bright = (midgray - midtmp) * 15.0 / (0.10833 - 0.0833 * midtmp);
    }

    bright = 0.25 */*(median/ave)*(hidev/lodev)*/max (0, bright);

    //compute contrast that spreads the average spacing of octiles
    contr = (int) 50.0f * (1.1f - ospread);
    contr = max (0, min (100, contr));
    //take gamma into account
    double whiteclipg = (int) (CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);

    float gavg = 0.;

    float val = 0.f;
    float increment = corr * (1 << histcompr);

    for (int i = 0; i < 65536 >> histcompr; i++) {
        gavg += histogram[i] * Color::gamma2curve[val];
        val += increment;
    }

    gavg /= sum;

    if (black < gavg) {
        int maxwhiteclip = (gavg - black) * 4 / 3 + black; // don't let whiteclip be such large that the histogram average goes above 3/4

        if (whiteclipg < maxwhiteclip) {
            whiteclipg = maxwhiteclip;
        }
    }

    whiteclipg = CurveFactory::igamma2 ((float) (whiteclipg / 65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

    //corection with gamma
    black = (int) ((65535 * black) / whiteclipg);
    //expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

    //diagnostics
    //printf ("**************** AUTO LEVELS ****************\n");
    /*
    if (settings->verbose) {
        printf("sum=%i clip=%f clippable=%i  clipWh=%i  clipBl=%i\n",somm, clip, clippable,clipwh, clipbl);
        printf ("expcomp1= %f   expcomp2= %f gain= %f  expcomp=%f\n",expcomp1,expcomp2,gain,expcomp);
        printf ("expo=%f\n",expo);
        printf ("median: %i  average: %f    median/average: %f\n",median,ave, median/ave);
        printf ("average: %f\n",ave);
        printf("comp=%f hlcomp: %i\n",comp, hlcompr);
        printf ("median/average: %f\n",median/ave);
        printf ("lodev: %f   hidev: %f      hidev/lodev: %f\n",lodev,hidev,hidev/lodev);
        printf ("lodev: %f\n",lodev);
        printf ("hidev: %f\n",hidev);
        printf ("imax=%d rawmax= %d  whiteclip= %d  gain= %f\n",imax,rawmax,whiteclip,gain);
        printf ("octiles: %f %f %f %f %f %f %f %f\n",octile[0],octile[1],octile[2],octile[3],octile[4],octile[5],octile[6],octile[7]);
        printf ("ospread= %f\n",ospread);
        printf ("overexp= %i\n",overex);
    }
    */
    /*
     // %%%%%%%%%% LEGACY AUTOEXPOSURE CODE %%%%%%%%%%%%%
     // black point selection is based on the linear result (yielding better visual results)
     black = (int)(shc * corr);
     // compute the white point of the exp. compensated gamma corrected image
     double whiteclipg = (int)(CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);

     // compute average intensity of the exp compensated, gamma corrected image
     double gavg = 0;
     for (int i=0; i<65536>>histcompr; i++)
     gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;

     if (black < gavg) {
     int maxwhiteclip = (gavg - black) * 4 / 3 + black; // don't let whiteclip be such large that the histogram average goes above 3/4
     //double mavg = 65536.0 / (whiteclipg-black) * (gavg - black);
     if (whiteclipg < maxwhiteclip)
     whiteclipg = maxwhiteclip;
     }

     whiteclipg = CurveFactory::igamma2 ((float)(whiteclipg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

     black = (int)((65535*black)/whiteclipg);
     expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

     if (expcomp<0.0)   expcomp = 0.0;*/
    if (expcomp < -5.0) {
        expcomp = -5.0;
    }

    if (expcomp > 12.0) {
        expcomp = 12.0;
    }

    bright = max (-100, min (bright, 100));

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double ImProcFunctions::getAutoDistor  (const Glib::ustring &fname, int thumb_size)
{
    if (fname != "") {
        int w_raw = -1, h_raw = thumb_size;
        int w_thumb = -1, h_thumb = thumb_size;

        eSensorType sensorType = rtengine::ST_NONE;
        Thumbnail* thumb = rtengine::Thumbnail::loadQuickFromRaw (fname, sensorType, w_thumb, h_thumb, 1, FALSE);

        if (!thumb) {
            return 0.0;
        }

        Thumbnail* raw =   rtengine::Thumbnail::loadFromRaw(fname, sensorType, w_raw, h_raw, 1, 1.0, FALSE);

        if (!raw) {
            delete thumb;
            return 0.0;
        }

        if (h_thumb != h_raw) {
            delete thumb;
            delete raw;
            return 0.0;
        }

        int width;

        if (w_thumb > w_raw) {
            width = w_raw;
        } else {
            width = w_thumb;
        }

        unsigned char* thumbGray;
        unsigned char* rawGray;
        thumbGray = thumb->getGrayscaleHistEQ (width);
        rawGray = raw->getGrayscaleHistEQ (width);

        if (!thumbGray || !rawGray) {
            if (thumbGray) {
                delete thumbGray;
            }

            if (rawGray) {
                delete rawGray;
            }

            delete thumb;
            delete raw;
            return 0.0;
        }

        double dist_amount;
        int dist_result = calcDistortion (thumbGray, rawGray, width, h_thumb, 1, dist_amount);

        if (dist_result == -1) { // not enough features found, try increasing max. number of features by factor 4
            calcDistortion (thumbGray, rawGray, width, h_thumb, 4, dist_amount);
        }

        delete[] thumbGray;
        delete[] rawGray;
        delete thumb;
        delete raw;
        return dist_amount;
    } else {
        return 0.0;
    }
}

void ImProcFunctions::rgb2lab (Imagefloat &src, LabImage &dst, const Glib::ustring &workingSpace)
{
    src.assignColorSpace(workingSpace);
    src.toLab(dst, true);
}

void ImProcFunctions::lab2rgb (const LabImage &src, Imagefloat &dst, const Glib::ustring &workingSpace)
{
    dst.assignColorSpace(workingSpace);
    dst.assignMode(Imagefloat::Mode::RGB);
    
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix ( workingSpace );
    const float wip[3][3] = {
        {static_cast<float> (wiprof[0][0]), static_cast<float> (wiprof[0][1]), static_cast<float> (wiprof[0][2])},
        {static_cast<float> (wiprof[1][0]), static_cast<float> (wiprof[1][1]), static_cast<float> (wiprof[1][2])},
        {static_cast<float> (wiprof[2][0]), static_cast<float> (wiprof[2][1]), static_cast<float> (wiprof[2][2])}
    };

    const int W = dst.getWidth();
    const int H = dst.getHeight();
#ifdef __SSE2__
    vfloat wipv[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            wipv[i][j] = F2V (wiprof[i][j]);
        }
    }

#endif

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = 0; i < H; i++) {
        int j = 0;
#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            vfloat X, Y, Z;
            vfloat R, G, B;
            Color::Lab2XYZ (LVFU (src.L[i][j]), LVFU (src.a[i][j]), LVFU (src.b[i][j]), X, Y, Z);
            Color::xyz2rgb (X, Y, Z, R, G, B, wipv);
            STVFU (dst.r (i, j), R);
            STVFU (dst.g (i, j), G);
            STVFU (dst.b (i, j), B);
        }

#endif

        for (; j < W; j++) {
            float X, Y, Z;
            Color::Lab2XYZ (src.L[i][j], src.a[i][j], src.b[i][j], X, Y, Z);
            Color::xyz2rgb (X, Y, Z, dst.r (i, j), dst.g (i, j), dst.b (i, j), wip);
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void ImProcFunctions::setViewport(int ox, int oy, int fw, int fh)
{
    offset_x = ox;
    offset_y = oy;
    full_width = fw;
    full_height = fh;
}


void ImProcFunctions::setOutputHistograms(LUTu *histToneCurve, LUTu *histCCurve, LUTu *histLCurve)
{
    this->histToneCurve = histToneCurve;
    this->histCCurve = histCCurve;
    this->histLCurve = histLCurve;
}


void ImProcFunctions::setShowSharpeningMask(bool yes)
{
    show_sharpening_mask = yes;
}


namespace {

constexpr int NUM_PIPELINE_STEPS = 23;

} // namespace

void ImProcFunctions::setProgressListener(ProgressListener *pl, int num_previews)
{
    plistener = pl;
    progress_step = 0;
    progress_end = NUM_PIPELINE_STEPS * std::max(num_previews, 1);
    if (plistener) {
        plistener->setProgressStr("PROGRESSBAR_PROCESSING");
        plistener->setProgress(0);
    }
}


template <class Ret, class Method>
Ret ImProcFunctions::apply(Method op, Imagefloat *img)
{
    if (plistener) {
        float percent = float(++progress_step) / float(progress_end);
        plistener->setProgress(percent);
    }
    return (this->*op)(img);
}


bool ImProcFunctions::process(Pipeline pipeline, Stage stage, Imagefloat *img)
{
    bool stop = false;
    cur_pipeline = pipeline;

#define STEP_(op) apply<void>(&ImProcFunctions::op, img)
#define STEP_s_(op) apply<bool>(&ImProcFunctions::op, img)
        
    switch (stage) {
    case Stage::STAGE_0:
        STEP_(dehaze);
        STEP_(dynamicRangeCompression);
        break;
    case Stage::STAGE_1:
        STEP_(channelMixer);
        STEP_(exposure);
        STEP_(hslEqualizer);
        stop = STEP_s_(toneEqualizer);
        if (params->icm.workingProfile == "ProPhoto") {
            proPhotoBlue(img, multiThread);
        }
        break;
    case Stage::STAGE_2:
        if (pipeline == Pipeline::OUTPUT ||
            (pipeline == Pipeline::PREVIEW /*&& scale == 1*/)) {
            stop = STEP_s_(sharpening);
            if (!stop) {
                STEP_(impulsedenoise);
                STEP_(defringe);
            }
        }
        stop = stop || STEP_s_(colorCorrection);
        stop = stop || STEP_s_(guidedSmoothing);
        break;
    case Stage::STAGE_3:
        STEP_(creativeGradients);
        stop = stop || STEP_s_(textureBoost);
        if (!stop) { 
            STEP_(logEncoding);
            STEP_(saturationVibrance);
            dcpProfile(img, dcpProf, dcpApplyState, multiThread);
            STEP_(toneCurve);
            STEP_(rgbCurves);
            STEP_(labAdjustments);
            // stop = stop || STEP_s_(textureBoost);
            STEP_(softLight);
        }
        stop = stop || STEP_s_(localContrast);
        if (!stop) {
            STEP_(filmSimulation);
            STEP_(blackAndWhite);
            STEP_(filmGrain);
        }
        if (pipeline == Pipeline::PREVIEW && params->prsharpening.enabled) {
            double s = scale;
            int fw = full_width * s, fh = full_height * s;
            int imw, imh;
            double s2 = resizeScale(params, fw, fh, imw, imh);
            scale = s * s2;
            STEP_s_(prsharpening);
            scale = s;
        }
        break;
    }
    return stop;
}


int ImProcFunctions::setDeltaEData(EditUniqueID id, double x, double y)
{
    deltaE.ok = false;
    deltaE.x = x;
    deltaE.y = y;
    deltaE.L = 0;
    deltaE.C = 0;
    deltaE.H = 0;

    switch (id) {
    case EUID_LabMasks_DE1:
        return LUMINANCECURVE | M_LUMACURVE;
    case EUID_LabMasks_DE2:
        return DISPLAY;
    case EUID_LabMasks_DE3:
        return LUMINANCECURVE | M_LUMACURVE;
    case EUID_LabMasks_DE4:
        return DISPLAY;
    default:
        return 0;
    }
}


} // namespace rtengine
