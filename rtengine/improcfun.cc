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
#include "EdgePreservingDecomposition.h"
#include "improccoordinator.h"
#include "clutstore.h"
#include "StopWatch.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/guiutils.h"

#undef CLIPD
#define CLIPD(a) ((a)>0.0f?((a)<1.0f?(a):1.0f):0.0f)

namespace {

using namespace rtengine;


// begin of helper function for rgbProc()
void shadowToneCurve(const LUTf &shtonecurve, Imagefloat *rgb, bool multiThread)
{
    float **rtemp = rgb->r.ptrs;
    float **gtemp = rgb->g.ptrs;
    float **btemp = rgb->b.ptrs;
    int W = rgb->getWidth();
    int H = rgb->getHeight();

#if defined( __SSE2__ ) && defined( __x86_64__ )
    vfloat cr = F2V(0.299f);
    vfloat cg = F2V(0.587f);
    vfloat cb = F2V(0.114f);
#endif

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        for (; x < W - 3; x += 4) {
            vfloat rv = LVF(rtemp[y][x]);
            vfloat gv = LVF(gtemp[y][x]);
            vfloat bv = LVF(btemp[y][x]);

            //shadow tone curve
            vfloat Yv = cr * rv + cg * gv + cb * bv;
            vfloat tonefactorv = shtonecurve(Yv);
            STVF(rtemp[y][x], rv * tonefactorv);
            STVF(gtemp[y][x], gv * tonefactorv);
            STVF(btemp[y][x], bv * tonefactorv);
        }
#endif

        for (; x < W; ++x) {
            float r = rtemp[y][x];
            float g = gtemp[y][x];
            float b = btemp[y][x];

            //shadow tone curve
            float Y = (0.299f * r + 0.587f * g + 0.114f * b);
            float tonefactor = shtonecurve[Y];
            rtemp[y][x] = rtemp[y][x] * tonefactor;
            gtemp[y][x] = gtemp[y][x] * tonefactor;
            btemp[y][x] = btemp[y][x] * tonefactor;
        }
    }
}

void highlightToneCurve(const LUTf &hltonecurve, Imagefloat *rgb, float exp_scale, float comp, float hlrange, bool multiThread)
{
    float **rtemp = rgb->r.ptrs;
    float **gtemp = rgb->g.ptrs;
    float **btemp = rgb->b.ptrs;
    int W = rgb->getWidth();
    int H = rgb->getHeight();

#if defined( __SSE2__ ) && defined( __x86_64__ )
    vfloat threev = F2V(3.f);
    vfloat maxvalfv = F2V(MAXVALF);
#endif

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        for (; x < W - 3; x += 4) {

            vfloat rv = LVF(rtemp[y][x]);
            vfloat gv = LVF(gtemp[y][x]);
            vfloat bv = LVF(btemp[y][x]);

            //TODO: proper treatment of out-of-gamut colors
            //float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
            vmask maxMask = vmaskf_ge(vmaxf(rv, vmaxf(gv, bv)), maxvalfv);

            if (_mm_movemask_ps((vfloat)maxMask)) {
                for (int k = 0; k < 4; ++k) {
                    float r = rtemp[y][x + k];
                    float g = gtemp[y][x + k];
                    float b = btemp[y][x + k];
                    float tonefactor = ((r < MAXVALF ? hltonecurve[r] : CurveFactory::hlcurve(exp_scale, comp, hlrange, r)) +
                                        (g < MAXVALF ? hltonecurve[g] : CurveFactory::hlcurve(exp_scale, comp, hlrange, g)) +
                                        (b < MAXVALF ? hltonecurve[b] : CurveFactory::hlcurve(exp_scale, comp, hlrange, b))) / 3.0;

                    // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
                    rtemp[y][x + k] = r * tonefactor;
                    gtemp[y][x + k] = g * tonefactor;
                    btemp[y][x + k] = b * tonefactor;
                }
            } else {
                vfloat tonefactorv = (hltonecurve.cb(rv) + hltonecurve.cb(gv) + hltonecurve.cb(bv)) / threev;
                // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
                STVF(rtemp[y][x], rv * tonefactorv);
                STVF(gtemp[y][x], gv * tonefactorv);
                STVF(btemp[y][x], bv * tonefactorv);
            }
        }

#endif

        for (; x < W; ++x) {
            float r = rtemp[y][x];
            float g = gtemp[y][x];
            float b = btemp[y][x];

            //TODO: proper treatment of out-of-gamut colors
            //float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
            float tonefactor = ((r < MAXVALF ? hltonecurve[r] : CurveFactory::hlcurve(exp_scale, comp, hlrange, r)) +
                                (g < MAXVALF ? hltonecurve[g] : CurveFactory::hlcurve(exp_scale, comp, hlrange, g)) +
                                (b < MAXVALF ? hltonecurve[b] : CurveFactory::hlcurve(exp_scale, comp, hlrange, b))) / 3.0;

            // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
            rtemp[y][x] = r * tonefactor;
            gtemp[y][x] = g * tonefactor;
            btemp[y][x] = b * tonefactor;
        }
    }
}

void proPhotoBlue(float *rtemp, float *gtemp, float *btemp, int istart, int tH, int jstart, int tW, int tileSize) {
    // this is a hack to avoid the blue=>black bug (Issue 2141)
    for (int i = istart, ti = 0; i < tH; i++, ti++) {
        int j = jstart, tj = 0;
#ifdef __SSE2__
        for (; j < tW - 3; j+=4, tj+=4) {
            vfloat rv = LVF(rtemp[ti * tileSize + tj]);
            vfloat gv = LVF(gtemp[ti * tileSize + tj]);
            vmask zeromask = vorm(vmaskf_eq(rv, ZEROV), vmaskf_eq(gv, ZEROV));
            if(_mm_movemask_ps((vfloat)zeromask)) {
                for (int k = 0; k < 4; ++k) {
                    float r = rtemp[ti * tileSize + tj + k];
                    float g = gtemp[ti * tileSize + tj + k];
                    float b = btemp[ti * tileSize + tj + k];
                    
                    if ((r == 0.0f || g == 0.0f) && rtengine::min(r, g, b) >= 0.f) {
                        float h, s, v;
                        Color::rgb2hsv (r, g, b, h, s, v);
                        s *= 0.99f;
                        Color::hsv2rgb (h, s, v, rtemp[ti * tileSize + tj + k], gtemp[ti * tileSize + tj + k], btemp[ti * tileSize + tj + k]);
                    }
                }
            }
        }
#endif
        for (; j < tW; j++, tj++) {
            float r = rtemp[ti * tileSize + tj];
            float g = gtemp[ti * tileSize + tj];
            float b = btemp[ti * tileSize + tj];

            if ((r == 0.0f || g == 0.0f) && rtengine::min(r, g, b) >= 0.f) {
                float h, s, v;
                Color::rgb2hsv (r, g, b, h, s, v);
                s *= 0.99f;
                Color::hsv2rgb (h, s, v, rtemp[ti * tileSize + tj], gtemp[ti * tileSize + tj], btemp[ti * tileSize + tj]);
            }
        }
    }
}


} // namespace

namespace rtengine
{

using namespace procparams;

extern const Settings* settings;

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


void ImProcFunctions::updateColorProfiles (const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, bool gamutCheck)
{
    // set up monitor transform
    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }
    gamutWarning.reset(nullptr);

    monitorTransform = nullptr;

    cmsHPROFILE monitor = nullptr;

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
        cmsHPROFILE iprof  = cmsCreateLab4Profile (nullptr);
        cmsHPROFILE gamutprof = nullptr;
        cmsUInt32Number gamutbpc = 0;
        RenderingIntent gamutintent = RI_RELATIVE;

        bool softProofCreated = false;

        if (softProof) {
            cmsHPROFILE oprof = nullptr;
            RenderingIntent outIntent;
            
            flags = cmsFLAGS_SOFTPROOFING | cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (!settings->printerProfile.empty()) {
                oprof = ICCStore::getInstance()->getProfile (settings->printerProfile);
                if (settings->printerBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }
                outIntent = settings->printerIntent;
            } else {
                oprof = ICCStore::getInstance()->getProfile(params->icm.outputProfile);
                if (params->icm.outputBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }
                outIntent = params->icm.outputIntent;
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
                                       iprof, TYPE_Lab_FLT,
                                       monitor, TYPE_RGB_FLT,
                                       softproof, //oprof,
                                       monitorIntent, outIntent,
                                       flags
                                   );

                if (softproof) {
                    cmsCloseProfile(softproof);
                }

                if (monitorTransform) {
                    softProofCreated = true;
                }

                if (gamutCheck) {
                    gamutprof = oprof;
                    if (params->icm.outputBPC) {
                        gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
                    }
                    gamutintent = outIntent;
                }
            }
        } else if (gamutCheck) {
            // flags = cmsFLAGS_GAMUTCHECK | cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;
            // if (settings->monitorBPC) {
            //     flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            // }

            // monitorTransform = cmsCreateProofingTransform(iprof, TYPE_Lab_FLT, monitor, TYPE_RGB_8, monitor, monitorIntent, monitorIntent, flags);

            // if (monitorTransform) {
            //     softProofCreated = true;
            // }
            gamutprof = monitor;
            if (settings->monitorBPC) {
                gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
            }
            gamutintent = monitorIntent;
        }

        if (!softProofCreated) {
            flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (settings->monitorBPC) {
                flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            }

            monitorTransform = cmsCreateTransform (iprof, TYPE_Lab_FLT, monitor, TYPE_RGB_FLT, monitorIntent, flags);
        }

        if (gamutCheck && gamutprof) {
            gamutWarning.reset(new GamutWarning(iprof, gamutprof, gamutintent, gamutbpc));
        }

        cmsCloseProfile (iprof);
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


void ImProcFunctions::moyeqt (Imagefloat* working, float &moyS, float &eqty)
{
    BENCHFUN

    int tHh = working->getHeight();
    int tWw = working->getWidth();
    double moy = 0.0;
    double sqrs = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:moy,sqrs) schedule(dynamic,16)
#endif

    for (int i = 0; i < tHh; i++) {
        for (int j = 0; j < tWw; j++) {
            float s = Color::rgb2s (CLIP (working->r (i, j)), CLIP (working->g (i, j)), CLIP (working->b (i, j)));
            moy += s;
            sqrs += SQR (s);
        }
    }

    moy /= (tHh * tWw);
    sqrs /= (tHh * tWw);
    eqty = sqrt (sqrs - SQR (moy));
    moyS = moy;
}


// Process RGB image and convert to LAB space
void ImProcFunctions::rgbProc(Imagefloat *working)
{
    BENCHFUN

    working->setMode(Imagefloat::Mode::RGB, multiThread);
        
    constexpr int TS = 112;
    
    LUTf hltonecurve(65536);
    LUTf shtonecurve(65536);
    LUTf rCurve, gCurve, bCurve;

    double expcomp = params->exposure.enabled ? params->exposure.expcomp : 0.0;
    int hlcompr = params->exposure.enabled ? params->exposure.hlcompr : 0;
    int hlcomprthresh = params->exposure.hlcomprthresh;

    {
        int black = params->exposure.enabled ? params->exposure.black : 0;
        int shcompr = params->exposure.enabled ? params->exposure.shcompr : 0;
        
        LUTf tonecurve(65536);
        LUTu vhist16(65536), histToneCurve(256);
        ToneCurve customToneCurve1, customToneCurve2;
        
        CurveFactory::complexCurve(expcomp, black / 65535.0,
                                   hlcompr, hlcomprthresh,
                                   shcompr, 0, 0, 
                                   { DCT_Linear }, { DCT_Linear },
                                   vhist16, hltonecurve, shtonecurve, tonecurve,
                                   histToneCurve, customToneCurve1,
                                   customToneCurve2, scale);
    }

    // std::unique_ptr<Imagefloat> workimage(working->copy());
    // working = workimage.get();
    
    Imagefloat *tmpImage = nullptr;

    Imagefloat* editImgFloat = nullptr;
    PlanarWhateverData<float>* editWhatever = nullptr;
    EditUniqueID editID = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if (editID != EUID_None) {
        switch  (pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType()) {
            case (BT_IMAGEFLOAT):
                editImgFloat = pipetteBuffer->getImgFloatBuffer();
                break;

            case (BT_LABIMAGE):
                break;

            case (BT_SINGLEPLANE_FLOAT):
                editWhatever = pipetteBuffer->getSinglePlaneBuffer();
                break;
        }
    }

    if (params->rgbCurves.enabled) {
        CurveFactory::RGBCurve(params->rgbCurves.rcurve, rCurve, scale);
        CurveFactory::RGBCurve(params->rgbCurves.gcurve, gCurve, scale);
        CurveFactory::RGBCurve(params->rgbCurves.bcurve, bCurve, scale);
    }

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params->icm.workingProfile);
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.workingProfile);

    float toxyz[3][3] = {
        {
            static_cast<float> ( wprof[0][0] / Color::D50x),
            static_cast<float> ( wprof[0][1] / Color::D50x),
            static_cast<float> ( wprof[0][2] / Color::D50x)
        }, {
            static_cast<float> ( wprof[1][0]),
            static_cast<float> ( wprof[1][1]),
            static_cast<float> ( wprof[1][2])
        }, {
            static_cast<float> ( wprof[2][0] / Color::D50z),
            static_cast<float> ( wprof[2][1] / Color::D50z),
            static_cast<float> ( wprof[2][2] / Color::D50z)
        }
    };
    float maxFactorToxyz = max (toxyz[1][0], toxyz[1][1], toxyz[1][2]);
    float equalR = maxFactorToxyz / toxyz[1][0];
    float equalG = maxFactorToxyz / toxyz[1][1];
    float equalB = maxFactorToxyz / toxyz[1][2];

    //inverse matrix user select
    double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };

    double wp[3][3] = {
        {wprof[0][0], wprof[0][1], wprof[0][2]},
        {wprof[1][0], wprof[1][1], wprof[1][2]},
        {wprof[2][0], wprof[2][1], wprof[2][2]}
    };

    bool mixchannels = params->chmixer.enabled &&
        (params->chmixer.red[0] != 100 || params->chmixer.red[1] != 0     || params->chmixer.red[2] != 0   ||
                        params->chmixer.green[0] != 0 || params->chmixer.green[1] != 100 || params->chmixer.green[2] != 0 ||
                        params->chmixer.blue[0] != 0  || params->chmixer.blue[1] != 0    || params->chmixer.blue[2] != 100);

    FlatCurve* bwlCurve = nullptr;

    FlatCurveType bwlCurveType = (FlatCurveType)params->blackwhite.luminanceCurve.at (0);
    bool bwlCurveEnabled = bwlCurveType > FCT_Linear;

    if (bwlCurveEnabled) {
        bwlCurve = new FlatCurve (params->blackwhite.luminanceCurve);

        if (bwlCurve->isIdentity()) {
            delete bwlCurve;
            bwlCurve = nullptr;
            bwlCurveEnabled = false;
        }
    }

    const float exp_scale = pow (2.0, expcomp);
    const float comp = (max (0.0, expcomp) + 1.0) * hlcompr / 100.0;
    const float shoulder = ((65536.0 / max (1.0f, exp_scale)) * (hlcomprthresh / 200.0)) + 0.1;
    const float hlrange = 65536.0 - shoulder;
    const bool isProPhoto = (params->icm.workingProfile == "ProPhoto");
    bool highlight = params->exposure.enabled && params->exposure.hrmode != procparams::ExposureParams::HR_OFF;

    float chMixRR = float (params->chmixer.red[0])/10.f;
    float chMixRG = float (params->chmixer.red[1])/10.f;
    float chMixRB = float (params->chmixer.red[2])/10.f;
    float chMixGR = float (params->chmixer.green[0])/10.f;
    float chMixGG = float (params->chmixer.green[1])/10.f;
    float chMixGB = float (params->chmixer.green[2])/10.f;
    float chMixBR = float (params->chmixer.blue[0])/10.f;
    float chMixBG = float (params->chmixer.blue[1])/10.f;
    float chMixBB = float (params->chmixer.blue[2])/10.f;

    bool blackwhite = params->blackwhite.enabled;
    float bwr = float (params->blackwhite.mixerRed);
    float bwg = float (params->blackwhite.mixerGreen);
    float bwb = float (params->blackwhite.mixerBlue);
    float bwrgam = float (params->blackwhite.gammaRed);
    float bwggam = float (params->blackwhite.gammaGreen);
    float bwbgam = float (params->blackwhite.gammaBlue);
    int algm = 0;

    if     (params->blackwhite.method == "Desaturation") {
        algm = 0;
    } else if (params->blackwhite.method == "LumEqualizer") {
        algm = 1;
    } else if (params->blackwhite.method == "ChannelMixer") {
        algm = 2;
    }

    float kcorec = 1.f;
    //gamma correction of each channel
    float gamvalr = 125.f;
    float gamvalg = 125.f;
    float gamvalb = 125.f;

    if (bwrgam < 0) {
        gamvalr = 100.f;
    }

    if (bwggam < 0) {
        gamvalg = 100.f;
    }

    if (bwbgam < 0) {
        gamvalb = 100.f;
    }

    float gammabwr = 1.f;
    float gammabwg = 1.f;
    float gammabwb = 1.f;
    //if     (params->blackwhite.setting=="Ma" || params->blackwhite.setting=="Mr" || params->blackwhite.setting=="Fr" || params->blackwhite.setting=="Fa")  {
    {
        gammabwr = 1.f - bwrgam / gamvalr;
        gammabwg = 1.f - bwggam / gamvalg;
        gammabwb = 1.f - bwbgam / gamvalb;
    }
    bool hasgammabw = gammabwr != 1.f || gammabwg != 1.f || gammabwb != 1.f;

    if (blackwhite) {// || (params->dirpyrequalizer.cbdlMethod == "bef" && params->dirpyrequalizer.enabled)) {
        tmpImage = new Imagefloat (working->getWidth(), working->getHeight());
    }

    if (mixchannels) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < working->getHeight(); ++y) {
            for (int x = 0; x < working->getWidth(); ++x) {
                float r = working->r(y, x);
                float g = working->g(y, x);
                float b = working->b(y, x);

                float rmix = (r * chMixRR + g * chMixRG + b * chMixRB) / 100.f;
                float gmix = (r * chMixGR + g * chMixGG + b * chMixGB) / 100.f;
                float bmix = (r * chMixBR + g * chMixBG + b * chMixBB) / 100.f;

                working->r(y, x) = rmix;
                working->g(y, x) = gmix;
                working->b(y, x) = bmix;
            }
        }
    }

    highlightToneCurve(hltonecurve, working, exp_scale, comp, hlrange, multiThread);
    shadowToneCurve(shtonecurve, working, multiThread);
    toneEqualizer(working);

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        size_t perChannelSizeBytes = padToAlignment(sizeof (float) * TS * TS + 4 * 64);
        AlignedBuffer<float> buffer(3 * perChannelSizeBytes);
        char *editIFloatBuffer = nullptr;
        char *editWhateverBuffer = nullptr;

        float *rtemp = buffer.data;
        float *gtemp = &rtemp[perChannelSizeBytes / sizeof(float)];
        float *btemp = &gtemp[perChannelSizeBytes / sizeof(float)];
        int istart;
        int jstart;
        int tW;
        int tH;

        // zero out the buffers
        memset(rtemp, 0, 3 * perChannelSizeBytes);

        float *editIFloatTmpR = nullptr, *editIFloatTmpG = nullptr, *editIFloatTmpB = nullptr, *editWhateverTmp = nullptr;

        if (editImgFloat) {
            editIFloatBuffer = (char *) malloc (3 * sizeof (float) * TS * TS + 20 * 64 + 63);
            char *data = (char*) ( ( uintptr_t (editIFloatBuffer) + uintptr_t (63)) / 64 * 64);

            editIFloatTmpR = (float (*))data;
            editIFloatTmpG = (float (*))         ((char*)editIFloatTmpR + sizeof (float) * TS * TS + 4 * 64);
            editIFloatTmpB = (float (*))         ((char*)editIFloatTmpG + sizeof (float) * TS * TS + 8 * 64);
        }

        if (editWhatever) {
            editWhateverBuffer = (char *) malloc (sizeof (float) * TS * TS + 20 * 64 + 63);
            char *data = (char*) ( ( uintptr_t (editWhateverBuffer) + uintptr_t (63)) / 64 * 64);

            editWhateverTmp = (float (*))data;
        }

        // float out_rgbx[4 * TS] ALIGNED16; // Line buffer for CLUT
        // float clutr[TS] ALIGNED16;
        // float clutg[TS] ALIGNED16;
        // float clutb[TS] ALIGNED16;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif
        for (int ii = 0; ii < working->getHeight(); ii += TS) {
            for (int jj = 0; jj < working->getWidth(); jj += TS) {
                istart = ii;
                jstart = jj;
                tH = min (ii + TS, working->getHeight());
                tW = min (jj + TS, working->getWidth());


                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                        rtemp[ti * TS + tj] = working->r (i, j);
                        gtemp[ti * TS + tj] = working->g (i, j);
                        btemp[ti * TS + tj] = working->b (i, j);
                    }
                }
                
                if (editID == EUID_RGB_R) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[rtemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                } else if (editID == EUID_RGB_G) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[gtemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                } else if (editID == EUID_RGB_B) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[btemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                }

                if (params->rgbCurves.enabled && (rCurve || gCurve || bCurve)) { // if any of the RGB curves is engaged
                    if (!params->rgbCurves.lumamode) { // normal RGB mode

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                // individual R tone curve
                                if (rCurve) {
                                    setUnlessOOG(rtemp[ti * TS + tj], rCurve[ rtemp[ti * TS + tj] ]);
                                }

                                // individual G tone curve
                                if (gCurve) {
                                    setUnlessOOG(gtemp[ti * TS + tj], gCurve[ gtemp[ti * TS + tj] ]);
                                }

                                // individual B tone curve
                                if (bCurve) {
                                    setUnlessOOG(btemp[ti * TS + tj], bCurve[ btemp[ti * TS + tj] ]);
                                }
                            }
                        }
                    } else { //params->rgbCurves.lumamode==true (Luminosity mode)
                        // rCurve.dump("r_curve");//debug

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                // rgb values before RGB curves
                                float r = rtemp[ti * TS + tj] ;
                                float g = gtemp[ti * TS + tj] ;
                                float b = btemp[ti * TS + tj] ;
                                //convert to Lab to get a&b before RGB curves
                                float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
                                float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                                float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;

                                float fx = x < MAXVALF ? Color::cachef[x] : 327.68f * std::cbrt (x / MAXVALF);
                                float fy = y < MAXVALF ? Color::cachef[y] : 327.68f * std::cbrt (y / MAXVALF);
                                float fz = z < MAXVALF ? Color::cachef[z] : 327.68f * std::cbrt (z / MAXVALF);

                                float a_1 = 500.0f * (fx - fy);
                                float b_1 = 200.0f * (fy - fz);

                                // rgb values after RGB curves
                                if (rCurve) {
                                    float rNew = rCurve[r];
                                    r += (rNew - r) * equalR;
                                }

                                if (gCurve) {
                                    float gNew = gCurve[g];
                                    g += (gNew - g) * equalG;
                                }

                                if (bCurve) {
                                    float bNew = bCurve[b];
                                    b += (bNew - b) * equalB;
                                }

                                // Luminosity after
                                // only Luminance in Lab
                                float newy = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                                float L_2 = newy <= MAXVALF ? Color::cachefy[newy] : 327.68f * (116.f * xcbrtf(newy / MAXVALF) - 16.f);

                                //gamut control
                                if (settings->rgbcurveslumamode_gamut) {
                                    float Lpro = L_2 / 327.68f;
                                    float Chpro = sqrtf (SQR (a_1) + SQR (b_1)) / 327.68f;
                                    float HH = NAN; // we set HH to NAN, because then it will be calculated in Color::gamutLchonly only if needed
//                                    float HH = xatan2f(b_1, a_1);
                                    // According to mathematical laws we can get the sin and cos of HH by simple operations even if we don't calculate HH
                                    float2 sincosval;

                                    if (Chpro == 0.0f) {
                                        sincosval.y = 1.0f;
                                        sincosval.x = 0.0f;
                                    } else {
                                        sincosval.y = a_1 / (Chpro * 327.68f);
                                        sincosval.x = b_1 / (Chpro * 327.68f);
                                    }

#ifdef _DEBUG
                                    bool neg = false;
                                    bool more_rgb = false;
                                    //gamut control : Lab values are in gamut
                                    Color::gamutLchonly (HH, sincosval, Lpro, Chpro, r, g, b, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                    //gamut control : Lab values are in gamut
                                    Color::gamutLchonly (HH, sincosval, Lpro, Chpro, r, g, b, wip, highlight, 0.15f, 0.96f);
#endif
                                    //end of gamut control
                                } else {
                                    float x_, y_, z_;
                                    //calculate RGB with L_2 and old value of a and b
                                    Color::Lab2XYZ (L_2, a_1, b_1, x_, y_, z_) ;
                                    Color::xyz2rgb (x_, y_, z_, r, g, b, wip);
                                }

                                setUnlessOOG(rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], r, g, b);
                            }
                        }
                    }
                }

                if (editID == EUID_HSV_H || editID == EUID_HSV_S || editID == EUID_HSV_V) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float h, s, v;
                            Color::rgb2hsv (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], h, s, v);
                            editWhateverTmp[ti * TS + tj] = h;
                        }
                    }
                }

                if (isProPhoto) { // this is a hack to avoid the blue=>black bug (Issue 2141)
                    proPhotoBlue(rtemp, gtemp, btemp, istart, tH, jstart, tW, TS);
                }

                // filling the pipette buffer
                if (editID == EUID_BlackWhiteLuminance) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float X, Y, Z, L, aa, bb;
                            //rgb=>lab
                            Color::rgbxyz (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], X, Y, Z, wp);
                            //convert Lab
                            Color::XYZ2Lab (X, Y, Z, L, aa, bb);
                            //end rgb=>lab
                            float HH = xatan2f (bb, aa); // HH hue in -3.141  +3.141

                            editWhateverTmp[ti * TS + tj] = float (Color::huelab_to_huehsv2 (HH));
                        }
                    }
                }

                //black and white
                if (blackwhite) {
                    if (algm == 0) { //lightness
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                                float r = rtemp[ti * TS + tj];
                                float g = gtemp[ti * TS + tj];
                                float b = btemp[ti * TS + tj];

                                // --------------------------------------------------

                                // Method 1: Luminosity (code taken from Gimp)
                                /*
                                float maxi = max(r, g, b);
                                float mini = min(r, g, b);
                                r = g = b = (maxi+mini)/2;
                                */

                                // Method 2: Luminance (former RT code)
                                r = g = b = (0.299f * r + 0.587f * g + 0.114f * b);

                                // --------------------------------------------------

#ifndef __SSE2__

                                //gamma correction: pseudo TRC curve
                                if (hasgammabw) {
                                    Color::trcGammaBW (r, g, b, gammabwr, gammabwg, gammabwb);
                                }

#endif
                                rtemp[ti * TS + tj] = r;
                                gtemp[ti * TS + tj] = g;
                                btemp[ti * TS + tj] = b;
                            }

#ifdef __SSE2__

                            if (hasgammabw) {
                                //gamma correction: pseudo TRC curve
                                Color::trcGammaBWRow (&rtemp[ti * TS], &gtemp[ti * TS], &btemp[ti * TS], tW - jstart, gammabwr, gammabwg, gammabwb);
                            }

#endif

                        }
                    } else if (algm == 1) { //Luminance mixer in Lab mode to avoid artifacts
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                //rgb => xyz
                                float X, Y, Z;
                                Color::rgbxyz (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], X, Y, Z, wp);
                                //xyz => Lab
                                float L, aa, bb;
                                Color::XYZ2Lab (X, Y, Z, L, aa, bb);
                                float CC = sqrtf (SQR (aa) + SQR (bb)) / 327.68f; //CC chromaticity in 0..180 or more
                                float HH = xatan2f (bb, aa); // HH hue in -3.141  +3.141
                                float2 sincosval;

                                if (CC == 0.0f) {
                                    sincosval.y = 1.f;
                                    sincosval.x = 0.0f;
                                } else {
                                    sincosval.y = aa / (CC * 327.68f);
                                    sincosval.x = bb / (CC * 327.68f);
                                }

                                if (bwlCurveEnabled) {
                                    L /= 32768.f;
                                    double hr = Color::huelab_to_huehsv2 (HH);
                                    float valparam = float ((bwlCurve->getVal (hr) - 0.5f) * 2.0f); //get l_r=f(H)
                                    float kcc = (CC / 70.f); //take Chroma into account...70 "middle" of chromaticity (arbitrary and simple), one can imagine other algorithme
                                    //reduct action for low chroma and increase action for high chroma
                                    valparam *= kcc;

                                    if (valparam > 0.f) {
                                        L = (1.f - valparam) * L + valparam * (1.f - SQR (SQR (SQR (SQR (1.f - min (L, 1.0f)))))); // SQR (SQR((SQR)  to increase action in low light
                                    } else {
                                        L *= (1.f + valparam);    //for negative
                                    }

                                    L *= 32768.f;
                                }

                                float RR, GG, BB;
                                L /= 327.68f;
#ifdef _DEBUG
                                bool neg = false;
                                bool more_rgb = false;
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (HH, sincosval, L, CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (HH, sincosval, L, CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f);
#endif
                                L *= 327.68f;
                                //convert l => rgb
                                Color::L2XYZ (L, X, Y, Z);
                                float newRed; // We use the red channel for bw
                                Color::xyz2r (X, Y, Z, newRed, wip);
                                rtemp[ti * TS + tj] = gtemp[ti * TS + tj] = btemp[ti * TS + tj] = newRed;
#ifndef __SSE2__

                                if (hasgammabw) {
                                    //gamma correction: pseudo TRC curve
                                    Color::trcGammaBW (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], gammabwr, gammabwg, gammabwb);
                                }

#endif
                            }

#ifdef __SSE2__

                            if (hasgammabw) {
                                //gamma correction: pseudo TRC curve
                                Color::trcGammaBWRow (&rtemp[ti * TS], &gtemp[ti * TS], &btemp[ti * TS], tW - jstart, gammabwr, gammabwg, gammabwb);
                            }

#endif
                        }
                    }
                }

                if (!blackwhite) {
                    if (editImgFloat || editWhatever) {
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                                // filling the pipette buffer by the content of the temp pipette buffers
                                if (editImgFloat) {
                                    editImgFloat->r (i, j) = editIFloatTmpR[ti * TS + tj];
                                    editImgFloat->g (i, j) = editIFloatTmpG[ti * TS + tj];
                                    editImgFloat->b (i, j) = editIFloatTmpB[ti * TS + tj];
                                } else if (editWhatever) {
                                    editWhatever->v (i, j) = editWhateverTmp[ti * TS + tj];
                                }
                            }
                        }
                    }
                    // ready, fill lab
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            int idx = ti * TS + tj;
                            working->r(i, j) = rtemp[idx];
                            working->g(i, j) = gtemp[idx];
                            working->b(i, j) = btemp[idx];
                        }
                    }
                } else { // black & white
                    // Auto channel mixer needs whole image, so we now copy to tmpImage and close the tiled processing
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            // filling the pipette buffer by the content of the temp pipette buffers
                            if (editImgFloat) {
                                editImgFloat->r (i, j) = editIFloatTmpR[ti * TS + tj];
                                editImgFloat->g (i, j) = editIFloatTmpG[ti * TS + tj];
                                editImgFloat->b (i, j) = editIFloatTmpB[ti * TS + tj];
                            } else if (editWhatever) {
                                editWhatever->v (i, j) = editWhateverTmp[ti * TS + tj];
                            }

                            tmpImage->r (i, j) = rtemp[ti * TS + tj];
                            tmpImage->g (i, j) = gtemp[ti * TS + tj];
                            tmpImage->b (i, j) = btemp[ti * TS + tj];
                        }
                    }
                }
            }
        }

        if (editIFloatBuffer) {
            free (editIFloatBuffer);
        }

        if (editWhateverBuffer) {
            free (editWhateverBuffer);
        }
    }

    // starting a new tile processing with a 'reduction' clause for the auto mixer computing
    if (blackwhite) {//channel-mixer
        int tW = working->getWidth();
        int tH = working->getHeight();

        if (algm == 2) { //channel-mixer
            float filcor;
            double rrm, ggm, bbm;
            Color::computeBWMixerConstants (params->blackwhite.setting, params->blackwhite.filter, "", filcor,
                                            bwr, bwg, bwb, kcorec, rrm, ggm, bbm);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 16)
#endif

            for (int i = 0; i < tH; i++) {
                for (int j = 0; j < tW; j++) {

                    //mix channel
                    tmpImage->r (i, j) = tmpImage->g (i, j) = tmpImage->b (i, j) = /*CLIP*/ ((bwr * tmpImage->r (i, j) + bwg * tmpImage->g (i, j) + bwb * tmpImage->b (i, j)) * kcorec);

#ifndef __SSE2__

                    //gamma correction: pseudo TRC curve
                    if (hasgammabw) {
                        Color::trcGammaBW (tmpImage->r (i, j), tmpImage->g (i, j), tmpImage->b (i, j), gammabwr, gammabwg, gammabwb);
                    }

#endif
                }

#ifdef __SSE2__

                if (hasgammabw) {
                    //gamma correction: pseudo TRC curve
                    Color::trcGammaBWRow (tmpImage->r (i), tmpImage->g (i), tmpImage->b (i), tW, gammabwr, gammabwg, gammabwb);
                }

#endif
            }
        }

        // ready, fill lab (has to be the same code than the "fill lab" above!)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 5)
#endif

        for (int i = 0; i < tH; i++) {
            for (int j = 0; j < tW; j++) {
                working->r(i, j) = tmpImage->r(i, j);
                working->g(i, j) = tmpImage->g(i, j);
                working->b(i, j) = tmpImage->b(i, j);
            }
        }
    }

    if (tmpImage) {
        delete tmpImage;
    }
}

/**
* @brief retreave RGB value with maximum saturation
* @param r red input and in exit new r
* @param g green input and in exit new g
* @param b blue input and in exit new b
**/
void ImProcFunctions::retreavergb (float &r, float &g, float &b)
{
    float mini = min (r, g, b);
    float maxi = max (r, g, b);
    float kkm = 65535.f / maxi;

    if (b == mini && r == maxi) {
        r = 65535.f;
        g = kkm * (g - b);
        b = 0.f;
    } else if (b == mini && g == maxi) {
        g = 65535.f;
        r = kkm * (r - b);
        b = 0.f;
    } else if (g == mini && r == maxi) {
        r = 65535.f;
        b = kkm * (b - g);
        g = 0.f;
    } else if (g == mini && b == maxi) {
        b = 65535.f;
        r = kkm * (r - g);
        g = 0.f;
    } else if (r == mini && b == maxi) {
        b = 65535.f;
        g = kkm * (g - r);
        r = 0.f;
    } else if (r == mini && g == maxi) {
        g = 65535.f;
        b = kkm * (b - r);
        r = 0.f;
    }
}

/**
* @brief Interpolate by decreasing with a parabol k = aa*v*v + bb*v +c  v[0..1]
* @param reducac val ue of the reduction in the middle of the range
* @param vinf value [0..1] for beginning decrease
* @param aa second degree parameter
* @param bb first degree parameter
* @param cc third parameter
**/
void ImProcFunctions::secondeg_end (float reducac, float vinf, float &aa, float &bb, float &cc)
{
    float zrd = reducac; //value at me  linear =0.5
    float v0 = vinf; //max shadows
    float me = (1.f + v0) / 2.f; //"median" value = (v0 + 1.=/2)
    //float a1=1.f-v0;
    float a2 = me - v0;
    float a3 = 1.f - v0 * v0;
    float a4 = me * me - v0 * v0;
    aa = (1.f + (zrd - 1.f) * (1 - v0) / a2) / (a4 * (1.f - v0) / a2 - a3);
    bb = - (1.f + a3 * aa) / (1.f - v0);
    cc = - (aa + bb);
}

/**
* @brief Interpolate by increasing with a parabol k = aa*v*v + bb*v  v[0..1]
* @param reducac val ue of the reduction in the middle of the range
* @param vend value [0..1] for beginning increase
* @param aa second degree parameter
* @param bb first degree parameter
**/
void ImProcFunctions::secondeg_begin (float reducac, float vend, float &aam, float &bbm)
{
    aam = (2.f - 4.f * reducac) / (vend * vend);
    bbm = 1.f / vend - aam * vend;
}


/**
* @brief color toning with 9 sliders shadows middletones highlight
* @param r red input values [0..65535]
* @param g green input values [0..65535]
* @param b blue input values [0..65535]
* @param ro red output values [0..65535]
* @param go green output values [0..65535]
* @param bo blue output values [0..65535]
* @param RedLow    [-1..1] value after transformations of sliders [-100..100] for shadows
* @param GreenLow  [-1..1] value after transformations of sliders [-100..100] for shadows
* @param BlueLow   [-1..1] value after transformations of sliders [-100..100] for shadows
* @param RedMed    [-1..1] value after transformations of sliders [-100..100] for midtones
* @param GreenMed  [-1..1] value after transformations of sliders [-100..100] for midtones
* @param BlueMed   [-1..1] value after transformations of sliders [-100..100] for midtones
* @param RedHigh   [-1..1] value after transformations of sliders [-100..100] for highlights
* @param GreenHigh [-1..1] value after transformations of sliders [-100..100] for highlights
* @param BlueHigh  [-1..1] value after transformations of sliders [-100..100] for highlights
* @param reducac value of the reduction in the middle of the range for second degree increse or decrease action
* @param mode 0 = colour, 1 = Black and White
* @param strProtect ?
**/
void ImProcFunctions::toningsmh(float r, float g, float b, float &ro, float &go, float &bo, float RedLow, float GreenLow, float BlueLow, float RedMed, float GreenMed, float BlueMed, float RedHigh, float GreenHigh, float BlueHigh, float reducac, int mode, float strProtect)
{
    const float v = max(r, g, b) / 65535.f;
    float kl = 1.f;
    float rlo; //0.4  0.5
    float rlm; //1.1
    float rlh; //1.1
    float rlob; //for BW old mode

    if (mode == 0) { //colour
        rlo = rlob = strProtect; //0.5 ==>  0.75
        rlh = 2.2f * strProtect;
        rlm = 1.5f * strProtect;
        constexpr float v0 = 0.15f;
        //second degree

        if (v > v0) {
            float aa, bb, cc;
            secondeg_end (reducac, v0, aa, bb, cc);
            kl = aa * v * v + bb * v + cc;    //verified ==> exact
        } else {
            float aab, bbb;
            secondeg_begin (0.7f, v0, aab, bbb);
            kl = aab * v * v + bbb * v;
        }
    } else { //bw coefficient to preserve same results as before for satlimtopacity = 0.5 (default)
        rlo = strProtect * 0.8f; //0.4
        rlob = strProtect; //0.5
        rlm = strProtect * 2.2f; //1.1
        rlh = strProtect * 2.4f; //1.2
        if (v > 0.15f) {
            kl = (-1.f / 0.85f) * v + 1.f / 0.85f;    //Low light ==> decrease action after v=0.15
        }
    }

    {
        const float corr = 20000.f * RedLow * kl * rlo;
        if (RedLow > 0.f) {
            g -= corr;
            b -= corr;
        } else {
            r += corr;
        }

        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    {
        const float corr = 20000.f * GreenLow * kl * rlo;
        if (GreenLow > 0.f) {
            r -= corr;
            b -= corr;
        } else {
            g += corr;
        }

        // r = CLIP(r);
        // b = CLIP(b);
        // g = CLIP(g);
    }


    {
        const float corr = 20000.f * BlueLow * kl * rlob;

        if (BlueLow > 0.f) {
            r -= corr;
            g -= corr;
        } else {
            b += corr;
        }

        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    // mid tones
    float km;
    constexpr float v0m = 0.5f; //max action

    if (v < v0m) {
        float aam, bbm;
        float vend = v0m;
        secondeg_begin (reducac, vend, aam, bbm);
        km = aam * v * v + bbm * v; //verification = good
    } else {
        float v0mm = 0.5f; //max
        float aamm, bbmm, ccmm;
        secondeg_end (reducac, v0mm, aamm, bbmm, ccmm);
        km = aamm * v * v + bbmm * v + ccmm; //verification good
    }

    {
        const float RedM = RedMed * km * rlm;

        if (RedMed > 0.f) {
            r += 20000.f * RedM;
            g -= 10000.f * RedM;
            b -= 10000.f * RedM;
        } else {
            r += 10000.f * RedM;
            g -= 20000.f * RedM;
            b -= 20000.f * RedM;
        }
        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    {
        const float GreenM = GreenMed * km * rlm;

        if (GreenMed > 0.f) {
            r -= 10000.f * GreenM;
            g += 20000.f * GreenM;
            b -= 10000.f * GreenM;
        } else {
            r -= 20000.f * GreenM;
            g += 10000.f * GreenM;
            b -= 20000.f * GreenM;
        }
        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    {
        const float BlueM = BlueMed * km * rlm;

        if (BlueMed > 0.f) {
            r -= 10000.f * BlueM;
            g -= 10000.f * BlueM;
            b += 20000.f * BlueM;
        } else {
            r -= 20000.f * BlueM;
            g -= 20000.f * BlueM;
            b += 10000.f * BlueM;
        }
        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    //high tones
    constexpr float v00 = 0.8f; //max action
    float aa0, bb0;
    secondeg_begin (reducac, v00, aa0, bb0);

    float kh;
    if (v > v00) { //max action
        kh = (1.f - v) / (1.f - v00);    //High tones
    } else {
        kh = v * (aa0 * v + bb0);    //verification = good
    }

    {
        const float corr = 20000.f * RedHigh * kh * rlh; //1.2

        if (RedHigh > 0.f) {
            r += corr;
        } else {
            g -= corr;
            b -= corr;
        }

        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    {
        const float corr = 20000.f * GreenHigh * kh * rlh; //1.2

        if (GreenHigh > 0.f) {
            g += corr;
        } else {
            r -= corr;
            b -= corr;
        }

        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    {
        const float corr = 20000.f * BlueHigh * kh * rlh; //1.2

        if (BlueHigh > 0.f) {
            b += corr;
        } else {
            r -= corr;
            g -= corr;
        }

        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    ro = r;
    go = g;
    bo = b;
}

/**
* @brief color toning with 2 colors - 2 sliders saturation shadows and highlight and one balance
* @param r g b input values [0..65535]
* @param ro go bo output values [0..65535]
* @param iplow iphigh [0..1] from curve color - value of luminance shadows and highlights
* @param rl gl bl [0..65535] - color of reference shadow
* @param rh gh bh [0..65535] - color of reference highlight
* @param SatLow SatHigh [0..1] from sliders saturation shadows and highlight
* @param balanS [0..1] balance for shadows (one slider)
* @param balanH [0..1] balance for highlights (same slider than for balanS)
* @param reducac value of the reduction in the middle of the range for second degree, increase or decrease action
**/
void ImProcFunctions::toning2col (float r, float g, float b, float &ro, float &go, float &bo, float iplow, float iphigh, float krl, float kgl, float kbl, float krh, float kgh, float kbh, float SatLow, float SatHigh, float balanS, float balanH, float reducac, int mode, int preser, float strProtect)
{
    const float lumbefore = 0.299f * r + 0.587f * g + 0.114f * b;
    const float v = max(r, g, b) / 65535.f;

    const float rlo = strProtect;  //0.5 ==> 0.75  transferred value for more action
    const float rlh = 2.2f * strProtect;

    //low tones
    //second degree
    float aa, bb, cc;
    //fixed value of reducac =0.4;
    secondeg_end (reducac, iplow, aa, bb, cc);

    float aab, bbb;
    secondeg_begin (0.7f, iplow, aab, bbb);

    if (SatLow > 0.f) {
        float kl = 1.f;
        if (v > iplow) {
            kl = aa * v * v + bb * v + cc;
        } else if (mode == 0) {
            kl = aab * v * v + bbb * v;
        }
        const float kmgb = min(r, g, b);
        if (kmgb < 20000.f) {
            //I have tested ...0.85 compromise...
            kl *= pow_F ((kmgb / 20000.f), 0.85f);
        }

        const float factor = 20000.f * SatLow * kl * rlo * balanS;

        if (krl > 0.f) {
            g -= factor * krl;
            b -= factor * krl;
        }

        // g = CLIP(g);
        // b = CLIP(b);

        if (kgl > 0.f) {
            r -= factor * kgl;
            b -= factor * kgl;
        }

        // r = CLIP(r);
        // b = CLIP(b);

        if (kbl > 0.f) {
            r -= factor * kbl;
            g -= factor * kbl;
        }

        // r = CLIP(r);
        // g = CLIP(g);
    }

    //high tones
    float aa0, bb0;
    //fixed value of reducac ==0.4;
    secondeg_begin (reducac, iphigh, aa0, bb0);

    if (SatHigh > 0.f) {
        float kh = 1.f;
        if (v > iphigh) {
            kh = (1.f - v) / (1.f - iphigh);    //Low light ==> decrease action after iplow
        } else {
            kh = aa0 * v * v + bb0 * v;
        }

        const float kmgb = max(r, g, b);
        if (kmgb > 45535.f) {
            constexpr float cora = 1.f / (45535.f - 65535.f);
            constexpr float corb = 1.f - cora * 45535.f;
            kh *= kmgb * cora + corb;
        }
        const float factor = 20000.f * SatHigh * kh * rlh * balanH;
        r += factor * (krh > 0.f ? krh : 0.f);
        g += factor * (kgh > 0.f ? kgh : 0.f);
        b += factor * (kbh > 0.f ? kbh : 0.f);

        // r = CLIP(r);
        // g = CLIP(g);
        // b = CLIP(b);
    }

    float preserv = 1.f;
    if (preser == 1) {
        float lumafter = 0.299f * r + 0.587f * g + 0.114f * b;
        preserv = lumbefore / lumafter;
    }

    setUnlessOOG(ro, go, bo, CLIP(r * preserv), CLIP(g * preserv), CLIP(b * preserv));
}

/**
* @brief color toning with interpolation in mode Lab
* @param r g b input values [0..65535]
* @param ro go bo output values [0..65535]
* @param algm  metchrom twoc - methods
* @param ctColorCurve curve 500 colors
* @param ctOpacityCurve curve standard 'ab'
* @param clToningcurve  curve special 'ab' and 'a'
* @param cl2Toningcurve curve special 'b'
* @param iplow iphigh [0..1] luminance
* @param wp wip 3x3 matrix and inverse conversion rgb XYZ
**/
void ImProcFunctions::labtoning (float r, float g, float b, float &ro, float &go, float &bo, int algm, int metchrom, int twoc, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, LUTf & clToningcurve, LUTf & cl2Toningcurve, float iplow, float iphigh, double wp[3][3], double wip[3][3]  )
{
    ro = CLIP(r);
    go = CLIP(g);
    bo = CLIP(b);
    
    float realL;
    float h, s, l;
    Color::rgb2hsl (ro, go, bo, h, s, l);
    float x2, y2, z2;
    float xl, yl, zl;

    if (twoc != 1) {
        l = (Color::gammatab_13_2[     l * 65535.f]) / 65535.f; //to compensate L from Lab
        iphigh = (Color::gammatab_13_2[iphigh * 65535.f]) / 65535.f;
        iplow  = (Color::gammatab_13_2[ iplow * 65535.f]) / 65535.f;
    }

    if (twoc == 1) {
        ctColorCurve.getVal (l, x2, y2, z2);
    } else {
        ctColorCurve.getVal (iphigh, x2, y2, z2);
        ctColorCurve.getVal (iplow, xl, yl, zl);
    }

    realL = l;


    //float opacity = ctOpacityCurve.lutOpacityCurve[l*500.f];
    //if(params->blackwhite.enabled){satLimit=80.f;satLimitOpacity=30.f;}//force BW

    // get the opacity and tweak it to preserve saturated colors
    //float l_ = Color::gamma_srgb(l*65535.f)/65535.f;
    float opacity = (1.f - min<float> (s / satLimit, 1.f) * (1.f - satLimitOpacity)) * ctOpacityCurve.lutOpacityCurve[l * 500.f];
    float opacity2 = (1.f - min<float> (s / satLimit, 1.f) * (1.f - satLimitOpacity));

    l *= 65535.f;
    float chromat = 0.f, luma = 0.f;

    if (clToningcurve[l] < l) {
        chromat = clToningcurve[l] / l - 1.f;  //special effect
    } else if (clToningcurve[l] > l) {
        chromat = 1.f - SQR(SQR(l / clToningcurve[l])); //apply C=f(L) acts  on 'a' and 'b'
    }

    if (cl2Toningcurve[l] < l) {
        luma = cl2Toningcurve[l] / l - 1.f;  //special effect
    } else if (cl2Toningcurve[l] > l) {
        luma = 1.f - SQR(SQR(l / cl2Toningcurve[l])); //apply C2=f(L) acts only on 'b'
    }

    if (algm == 1) {
        Color::interpolateRGBColor (realL, iplow, iphigh, algm, opacity, twoc, metchrom, chromat, luma, r, g, b, xl, yl, zl, x2, y2, z2, wp, wip, ro, go, bo);
    } else {
        Color::interpolateRGBColor (realL, iplow, iphigh, algm, opacity2, twoc, metchrom, chromat, luma, r, g, b, xl, yl, zl, x2, y2, z2, wp, wip, ro, go, bo);
    }
}


void ImProcFunctions::luminanceCurve (LabImage* lold, LabImage* lnew, LUTf & curve)
{

    int W = lold->W;
    int H = lold->H;

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            float Lin = lold->L[i][j];
            //if (Lin>0 && Lin<65535)
            lnew->L[i][j] = curve[Lin];
        }
}





//#include "cubic.cc"

//void ImProcFunctions::colorCurve (LabImage* lold, LabImage* lnew)
//{

    /*    LUT<double> cmultiplier(181021);

        double boost_a = ((float)params->colorBoost.amount + 100.0) / 100.0;
        double boost_b = ((float)params->colorBoost.amount + 100.0) / 100.0;

        double c, amul = 1.0, bmul = 1.0;
        if (boost_a > boost_b) {
            c = boost_a;
            if (boost_a > 0)
                bmul = boost_b / boost_a;
        }
        else {
            c = boost_b;
            if (boost_b > 0)
                amul = boost_a / boost_b;
        }

        if (params->colorBoost.enable_saturationlimiter && c>1.0) {
            // re-generate color multiplier lookup table
            double d = params->colorBoost.saturationlimit / 3.0;
            double alpha = 0.5;
            double threshold1 = alpha * d;
            double threshold2 = c*d*(alpha+1.0) - d;
            for (int i=0; i<=181020; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
                double chrominance = (double)i/4.0;
                if (chrominance < threshold1)
                    cmultiplier[i] = c;
                else if (chrominance < d)
                    cmultiplier[i] = (c / (2.0*d*(alpha-1.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
                else if (chrominance < threshold2)
                    cmultiplier[i] = (1.0 / (2.0*d*(c*(alpha+1.0)-2.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
                else
                    cmultiplier[i] = 1.0;
            }
        }

        float eps = 0.001;
        double shift_a = params->colorShift.a + eps, shift_b = params->colorShift.b + eps;

        float** oa = lold->a;
        float** ob = lold->b;

        #pragma omp parallel for if (multiThread)
        for (int i=0; i<lold->H; i++)
            for (int j=0; j<lold->W; j++) {

                double wanted_c = c;
                if (params->colorBoost.enable_saturationlimiter && c>1) {
                    float chroma = (float)(4.0 * sqrt((oa[i][j]+shift_a)*(oa[i][j]+shift_a) + (ob[i][j]+shift_b)*(ob[i][j]+shift_b)));
                    wanted_c = cmultiplier [chroma];
                }

                double real_c = wanted_c;
                if (wanted_c >= 1.0 && params->colorBoost.avoidclip) {
                    double cclip = 100000.0;
                    double cr = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, 3.079935, -1.5371515, -0.54278342);
                    double cg = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, -0.92123418, 1.87599, 0.04524418);
                    double cb = tightestroot ((double)lnew->L[i][j]/655.35, (double)(oa[i][j]+shift_a)*amul, (double)(ob[i][j]+shift_b)*bmul, 0.052889682, -0.20404134, 1.15115166);
                    if (cr>1.0 && cr<cclip) cclip = cr;
                    if (cg>1.0 && cg<cclip) cclip = cg;
                    if (cb>1.0 && cb<cclip) cclip = cb;
                    if (cclip<100000.0) {
                        real_c = -cclip + 2.0*cclip / (1.0+exp(-2.0*wanted_c/cclip));
                        if (real_c<1.0)
                            real_c = 1.0;
                    }
                }

                float nna = ((oa[i][j]+shift_a) * real_c * amul);
                float nnb = ((ob[i][j]+shift_b) * real_c * bmul);
                lnew->a[i][j] = LIM(nna,-32000.0f,32000.0f);
                lnew->b[i][j] = LIM(nnb,-32000.0f,32000.0f);
            }
    */
    //delete [] cmultiplier;
//}

void ImProcFunctions::impulsedenoise(Imagefloat *rgb)
{

    if (params->impulseDenoise.enabled && rgb->getWidth() >= 8 && rgb->getHeight() >= 8)

    {
        rgb->setMode(Imagefloat::Mode::LAB, multiThread);
        impulse_nr(rgb, (float)params->impulseDenoise.thresh / 20.0 );
    }
}

void ImProcFunctions::defringe(Imagefloat *rgb)
{
    if (params->defringe.enabled && rgb->getWidth() >= 8 && rgb->getHeight() >= 8)

    {
        rgb->setMode(Imagefloat::Mode::LAB, multiThread);
        PF_correct_RT(rgb, params->defringe.radius, params->defringe.threshold);
    }
}

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

        delete thumbGray;
        delete rawGray;
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

}
