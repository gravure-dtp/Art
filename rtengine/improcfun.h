/* -*- C++ -*-
 *  
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
#ifndef _IMPROCFUN_H_
#define _IMPROCFUN_H_

#include "imagefloat.h"
#include "image16.h"
#include "image8.h"
#include "procparams.h"
#include "coord2d.h"
#include "color.h"
#include "labimage.h"
#include "LUT.h"
#include "lcp.h"
#include "dcp.h"
#include "curves.h"
#include "cplx_wavelet_dec.h"
#include "pipettebuffer.h"
#include "gamutwarning.h"

namespace rtengine {

using namespace procparams;

struct ImProcData {
    const ProcParams *params;
    double scale;
    bool multiThread;

    explicit ImProcData(const ProcParams *p=nullptr, double s=1.0, bool m=true):
        params(p), scale(s), multiThread(m) {}
};


class ImProcFunctions {
public:
    enum class Median {
        TYPE_3X3_SOFT,
        TYPE_3X3_STRONG,
        TYPE_5X5_SOFT,
        TYPE_5X5_STRONG,
        TYPE_7X7,
        TYPE_9X9
    };

    //----------------------------------------------------------------------
    // constructor/destructor and initialization/state manipulation
    //----------------------------------------------------------------------
    ImProcFunctions(const ProcParams* iparams, bool imultiThread=true);
    ~ImProcFunctions();
    
    void setScale(double iscale);

    void updateColorProfiles(const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, bool gamutCheck);

    void setDCPProfile(DCPProfile *dcp, const DCPProfile::ApplyState &as)
    {
        dcpProf = dcp;
        dcpApplyState = &as;
    }

    void setPipetteBuffer(PipetteBuffer *pb)
    {
        pipetteBuffer = pb;
    }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // image processing operations
    //----------------------------------------------------------------------
    void firstAnalysis(const Imagefloat* const working, const ProcParams &params, LUTu & vhist16);

    void labAdjustments(Imagefloat *rgb, LUTu *histCCurve=nullptr, LUTu *histLCurve=nullptr);
    bool sharpening(Imagefloat *rgb, const SharpeningParams &sharpenParam, bool showMask=false);
    void transform(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const FramesMetaData *metadata, int rawRotationDeg, bool fullImage);    
    void resize(Imagefloat* src, Imagefloat* dst, float dScale);
    void Lanczos(const LabImage* src, LabImage* dst, float scale);
    void Lanczos(Imagefloat *src, Imagefloat *dst, float scale);
    void impulsedenoise(Imagefloat *rgb);   //Emil's impulse denoise
    bool textureBoost(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);

    struct DenoiseInfoStore {
        DenoiseInfoStore () : chM (0), max_r{}, max_b{}, ch_M{}, valid (false)  {}
        float chM;
        float max_r[9];
        float max_b[9];
        float ch_M[9];
        bool valid;
    };
    void denoiseComputeParams(ImageSource *imgsrc, const ColorTemp &currWB, DenoiseInfoStore &store, procparams::DenoiseParams &dnparams);
    void denoise(ImageSource *imgsrc, const ColorTemp &currWB, Imagefloat *img, const DenoiseInfoStore &store, const procparams::DenoiseParams &dnparams);
    
    void defringe(Imagefloat *rgb);
    void dehaze(Imagefloat *rgb);
    void dynamicRangeCompression(Imagefloat *rgb);
    void localContrast(Imagefloat *rgb);
    void shadowsHighlights(Imagefloat *rgb);
    void toneEqualizer(Imagefloat *rgb);
    void softLight(Imagefloat *rgb);
    bool colorCorrection(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    void logEncoding(Imagefloat *rgb);
    bool contrastByDetailLevels(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    void filmGrain(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    bool guidedSmoothing(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    void hslEqualizer(Imagefloat *rgb);
    void channelMixer(Imagefloat *rgb);
    void exposure(Imagefloat *rgb);
    void rgbCurves(Imagefloat *rgb);
    void blackAndWhite(Imagefloat *rgb);
    void toneCurve(Imagefloat *img, LUTu *histToneCurve=nullptr);
    void brightnessContrastSaturation(Imagefloat *img);
    void filmSimulation(Imagefloat *img);

    enum class Stage {
        STAGE_0,
        STAGE_1,
        STAGE_2,
        STAGE_3
    };
    enum class Pipeline {
        THUMBNAIL,
        NAVIGATOR,
        PREVIEW,
        OUTPUT
    };
    bool process(Pipeline pipeline, Stage stage, Imagefloat *img);

    void setViewport(int ox, int oy, int fw, int fh);
    void setOutputHistograms(LUTu *histToneCurve, LUTu *histCCurve, LUTu *histLCurve);
    void setShowSharpeningMask(bool yes);
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // Lab/RGB conversion
    //----------------------------------------------------------------------
    void lab2monitorRgb(Imagefloat *img, Image8* image);
    
    Image8 *lab2rgb(Imagefloat *img, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, bool consider_histogram_settings = true);

    Imagefloat *lab2rgbOut(Imagefloat *img, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm);

    void rgb2lab(Imagefloat &src, LabImage &dst, const Glib::ustring &workingSpace);
    void rgb2lab(Imagefloat &src, LabImage &dst) { rgb2lab(src, dst, params->icm.workingProfile); }
    
    void lab2rgb(const LabImage &src, Imagefloat &dst, const Glib::ustring &workingSpace);    
    void lab2rgb(const LabImage &src, Imagefloat &dst) { lab2rgb(src, dst, params->icm.workingProfile); }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // misc helper functions for image processing ops
    //----------------------------------------------------------------------
    bool needsLuminanceOnly()
    {
        return !(needsCA() || needsDistortion() || needsRotation() || needsPerspective() || needsLCP() || needsLensfun()) && (needsVignetting() || needsPCVignetting() || needsGradient());
    }

    bool needsTransform();
    
    bool needsPCVignetting();
    
    float resizeScale(const ProcParams* params, int fw, int fh, int &imw, int &imh);
    
    void getAutoLog(ImageSource *imgsrc, procparams::LogEncodingParams &params);
    
    static void getAutoExp(const LUTu & histogram, int histcompr, double clip, double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh);
    
    static double getAutoDistor(const Glib::ustring& fname, int thumb_size);
    bool transCoord(int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr);
    bool transCoord(int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr);
    double getTransformAutoFill(int oW, int oH, const LensCorrection *pLCPMap = nullptr);
    //----------------------------------------------------------------------

private:
    cmsHTRANSFORM monitorTransform;
    std::unique_ptr<GamutWarning> gamutWarning;

    const ProcParams* params;
    double scale;
    bool multiThread;

    DCPProfile *dcpProf;
    const DCPProfile::ApplyState *dcpApplyState;

    PipetteBuffer *pipetteBuffer;

    double lumimul[3];

    int offset_x;
    int offset_y;
    int full_width;
    int full_height;

    LUTu *hist_tonecurve;
    LUTu *hist_ccurve;
    LUTu *hist_lcurve;

    bool show_sharpening_mask;
    
private:
    void impulse_nr(Imagefloat *lab, double thresh);

    void chromiLuminanceCurve(LabImage* lold, LabImage* lnew, LUTf &acurve, LUTf &bcurve, LUTf & satcurve, LUTf & satclcurve, LUTf &clcurve, LUTf &curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu *histCCurve, LUTu *histLurve);
    
    // pyramid wavelet
    void dirpyr_equalizer(float ** src, float ** dst, int srcwidth, int srcheight, float ** l_a, float ** l_b, const double * mult, const double dirpyrThreshold, const double skinprot, float b_l, float t_l, float t_r, float scale);    //Emil's directional pyramid wavelet
    void dirpyr_channel(float ** data_fine, float ** data_coarse, int width, int height, int level, int scale);
    void idirpyr_eq_channel(float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, float multi[6], const double dirpyrThreshold, float ** l_a_h, float ** l_b_c, const double skinprot, float b_l, float t_l, float t_r);

    void PF_correct_RT(Imagefloat *lab, double radius, int thresh);

    void calcVignettingParams(int oW, int oH, const VignettingParams& vignetting, double &w2, double &h2, double& maxRadius, double &v, double &b, double &mul);

    void transformLuminanceOnly(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH, int fW, int fH);
    void transformGeneral(bool highQuality, Imagefloat *original, Imagefloat *transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const LensCorrection *pLCPMap);
    void transformLCPCAOnly(Imagefloat *original, Imagefloat *transformed, int cx, int cy, const LensCorrection *pLCPMap);

    bool needsCA();
    bool needsDistortion();
    bool needsRotation();
    bool needsPerspective();
    bool needsGradient();
    bool needsVignetting();
    bool needsLCP();
    bool needsLensfun();
//   static cmsUInt8Number* Mempro = NULL;
};


} // namespace rtengine

#endif
