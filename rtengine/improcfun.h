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
#include "shmap.h"
#include "coord2d.h"
#include "color.h"
#include "labimage.h"
#include "cieimage.h"
#include "LUT.h"
#include "lcp.h"
#include "dcp.h"
#include "curves.h"
#include "cplx_wavelet_dec.h"
#include "pipettebuffer.h"
#include "gamutwarning.h"

namespace rtengine {

using namespace procparams;

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
    ImProcFunctions(const ProcParams* iparams, bool imultiThread = true)
        : monitorTransform(nullptr), params(iparams), scale(1), multiThread(imultiThread), dcpProf(nullptr), dcpApplyState(nullptr), pipetteBuffer(nullptr), lumimul{} {}
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

    void rgbProc(Imagefloat *working);

    void labAdjustments(Imagefloat *rgb, LUTu *histCCurve=nullptr, LUTu *histLCurve=nullptr);
    
    bool sharpening(Imagefloat *rgb, const SharpeningParams &sharpenParam, bool showMask=false);
    
    void transform(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const FramesMetaData *metadata, int rawRotationDeg, bool fullImage);    
    
    void resize(Imagefloat* src, Imagefloat* dst, float dScale);
    void Lanczos(const LabImage* src, LabImage* dst, float scale);
    void Lanczos(Imagefloat *src, Imagefloat *dst, float scale);

    void impulsedenoise(Imagefloat *rgb);   //Emil's impulse denoise

//    void dirpyrequalizer(LabImage* lab, int scale);  //Emil's wavelet
//    void EPDToneMap(LabImage *lab, unsigned int Iterates = 0, int skip = 1);
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
    //void shadowsHighlights(LabImage *lab);
    void shadowsHighlights(Imagefloat *rgb);
    void toneEqualizer(Imagefloat *rgb);
    void softLight(Imagefloat *rgb);
    bool colorCorrection(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    void logEncoding(Imagefloat *rgb, LUTu *histToneCurve=nullptr);
    bool contrastByDetailLevels(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    void filmGrain(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
    bool guidedSmoothing(Imagefloat *rgb, int offset_x=0, int offset_y=0, int full_width=-1, int full_height=-1);
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
    
    void moyeqt(Imagefloat* working, float &moyS, float &eqty);
    
    void Tile_calc(int tilesize, int overlap, int kall, int imwidth, int imheight, int &numtiles_W, int &numtiles_H, int &tilewidth, int &tileheight, int &tileWskip, int &tileHskip);
    
    void getAutoLog(ImageSource *imgsrc, procparams::LogEncodingParams &params);
    
    static void getAutoExp(const LUTu & histogram, int histcompr, double clip, double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh);
    
    static double getAutoDistor(const Glib::ustring& fname, int thumb_size);
    bool transCoord(int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr);
    bool transCoord(int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr);
    double getTransformAutoFill(int oW, int oH, const LensCorrection *pLCPMap = nullptr);

    void calcautodn_info(float &chaut, float &delta, int Nb, int levaut, float maxmax, float lumema, float chromina, int mode, int lissage, float redyel, float skinc, float nsknc);
    void Median_Denoise(float **src, float **dst, float upperBound, int width, int height, Median medianType, int iterations, int numThreads, float **buffer = nullptr);
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
    
private:
    void RGB_denoise(int kall, Imagefloat * src, Imagefloat * dst, Imagefloat * calclum, float * ch_M, float *max_r, float *max_b, bool isRAW, const procparams::DenoiseParams & dnparams, const double expcomp, const NoiseCurve & noiseLCurve, const NoiseCurve & noiseCCurve, float &nresi, float &highresi);
    void RGB_denoise_infoGamCurve(const procparams::DenoiseParams & dnparams, const bool isRAW, LUTf &gamcurve, float &gam, float &gamthresh, float &gamslope);    
    void RGB_denoise_info(Imagefloat * src, Imagefloat * provicalc, bool isRAW, LUTf &gamcurve, float gam, float gamthresh, float gamslope, const procparams::DenoiseParams & dnparams, const double expcomp, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float & maxblueaut, float &minredaut, float & minblueaut, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc, bool multiThread = false);

    void denoiseGuidedSmoothing(Imagefloat *rgb);

    void chromiLuminanceCurve(LabImage* lold, LabImage* lnew, LUTf &acurve, LUTf &bcurve, LUTf & satcurve, LUTf & satclcurve, LUTf &clcurve, LUTf &curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu *histCCurve, LUTu *histLurve);    
    
    void labtoning(float r, float g, float b, float &ro, float &go, float &bo, int algm, int metchrom, int twoc, float satLimit, float satLimitOpacity, const ColorGradientCurve & ctColorCurve, const OpacityCurve & ctOpacityCurve, LUTf & clToningcurve, LUTf & cl2Toningcurve, float iplow, float iphigh, double wp[3][3], double wip[3][3]);
    void toning2col(float r, float g, float b, float &ro, float &go, float &bo, float iplow, float iphigh, float rl, float gl, float bl, float rh, float gh, float bh, float SatLow, float SatHigh, float balanS, float balanH, float reducac, int mode, int preser, float strProtect);
    void toningsmh(float r, float g, float b, float &ro, float &go, float &bo, float RedLow, float GreenLow, float BlueLow, float RedMed, float GreenMed, float BlueMed, float RedHigh, float GreenHigh, float BlueHigh, float reducac, int mode, float strProtect);
    void toningsmh2(float r, float g, float b, float &ro, float &go, float &bo, float low[3], float satLow, float med[3], float satMed, float high[3], float satHigh, float reducac, int mode, int preser);
    void secondeg_begin(float reducac, float vend, float &aam, float &bbm);
    void secondeg_end(float reducac, float vinf, float &aa, float &bb, float &cc);

    void retreavergb(float &r, float &g, float &b);

    void luminanceCurve(LabImage* lold, LabImage* lnew, LUTf &curve);

    void impulse_nr(Imagefloat *lab, double thresh);

    void dirpyrdenoise(LabImage* src);    //Emil's pyramid denoise

//    void EPDToneMapResid(float * WavCoeffs_L0, unsigned int Iterates,  int skip, struct cont_params& cp, int W_L, int H_L, float max0, float min0);
    void CompressDR(float *Source, int W_L, int H_L, float Compression, float DetailBoost);
    void ContrastResid(float * WavCoeffs_L0, struct cont_params &cp, int W_L, int H_L, float max0, float min0);

    // pyramid denoise
    // procparams::DenoiseParams dnparams;
    // void dirpyr(LabImage* data_fine, LabImage* data_coarse, int level, LUTf &rangefn_L, LUTf &rangefn_ab,
    //             int pitch, int scale, const int luma, int chroma);
    // void idirpyr(LabImage* data_coarse, LabImage* data_fine, int level, LUTf &rangefn_L, LUTf & nrwt_l, LUTf & nrwt_ab,
    //              int pitch, int scale, const int luma, const int chroma/*, LUTf & Lcurve, LUTf & abcurve*/);

    void Median_Denoise(float **src, float **dst, int width, int height, Median medianType, int iterations, int numThreads, float **buffer = nullptr);
    void RGBtile_denoise(float * fLblox, int hblproc, float noisevar_Ldetail, float * nbrwt, float * blurbuffer);     //for DCT
    void RGBoutput_tile_row(float *bloxrow_L, float ** Ldetail, float ** tilemask_out, int height, int width, int top);
    bool WaveletDenoiseAllL(wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3], float * vari, int edge);
    bool WaveletDenoiseAllAB(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_ab, float *noisevarchrom, float madL[8][3], float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb);
    void WaveletDenoiseAll_info(int levwav, wavelet_decomposition &WaveletCoeffs_a,
                                wavelet_decomposition &WaveletCoeffs_b, float **noisevarlum, float **noisevarchrom, float **noisevarhue, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float & minblueaut, int schoice, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
                                float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb);

    bool WaveletDenoiseAll_BiShrinkL(wavelet_decomposition &WaveletCoeffs_L, float *noisevarlum, float madL[8][3]);
    bool WaveletDenoiseAll_BiShrinkAB(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_ab, float *noisevarchrom, float madL[8][3], float noisevar_ab,
                                      const bool useNoiseCCurve,  bool autoch, bool denoiseMethodRgb);
    void ShrinkAllL(wavelet_decomposition &WaveletCoeffs_L, float **buffer, int level, int dir, float *noisevarlum, float * madL, float * vari, int edge);
    void ShrinkAllAB(wavelet_decomposition &WaveletCoeffs_L, wavelet_decomposition &WaveletCoeffs_ab, float **buffer, int level, int dir,
                     float *noisevarchrom, float noisevar_ab, const bool useNoiseCCurve, bool autoch, bool denoiseMethodRgb, float * madL, float * madaab = nullptr, bool madCalculated = false);
    void ShrinkAll_info(float ** WavCoeffs_a, float ** WavCoeffs_b,
                        int W_ab, int H_ab, float **noisevarlum, float **noisevarchrom, float **noisevarhue, float &chaut, int &Nb, float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, int schoice, int lvl, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc,
                        float &maxchred, float &maxchblue, float &minchred, float &minchblue, int &nb, float &chau, float &chred, float &chblue, bool denoiseMethodRgb);
    void Noise_residualAB(wavelet_decomposition &WaveletCoeffs_ab, float &chresid, float &chmaxresid, bool denoiseMethodRgb);
    float MadMax(float * DataList, int &max, int datalen);
    float Mad(float * DataList, const int datalen);
    float MadRgb(float * DataList, const int datalen);

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
