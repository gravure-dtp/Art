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
#include "sleef.h"
#include "curves.h"
#include "linalgebra.h"

#include <fstream>

namespace rtengine {

namespace {

template <class Curve>
inline void apply(const Curve &c, Imagefloat *rgb, int W, int H, bool multithread)
{
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            c.Apply(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x));
        }
    }
}


void apply_tc(Imagefloat *rgb, const ToneCurve &tc, ToneCurveParams::TcMode curveMode, const Glib::ustring &working_profile, int perceptual_strength, float whitept, bool multithread)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    
    if (curveMode == ToneCurveParams::TcMode::PERCEPTUAL) {
        const PerceptualToneCurve &c = static_cast<const PerceptualToneCurve&>(tc);
        PerceptualToneCurveState state;
        c.initApplyState(state, working_profile);
        state.strength = LIM01(float(perceptual_strength) / 100.f);

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            c.BatchApply(0, W, rgb->r.ptrs[y], rgb->g.ptrs[y], rgb->b.ptrs[y], state);
        }
    } else if (curveMode == ToneCurveParams::TcMode::STD) {
        const StandardToneCurve &c = static_cast<const StandardToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::WEIGHTEDSTD) {
        const WeightedStdToneCurve &c = static_cast<const WeightedStdToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::FILMLIKE) {
        const AdobeToneCurve &c = static_cast<const AdobeToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::SATANDVALBLENDING) {
        const SatAndValueBlendingToneCurve &c = static_cast<const SatAndValueBlendingToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::LUMINANCE) {
        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);
        const LuminanceToneCurve &c = static_cast<const LuminanceToneCurve &>(tc);
//        apply(c, rgb, W, H, multithread);
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                c.Apply(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), ws);
            }
        }
    } else if (curveMode == ToneCurveParams::TcMode::NEUTRAL) {
        const NeutralToneCurve &c = static_cast<const NeutralToneCurve &>(tc);
        NeutralToneCurve::ApplyState state(working_profile, whitept);

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            c.BatchApply(0, W, rgb->r.ptrs[y], rgb->g.ptrs[y], rgb->b.ptrs[y], state);
        }
    }
}


class ContrastCurve: public Curve {
public:
    ContrastCurve(double a, double b, double w): a_(a), b_(b), w_(w) {}
    void getVal(const std::vector<double>& t, std::vector<double>& res) const {}
    bool isIdentity () const { return false; }
    
    double getVal(double x) const
    {
        double res = lin2log(std::pow(x/w_, a_), b_)*w_;
        return res;
    }

private:
    double a_;
    double b_;
    double w_;
};


void filmlike_clip(Imagefloat *rgb, float whitept, bool multithread)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    const float Lmax = 65535.f * whitept;

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            float &r = rgb->r(i, j);
            float &g = rgb->g(i, j);
            float &b = rgb->b(i, j);
            Color::filmlike_clip(&r, &g, &b, Lmax);
        }
    }
}


void legacy_contrast(Imagefloat *rgb, const ImProcData &im, int contrast, const Glib::ustring &working_profile, float whitept)
{
    if (contrast) {
        ToneCurve tc;
        auto &curve = tc.lutToneCurve;
        curve(65536);

        tc.Set(DiagonalCurve({DCT_Empty}));
        
        LUTu hist16(65536);
        ImProcFunctions ipf(im.params, im.multiThread);
        ipf.firstAnalysis(rgb, *im.params, hist16);

        CurveFactory::contrastCurve(contrast, hist16, curve, max(im.scale, 1.0));
        apply_tc(rgb, tc, ToneCurveParams::TcMode::STD, working_profile, 100, whitept, im.multiThread);
    }
}


std::unique_ptr<Curve> get_contrast_curve(Imagefloat *rgb, const ImProcData &im, int contrast, float whitept)
{
    std::unique_ptr<Curve> ccurve;
    
    if (contrast) {
        const double pivot = (im.params->logenc.enabled ? im.params->logenc.targetGray / 100.0 : 0.18) / whitept;
        const double c = std::pow(std::abs(contrast) / 100.0, 1.5) * 16.0;
        const double b = contrast > 0 ? (1 + c) : 1.0 / (1 + c);
        const double a = std::log((std::exp(std::log(b) * pivot) - 1) / (b - 1)) / std::log(pivot);

        ccurve.reset(new ContrastCurve(a, b, whitept));
    }
    return ccurve;
}


void satcurve_lut(const FlatCurve &curve, LUTf &sat, float whitept)
{
    sat(65536, LUT_CLIP_BELOW);
    sat[0] = curve.getVal(0) * 2.f;
    for (int i = 1; i < 65536; ++i) {
        float x = Color::gamma2curve[i] / 65535.f;
        float v = curve.getVal(x);
        sat[i] = v * 2.f;
    }
}


void apply_satcurve(Imagefloat *rgb, const FlatCurve &curve, const Glib::ustring &working_profile, float whitept, bool multithread)
{
    LUTf sat;
    satcurve_lut(curve, sat, whitept);

    // if (whitept > 1.f) {
        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);
        TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(working_profile);

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < rgb->getHeight(); ++y) {
            float X, Y, Z;
            float Jz, az, bz;
            for (int x = 0; x < rgb->getWidth(); ++x) {
                Color::rgbxyz(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), X, Y, Z, ws);
                Color::xyz2jzazbz(X, Y, Z, Jz, az, bz);
                float s = sat[Y];
                az *= s;
                bz *= s;
                Color::jzazbz2rgb(Jz, az, bz, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), iws);
            }
        }
//     } else {
//         rgb->setMode(Imagefloat::Mode::LAB, multithread);
// #ifdef _OPENMP
// #       pragma omp parallel for if (multithread)
// #endif
//         for (int y = 0; y < rgb->getHeight(); ++y) {
//             for (int x = 0; x < rgb->getWidth(); ++x) {
//                 float X, Y, Z;
//                 Color::L2XYZ(rgb->g(y, x), X, Y, Z);
//                 float s = sat[Y];
//                 rgb->r(y, x) *= s;
//                 rgb->b(y, x) *= s;
//             }
//         }
//         rgb->setMode(Imagefloat::Mode::RGB, multithread);
//     }
}


void fill_satcurve_pipette(Imagefloat *rgb, PlanarWhateverData<float>* editWhatever, const Glib::ustring &working_profile, float whitept, bool multithread)
{
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < rgb->getHeight(); ++y) {
        for (int x = 0; x < rgb->getWidth(); ++x) {
            float r = rgb->r(y, x), g = rgb->g(y, x), b = rgb->b(y, x);
            float Y = Color::rgbLuminance(r, g, b, ws);
            float s = Color::gamma2curve[Y] / 65535.f;
            editWhatever->v(y, x) = LIM01(s);
        }
    }
}


void update_tone_curve_histogram(Imagefloat *img, LUTu &hist, const Glib::ustring &profile, bool multithread)
{
    hist.clear();
    const int compression = log2(65536 / hist.getSize());

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(profile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < img->getHeight(); ++y) {
        for (int x = 0; x < img->getWidth(); ++x) {
            float r = CLIP(img->r(y, x));
            float g = CLIP(img->g(y, x));
            float b = CLIP(img->b(y, x));

            int y = CLIP<int>(Color::gamma2curve[Color::rgbLuminance(r, g, b, ws)]);//max(r, g, b)]);
            hist[y >> compression]++;
        }
    }

    // we make this log encoded
    int n = hist.getSize();
    float f = float(n);
    for (int i = 0; i < n; ++i) {
        hist[i] = xlin2log(float(hist[i]) / f, 2.f) * f;
    }
}

void fill_pipette(Imagefloat *img, Imagefloat *pipette, bool multithread)
{
    const int W = img->getWidth();
    const int H = img->getHeight();
    
#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            pipette->r(y, x) = Color::gamma2curve[CLIP(img->r(y, x))] / 65535.f;
            pipette->g(y, x) = Color::gamma2curve[CLIP(img->g(y, x))] / 65535.f;
            pipette->b(y, x) = Color::gamma2curve[CLIP(img->b(y, x))] / 65535.f;
        }
    }
}


class DoubleCurve: public Curve {
public:
    DoubleCurve(const Curve &c1, const Curve &c2):
        c1_(c1), c2_(c2) {}

    double getVal(double t) const override
    {
        return c2_.getVal(c1_.getVal(t));
    }
    
    void getVal(const std::vector<double>& t, std::vector<double>& res) const override
    {
        c1_.getVal(t, res);
        c2_.getVal(res, res);
    }

    bool isIdentity() const override
    {
        return c1_.isIdentity() && c2_.isIdentity();
    }
    
private:
    const Curve &c1_;
    const Curve &c2_;
};


void gamut_compression(Imagefloat *img, const Glib::ustring &outprofile, float whitept, bool multithread)
{
    Mat33<float> om;
    if (!ICCStore::getInstance()->getProfileMatrix(outprofile, om)) {
        return;
    }
    auto iom = inverse(om);
    auto ws = ICCStore::getInstance()->workingSpaceMatrix(img->colorSpace());
    auto to_out = dot_product(ws, iom);
    auto iws = ICCStore::getInstance()->workingSpaceInverseMatrix(img->colorSpace());
    auto to_work = dot_product(om, iws);

    const int W = img->getWidth();
    const int H = img->getHeight();
    const float Lmax = whitept * 65535.f;

    // from https://github.com/jedypod/gamut-compress

    // Distance limit: How far beyond the gamut boundary to compress
    //const Vec3<float> dl(1.147f, 1.264f, 1.312f); // original ACES values
    const Vec3<float> dl(1.1f, 1.2f, 1.5f); // hand-tuned
    // Amount of outer gamut to affect
    //const Vec3<float> th(0.815f, 0.803f, 0.88f); // original ACES values
    const Vec3<float> th(0.85f, 0.75f, 0.95f); // hand-tuned

    // Power or Parabolic compression functions: https://www.desmos.com/calculator/nvhp63hmtj
#if 0 // power compression
    constexpr float p = 1.2f;
    const auto scale =
        [](float l, float t, float p) -> float
        {
            return (l - t) / std::pow(std::pow((1.f - t)/(l - t), -p) - 1.f, 1.f/p);
        };
    const Vec3<float> s(scale(dl[0], th[0], p),
                        scale(dl[1], th[1], p),
                        scale(dl[2], th[2], p));


    const auto compr =
        [&](float x, int i) -> float
        {
            float t = (x - th[i])/s[i];
            return th[i] + s[i] * std::pow(t / (1.f + std::pow(t, p)), 1.f/p);
        };
#else // parabolic compression
    const auto scale =
        [](float l, float t) -> float
        {
            return (1.f - t) / std::sqrt(l-1.f);
        };
    Vec3<float> s(scale(dl[0], th[0]), scale(dl[1], th[1]), scale(dl[2], th[2]));

    const auto compr =
        [&](float x, int i) -> float
        {
            return s[i] * std::sqrt(x - th[i] + SQR(s[i])/4.0f) - s[i] * std::sqrt(SQR(s[i]) / 4.0f) + th[i];            
        };
#endif // power/parabolic compression
  
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Vec3<float> rgb(img->r(y, x), img->g(y, x), img->b(y, x));
            float iY = (rgb[0] + rgb[1] + rgb[2]) / 3.f;

            rgb = dot_product(to_out, rgb);

            // Achromatic axis
            float ac = max(rgb[0], rgb[1], rgb[2]);

            // Inverse RGB Ratios: distance from achromatic axis
            Vec3<float> d;
            float aac = std::abs(ac);
            if (ac != 0.f) {
                d[0] = (ac - rgb[0])/aac;
                d[1] = (ac - rgb[1])/aac;
                d[2] = (ac - rgb[2])/aac;
            }

            Vec3<float> cd; // Compressed distance
            for (int i = 0; i < 3; ++i) {
                cd[i] = d[i] < th[i] ? d[i] : compr(d[i], i);
            }

            // Inverse RGB Ratios to RGB
            rgb[0] = ac-cd[0]*aac;
            rgb[1] = ac-cd[1]*aac;
            rgb[2] = ac-cd[2]*aac;
          
            rgb = dot_product(to_work, rgb);

            float oY = (rgb[0] + rgb[1] + rgb[2]) / 3.f;
            if (oY > 0.f) {
                float f = iY / oY;
                rgb[0] *= f;
                rgb[1] *= f;
                rgb[2] *= f;
                Color::filmlike_clip(&rgb[0], &rgb[1], &rgb[2], Lmax);
            }
            
            img->r(y, x) = rgb[0];
            img->g(y, x) = rgb[1];
            img->b(y, x) = rgb[2];
        }
    }
}

} // namespace


void ImProcFunctions::toneCurve(Imagefloat *img)
{
    if (histToneCurve && *histToneCurve) {
        img->setMode(Imagefloat::Mode::RGB, multiThread);
        update_tone_curve_histogram(img, *histToneCurve, params->icm.workingProfile, multiThread);
    }

    Imagefloat *editImgFloat = nullptr;
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID editID = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((editID == EUID_ToneCurve1 || editID == EUID_ToneCurve2) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_IMAGEFLOAT) {
        editImgFloat = pipetteBuffer->getImgFloatBuffer();
    } else if (editID == EUID_ToneCurveSaturation && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }

    if (params->toneCurve.enabled) {
        img->setMode(Imagefloat::Mode::RGB, multiThread);


        const float whitept = params->toneCurve.hasWhitePoint() ? params->toneCurve.whitePoint : 1.f;

        const bool single_curve = params->toneCurve.curveMode == params->toneCurve.curveMode2;
        
        ImProcData im(params, scale, multiThread);
        if (single_curve && params->toneCurve.curveMode == ToneCurveParams::TcMode::NEUTRAL) {
            gamut_compression(img, params->icm.outputProfile, whitept, multiThread);
        } else {
            filmlike_clip(img, whitept, im.multiThread);
        }

        std::unique_ptr<Curve> ccurve;
        if (params->toneCurve.contrastLegacyMode) {
            legacy_contrast(img, im, params->toneCurve.contrast, params->icm.workingProfile, whitept);
        } else {
            ccurve = get_contrast_curve(img, im, params->toneCurve.contrast, whitept);
        }

        const auto expand =
            [whitept](double x) -> double
            {
                if (whitept <= 1.1f) {
                    return x;
                }

                double f = (pow_F(whitept, x) - 1) / (whitept - 1);
                double g = rtengine::intp(SQR(x)*x, f * whitept, x);
                return g;
            };
        
        const auto adjust =
            [whitept,&expand](std::vector<double> c) -> std::vector<double>
            {
                std::map<double, double> m;
                DiagonalCurveType tp = DiagonalCurveType(c[0]);
                bool add_c = (tp == DCT_CatmullRom || tp == DCT_Spline);
                DiagonalCurve curve(c);
                for (int i = 0; i < 25; ++i) {
                    double x = double(i)/100.0;
                    double v = Color::gammatab_srgb[x * 65535.0] / 65535.0;
                    double y = curve.getVal(v);
                    y = Color::igammatab_srgb[y * 65535.0] / 65535.0;
                    m[expand(x)] = expand(y);
                }
                for (int i = 25, j = 2; i < 100; ) {
                    double x = double(i)/100.0;
                    double v = Color::gammatab_srgb[x * 65535.0] / 65535.0;
                    double y = curve.getVal(v);
                    y = Color::igammatab_srgb[y * 65535.0] / 65535.0;
                    m[expand(x)] = expand(y);
                    i += j;
                    j *= 2;
                }
                if (add_c) {
                    for (size_t i = 0; i < (c.size()-2)/2; ++i) {
                        double x = c[2*i+1];
                        double v = Color::gammatab_srgb[x * 65535.0] / 65535.0;
                        double y = curve.getVal(v);
                        y = Color::igammatab_srgb[y * 65535.0] / 65535.0;
                        m[expand(x)] = expand(y);
                    }
                }
                m[expand(1.0)] = expand(curve.getVal(1.0));
                c = { DCT_CatmullRom };
                for (auto &p : m) {
                    c.push_back(p.first);
                    c.push_back(p.second);
                }
#if 0
                auto name = "/tmp/curve" + std::to_string(int(whitept)) + ".txt";
                std::ofstream out(name.c_str());
                for (size_t i = 1; i < c.size(); i += 2) {
                    out << c[i] << " " << c[i+1] << "\n";
                }
#endif // if 0
                return c;
            };

        ToneCurve tc;
        DiagonalCurve tcurve2(adjust(params->toneCurve.curve2), CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        DiagonalCurve tcurve1(adjust(params->toneCurve.curve), CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        DoubleCurve dcurve(tcurve1, tcurve2);
        std::unique_ptr<Curve> dccurve;
        Curve *tcurve = &dcurve;
        if (ccurve) {
            dccurve.reset(new DoubleCurve(*ccurve, dcurve));
            tcurve = dccurve.get();
        }

        if (single_curve) {
            if (editImgFloat && (editID == EUID_ToneCurve1 || editID == EUID_ToneCurve2)) {
                fill_pipette(img, editImgFloat, multiThread);
            }
            if (!tcurve->isIdentity()) {
                tc.Set(*tcurve, 65535.f * whitept);
                apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, params->toneCurve.perceptualStrength, whitept, multiThread);
            }
        } else {
            if (ccurve) {
                tc.Set(*ccurve, 65535.f * whitept);
                apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, 100, whitept, multiThread);
            }
            
            if (editImgFloat && editID == EUID_ToneCurve1) {
                fill_pipette(img, editImgFloat, multiThread);
            }
        
            if (!tcurve1.isIdentity()) {
                tc.Set(tcurve1, 65535.f * whitept);
                apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, params->toneCurve.perceptualStrength, whitept, multiThread);
            }

            if (editImgFloat && editID == EUID_ToneCurve2) {
                fill_pipette(img, editImgFloat, multiThread);
            }

            if (!tcurve2.isIdentity()) {
                tc.Set(tcurve2, 65535.f * whitept);
                apply_tc(img, tc, params->toneCurve.curveMode2, params->icm.workingProfile, params->toneCurve.perceptualStrength, whitept, multiThread);
            }
        }

        if (editWhatever) {
            fill_satcurve_pipette(img, editWhatever, params->icm.workingProfile, whitept, multiThread);
        }

        const FlatCurve satcurve(params->toneCurve.saturation, false, CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        if (!satcurve.isIdentity()) {
            apply_satcurve(img, satcurve, params->icm.workingProfile, whitept, multiThread);
        }
    } else if (editImgFloat) {
        const int W = img->getWidth();
        const int H = img->getHeight();

#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            std::fill(editImgFloat->r(y), editImgFloat->r(y)+W, 0.f);
            std::fill(editImgFloat->g(y), editImgFloat->g(y)+W, 0.f);
            std::fill(editImgFloat->b(y), editImgFloat->b(y)+W, 0.f);
        }
    } else if (editWhatever) {
        editWhatever->fill(0.f);
    }
}

} // namespace rtengine
