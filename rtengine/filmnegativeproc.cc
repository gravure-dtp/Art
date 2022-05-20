/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Romei <aldrop8@gmail.com>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <cmath>
#include <iostream>

#include "rawimage.h"
#include "rawimagesource.h"
#include "stdimagesource.h"
#include "imagefloat.h"

#include "coord.h"
#include "mytime.h"
#include "opthelper.h"
#include "pixelsmap.h"
#include "procparams.h"
#include "rt_algo.h"
#include "rtengine.h"
#include "sleef.h"
#include "settings.h"
//#define BENCHMARK
#include "StopWatch.h"

namespace rtengine {

extern const Settings *settings;

namespace {

bool channelsAvg(const rtengine::RawImage* ri, int width, int height, const float* cblacksom, rtengine::Coord spotPos, int spotSize, std::array<float, 3> &avgs)
{
    avgs = {}; // Channel averages

    if (settings->verbose > 1) {
        printf("Spot coord:  x=%d y=%d\n", spotPos.x, spotPos.y);
    }

    const int half_spot_size = spotSize / 2;

    const int& x1 = spotPos.x - half_spot_size;
    const int& x2 = spotPos.x + half_spot_size;
    const int& y1 = spotPos.y - half_spot_size;
    const int& y2 = spotPos.y + half_spot_size;

    if (x1 < 0 || x2 > width || y1 < 0 || y2 > height) {
        return false; // Spot goes outside bounds, bail out.
    }

    std::array<int, 3> pxCount = {}; // Per-channel sample counts

    if (ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS) {
        for (int c = x1; c < x2; ++c) {
            for (int r = y1; r < y2; ++r) {
                const int ch = ri->getSensorType() == ST_BAYER ? ri->FC(r, c) : ri->XTRANSFC(r, c);

                ++pxCount[ch];

                // Sample the original unprocessed values from RawImage, subtracting black levels.
                // Scaling is irrelevant, as we are only interested in the ratio between two spots.
                avgs[ch] += ri->data[r][c] - cblacksom[ch];
            }
        }
    } else if (ri->get_colors() == 3) {
        for (int c = x1; c < x2; ++c) {
            for (int r = y1; r < y2; ++r) {
                for (int ch = 0; ch < 3; ++ch) {
                    ++pxCount[ch];
                    avgs[ch] += ri->data[r][3*c + ch] - cblacksom[c];
                }
            }
        }
    } else if (ri->get_colors() == 1) {
        for (int c = x1; c < x2; ++c) {
            for (int r = y1; r < y2; ++r) {
                ++pxCount[0];
                avgs[0] += ri->data[r][c] - cblacksom[0];
            }
            pxCount[1] = pxCount[2] = pxCount[0];
            avgs[1] = avgs[2] = avgs[0];
        }        
    } else {
        assert(false);
        return false;
    }

    for (int ch = 0; ch < 3; ++ch) {
        avgs[ch] /= pxCount[ch];
    }

    return true;
}


bool channelsAvg(const rtengine::Imagefloat* img, int width, int height, rtengine::Coord spotPos, int spotSize, std::array<float, 3>& avgs)
{
    avgs = {}; // Channel averages

    if (settings->verbose > 1) {
        printf("Spot coord:  x=%d y=%d\n", spotPos.x, spotPos.y);
    }

    const int half_spot_size = spotSize / 2;

    const int& x1 = spotPos.x - half_spot_size;
    const int& x2 = spotPos.x + half_spot_size;
    const int& y1 = spotPos.y - half_spot_size;
    const int& y2 = spotPos.y + half_spot_size;

    if (x1 < 0 || x2 > width || y1 < 0 || y2 > height) {
        return false; // Spot goes outside bounds, bail out.
    }

    for (int c = x1; c < x2; ++c) {
        for (int r = y1; r < y2; ++r) {
            avgs[0] += img->r(r,c);
            avgs[1] += img->g(r,c);
            avgs[2] += img->b(r,c);
        }
    }

    for (int ch = 0; ch < 3; ++ch) {
        avgs[ch] /= (spotSize*spotSize);
    }

    return true;
}


void calcMedians(const rtengine::RawImage *ri, float **data, int x1, int y1, int x2, int y2, std::array<float, 3> &meds)
{
    MyTime t1, t2, t3;
    t1.set();

    // Channel vectors to calculate medians
    std::array<std::vector<float>, 3> cvs;

    // Sample one every 5 pixels, and push the value in the appropriate channel vector.
    // Choose an odd step, not a multiple of the CFA size, to get a chance to visit each channel.
    if (ri->getSensorType() == ST_BAYER) {
        for (int row = y1; row < y2; row += 5) {
            const int c0 = ri->FC(row, x1 + 0);
            const int c1 = ri->FC(row, x1 + 5);
            int col = x1;

            for (; col < x2 - 5; col += 10) {
                cvs[c0].push_back(data[row][col]);
                cvs[c1].push_back(data[row][col + 5]);
            }

            if (col < x2) {
                cvs[c0].push_back(data[row][col]);
            }
        }
    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
        for (int row = y1; row < y2; row += 5) {
            const std::array<unsigned int, 6> cs = {
                ri->XTRANSFC(row, x1 + 0),
                ri->XTRANSFC(row, x1 + 5),
                ri->XTRANSFC(row, x1 + 10),
                ri->XTRANSFC(row, x1 + 15),
                ri->XTRANSFC(row, x1 + 20),
                ri->XTRANSFC(row, x1 + 25)
            };
            int col = x1;

            for (; col < x2 - 25; col += 30) {
                for (int c = 0; c < 6; ++c) {
                    cvs[cs[c]].push_back(data[row][col + c * 5]);
                }
            }

            for (int c = 0; col < x2; col += 5, ++c) {
                cvs[cs[c]].push_back(data[row][col]);
            }
        }
    } else if (ri->get_colors() == 3) {
        for (int row = y1; row < y2; row += 5) {
            for (int col = x1; col < x2; col += 5) {
                for (int c = 0; c < 3; ++c) {
                    cvs[c].push_back(data[row][3*col+c]);
                }
            }
        }
    } else if (ri->get_colors() == 1) {
        for (int row = y1; row < y2; row += 5) {
            for (int col = x1; col < x2; col += 5) {
                for (int c = 0; c < 3; ++c) {
                    cvs[c].push_back(data[row][col]);
                }
            }
        }        
    }

    t2.set();

    if (settings->verbose > 1) {
        printf("Median vector fill loop time us: %d\n", t2.etime(t1));
    }

    t2.set();

    for (int c = 0; c < 3; ++c) {
        // Find median values for each channel
        if (!cvs[c].empty()) {
            rtengine::findMinMaxPercentile(cvs[c].data(), cvs[c].size(), 0.5f, meds[c], 0.5f, meds[c], true);
        }
    }

    t3.set();

    if (settings->verbose  > 1) {
        printf("Sample count: R=%zu, G=%zu, B=%zu\n", cvs[0].size(), cvs[1].size(), cvs[2].size());
        printf("Median calc time us: %d\n", t3.etime(t2));
    }

}

std::array<double, 3> calcWBMults(const rtengine::ColorTemp& wb, const rtengine::ImageMatrices& imatrices, const rtengine::RawImage *ri, const float ref_pre_mul[4])
{
    std::array<double, 3> wb_mul;
    double r, g, b;
    wb.getMultipliers(r, g, b);
    wb_mul[0] = imatrices.cam_rgb[0][0] * r + imatrices.cam_rgb[0][1] * g + imatrices.cam_rgb[0][2] * b;
    wb_mul[1] = imatrices.cam_rgb[1][0] * r + imatrices.cam_rgb[1][1] * g + imatrices.cam_rgb[1][2] * b;
    wb_mul[2] = imatrices.cam_rgb[2][0] * r + imatrices.cam_rgb[2][1] * g + imatrices.cam_rgb[2][2] * b;

    for (int c = 0; c < 3; ++c) {
        wb_mul[c] = ri->get_pre_mul(c) / wb_mul[c] / ref_pre_mul[c];
    }

    // Normalize max channel gain to 1.0
    float mg = rtengine::max(wb_mul[0], wb_mul[1], wb_mul[2]);

    for (int c = 0; c < 3; ++c) {
        wb_mul[c] /= mg;
    }

    return wb_mul;
}


void calcMedians(const Imagefloat* baseImg, const int x1, const int y1, const int x2, const int y2, int skip, std::array<float, 3> &meds)
{
    // Channel vectors to calculate medians
    std::vector<float> rv, gv, bv;

    const int sz = std::max(0, (y2 - y1) * (x2 - x1) / SQR(skip));
    rv.reserve(sz);
    gv.reserve(sz);
    bv.reserve(sz);

#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int i = y1; i < y2; i += skip) {
        for (int j = x1; j < x2; j += skip) {
            rv.push_back(baseImg->r(i, j));
            gv.push_back(baseImg->g(i, j));
            bv.push_back(baseImg->b(i, j));
        }
    }

    // Calculate channel medians from whole image
    findMinMaxPercentile(rv.data(), rv.size(), 0.5f, meds[0], 0.5f, meds[0], true);
    findMinMaxPercentile(gv.data(), gv.size(), 0.5f, meds[1], 0.5f, meds[1], true);
    findMinMaxPercentile(bv.data(), bv.size(), 0.5f, meds[2], 0.5f, meds[2], true);
}


} // namespace

bool RawImageSource::getFilmNegativeExponents(Coord2D spotA, Coord2D spotB, int tran, const procparams::FilmNegativeParams &currentParams, std::array<float, 3>& newExps)
{
    newExps = {
        static_cast<float>(currentParams.redRatio * currentParams.greenExp),
        static_cast<float>(currentParams.greenExp),
        static_cast<float>(currentParams.blueRatio * currentParams.greenExp)
    };

    constexpr int spotSize = 32; // TODO: Make this configurable?

    Coord spot;
    std::array<float, 3> clearVals;
    std::array<float, 3> denseVals;

    // Get channel averages in the two spots, sampling from the original ri->data buffer.
    // NOTE: rawData values might be affected by CA corection, FlatField, etc, so:
    //   rawData[y][x] == (ri->data[y][x] - cblacksom[c]) * scale_mul[c]
    // is not always true. To calculate exponents on the exact values, we should keep
    // a copy of the rawData buffer after preprocessing. Worth the memory waste?

    // Sample first spot
    transformPosition(spotA.x, spotA.y, tran, spot.x, spot.y);

    if (!channelsAvg(ri, W, H, cblacksom, spot, spotSize, clearVals)) {
        return false;
    }

    // Sample second spot
    transformPosition(spotB.x, spotB.y, tran, spot.x, spot.y);

    if (!channelsAvg(ri, W, H, cblacksom, spot, spotSize, denseVals)) {
        return false;
    }

    // Detect which one is the dense spot, based on green channel
    if (clearVals[1] < denseVals[1]) {
        std::swap(clearVals, denseVals);
    }

    if (settings->verbose) {
        printf("Clear film values: R=%g G=%g B=%g\n", static_cast<double>(clearVals[0]), static_cast<double>(clearVals[1]), static_cast<double>(clearVals[2]));
        printf("Dense film values: R=%g G=%g B=%g\n", static_cast<double>(denseVals[0]), static_cast<double>(denseVals[1]), static_cast<double>(denseVals[2]));
    }

    const float denseGreenRatio = clearVals[1] / denseVals[1];

    // Calculate logarithms in arbitrary base
    const auto logBase =
        [](float base, float num) -> float
        {
            return std::log(num) / std::log(base);
        };

    // Calculate exponents for each channel, based on the ratio between the bright and dark values,
    // compared to the ratio in the reference channel (green)
    for (int ch = 0; ch < 3; ++ch) {
        if (ch == 1) {
            newExps[ch] = 1.f;  // Green is the reference channel
        } else {
            newExps[ch] = LIM(logBase(clearVals[ch] / denseVals[ch], denseGreenRatio), 0.3f, 4.f);
        }
    }

    if (settings->verbose) {
        printf("New exponents:  R=%g G=%g B=%g\n", static_cast<double>(newExps[0]), static_cast<double>(newExps[1]), static_cast<double>(newExps[2]));
    }

    return true;
}

bool RawImageSource::getImageSpotValues(Coord2D spotCoord, int spotSize, int tran, const procparams::FilmNegativeParams &params, std::array<float, 3>& rawValues)
{
    Coord spot;
    transformPosition(spotCoord.x, spotCoord.y, tran, spot.x, spot.y);

    if (settings->verbose > 1) {
        printf("Transformed coords: %d,%d\n", spot.x, spot.y);
    }

    if (spotSize < 4) {
        return false;
    }

    // Calculate averages of raw unscaled channels
    if (!channelsAvg(ri, W, H, cblacksom, spot, spotSize, rawValues)) {
        return false;
    }

    if (settings->verbose > 1) {
        printf("Raw spot values: R=%g, G=%g, B=%g\n", rawValues[0], rawValues[1], rawValues[2]);
    }

    return true;
}

void RawImageSource::filmNegativeProcess(const procparams::FilmNegativeParams &params, std::array<float, 3>& filmBaseValues)
{
//    BENCHFUNMICRO

    if (!params.enabled) {
        return;
    }

    // Exponents are expressed as positive in the parameters, so negate them in order
    // to get the reciprocals.
    const std::array<float, 3> exps = {
        static_cast<float>(-params.redRatio * params.greenExp),
        static_cast<float>(-params.greenExp),
        static_cast<float>(-params.blueRatio * params.greenExp)
    };

    constexpr float MAX_OUT_VALUE = 65000.f;

    // Get multipliers for a known, fixed WB setting, that will be the starting point
    // for balancing the converted image.
    const std::array<double, 3> wb_mul = calcWBMults(
            ColorTemp(3500., 1., 1., "Custom"), imatrices, ri, ref_pre_mul);


    if (settings->verbose) {
        printf("Fixed WB mults: %g %g %g\n", wb_mul[0], wb_mul[1], wb_mul[2]);
    }

    // Compensation factor to restore the non-autoWB initialGain (see RawImageSource::load)
    const float autoGainComp = camInitialGain / initialGain;

    std::array<float, 3> mults;  // Channel normalization multipliers

    // If film base values are set in params, use those
    if (filmBaseValues[0] <= 0.f) {
        // ...otherwise, the film negative tool might have just been enabled on this image,
        // whithout any previous setting. So, estimate film base values from channel medians

        std::array<float, 3> medians;

        // Special value for backwards compatibility with profiles saved by RT 5.7
        const bool oldChannelScaling = filmBaseValues[0] == -1.f;

        // If using the old channel scaling method, get medians from the whole current image,
        // reading values from the already-scaled rawData buffer.
        if (oldChannelScaling) {
            calcMedians(ri, rawData, 0, 0, W, H, medians);
        } else {
            // Cut 20% border from medians calculation. It will probably contain outlier values
            // from the film holder, which will bias the median result.
            const int bW = W * 20 / 100;
            const int bH = H * 20 / 100;
            calcMedians(ri, rawData, bW, bH, W - bW, H - bH, medians);
        }

        // Un-scale rawData medians
        for (int c = 0; c < 3; ++c) {
            medians[c] /= scale_mul[c];
        }

        if (settings->verbose) {
            printf("Channel medians: R=%g, G=%g, B=%g\n", medians[0], medians[1], medians[2]);
        }

        for (int c = 0; c < 3; ++c) {

            // Estimate film base values, so that in the output data, each channel
            // median will correspond to 1/24th of MAX.
            filmBaseValues[c] = pow_F(24.f / 512.f, 1.f / exps[c]) * medians[c];

            if (oldChannelScaling) {
                // If using the old channel scaling method, apply WB multipliers here to undo their
                // effect later, since fixed wb compensation was not used in previous version.
                // Also, undo the effect of the autoGainComp factor (see below).
                filmBaseValues[c] /= pow_F((wb_mul[c] / autoGainComp), 1.f / exps[c]);// / autoGainComp;
            }

            // normalize to 0-1
            //filmBaseValues[c] /= 65535.f;
        }
    }


    // Calculate multipliers based on previously obtained film base input values.

    // Apply current scaling coefficients to raw, unscaled base values.
    std::array<float, 3> fb = {
        filmBaseValues[0] * scale_mul[0],
        filmBaseValues[1] * scale_mul[1],
        filmBaseValues[2] * scale_mul[2]
    };

    if (settings->verbose) {
        printf("Input film base values: %g %g %g\n", fb[0], fb[1], fb[2]);
    }

    for (int c = 0; c < 3; ++c) {
        // Apply channel exponents, to obtain the corresponding base values in the output data
        fb[c] = pow_F(max(fb[c], 1.f), exps[c]);

        // Determine the channel multiplier so that the film base value is 1/512th of max.
        mults[c] = (MAX_OUT_VALUE / 512.f) / fb[c];

        // Un-apply the fixed WB multipliers, to reverse their effect later in the WB tool.
        // This way, the output image will be adapted to this fixed WB setting
        mults[c] /= wb_mul[c];

        // Also compensate for the initialGain difference between the default scaling (forceAutoWB=true)
        // and the non-autoWB scaling (see camInitialGain).
        // This effectively restores camera scaling multipliers, and gives us stable multipliers
        // (not depending on image content).
        mults[c] *= autoGainComp;

    }

    if (settings->verbose) {
        printf("Output film base values: %g %g %g\n", static_cast<double>(fb[0]), static_cast<double>(fb[1]), static_cast<double>(fb[2]));
        printf("Computed multipliers: %g %g %g\n", static_cast<double>(mults[0]), static_cast<double>(mults[1]), static_cast<double>(mults[2]));
    }


    constexpr float CLIP_VAL = 65535.f;

    MyTime t1, t2, t3;

    t1.set();

    if (ri->getSensorType() == ST_BAYER) {
#ifdef __SSE2__
        const vfloat onev = F2V(1.f);
        const vfloat clipv = F2V(CLIP_VAL);
#endif

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif

        for (int row = 0; row < H; ++row) {
            int col = 0;
            // Avoid trouble with zeroes, minimum pixel value is 1.
            const float exps0 = exps[FC(row, col)];
            const float exps1 = exps[FC(row, col + 1)];
            const float mult0 = mults[FC(row, col)];
            const float mult1 = mults[FC(row, col + 1)];
#ifdef __SSE2__
            const vfloat expsv = _mm_setr_ps(exps0, exps1, exps0, exps1);
            const vfloat multsv = _mm_setr_ps(mult0, mult1, mult0, mult1);

            for (; col < W - 3; col += 4) {
                STVFU(rawData[row][col], vminf(multsv * pow_F(vmaxf(LVFU(rawData[row][col]), onev), expsv), clipv));
            }

#endif // __SSE2__

            for (; col < W - 1; col += 2) {
                rawData[row][col] = min(mult0 * pow_F(max(rawData[row][col], 1.f), exps0), CLIP_VAL);
                rawData[row][col + 1] = min(mult1 * pow_F(max(rawData[row][col + 1], 1.f), exps1), CLIP_VAL);
            }

            if (col < W) {
                rawData[row][col] = min(mult0 * pow_F(max(rawData[row][col], 1.f), exps0), CLIP_VAL);
            }
        }
    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef __SSE2__
        const vfloat onev = F2V(1.f);
        const vfloat clipv = F2V(CLIP_VAL);
#endif

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif

        for (int row = 0; row < H; row ++) {
            int col = 0;
            // Avoid trouble with zeroes, minimum pixel value is 1.
            const std::array<float, 6> expsc = {
                exps[ri->XTRANSFC(row, 0)],
                exps[ri->XTRANSFC(row, 1)],
                exps[ri->XTRANSFC(row, 2)],
                exps[ri->XTRANSFC(row, 3)],
                exps[ri->XTRANSFC(row, 4)],
                exps[ri->XTRANSFC(row, 5)]
            };
            const std::array<float, 6> multsc = {
                mults[ri->XTRANSFC(row, 0)],
                mults[ri->XTRANSFC(row, 1)],
                mults[ri->XTRANSFC(row, 2)],
                mults[ri->XTRANSFC(row, 3)],
                mults[ri->XTRANSFC(row, 4)],
                mults[ri->XTRANSFC(row, 5)]
            };
#ifdef __SSE2__
            const vfloat expsv0 = _mm_setr_ps(expsc[0], expsc[1], expsc[2], expsc[3]);
            const vfloat expsv1 = _mm_setr_ps(expsc[4], expsc[5], expsc[0], expsc[1]);
            const vfloat expsv2 = _mm_setr_ps(expsc[2], expsc[3], expsc[4], expsc[5]);
            const vfloat multsv0 = _mm_setr_ps(multsc[0], multsc[1], multsc[2], multsc[3]);
            const vfloat multsv1 = _mm_setr_ps(multsc[4], multsc[5], multsc[0], multsc[1]);
            const vfloat multsv2 = _mm_setr_ps(multsc[2], multsc[3], multsc[4], multsc[5]);

            for (; col < W - 11; col += 12) {
                STVFU(rawData[row][col], vminf(multsv0 * pow_F(vmaxf(LVFU(rawData[row][col]), onev), expsv0), clipv));
                STVFU(rawData[row][col + 4], vminf(multsv1 * pow_F(vmaxf(LVFU(rawData[row][col + 4]), onev), expsv1), clipv));
                STVFU(rawData[row][col + 8], vminf(multsv2 * pow_F(vmaxf(LVFU(rawData[row][col + 8]), onev), expsv2), clipv));
            }

#endif // __SSE2__

            for (; col < W - 5; col += 6) {
                for (int c = 0; c < 6; ++c) {
                    rawData[row][col + c] = min(multsc[c] * pow_F(max(rawData[row][col + c], 1.f), expsc[c]), CLIP_VAL);
                }
            }

            for (int c = 0; col < W; col++, c++) {
                rawData[row][col + c] = min(multsc[c] * pow_F(max(rawData[row][col + c], 1.f), expsc[c]), CLIP_VAL);
            }
        }
    } else if (ri->get_colors() == 3) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int row = 0; row < H; ++row) {
            for (int col = 0; col < W; ++col) {
                for (int c = 0; c < 3; ++c) {
                    rawData[row][3*col+c] = min(mults[c] * pow_F(max(rawData[row][3*col+c], 1.f), exps[c]), CLIP_VAL);
                }
            }
        }
    } else if (ri->get_colors() == 1) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int row = 0; row < H; ++row) {
            for (int col = 0; col < W; ++col) {
                rawData[row][col] = min(mults[1] * pow_F(max(rawData[row][col], 1.f), exps[1]), CLIP_VAL);
            }
        }
    }

    t2.set();

    if (settings->verbose > 1) {
        printf("Pow loop time us: %d\n", t2.etime(t1));
    }

    t2.set();

    PixelsMap bitmapBads(W, H);

    int totBP = 0; // Hold count of bad pixels to correct

    if (ri->getSensorType() == ST_BAYER) {
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:totBP) schedule(dynamic,16)
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                if (rawData[i][j] >= MAX_OUT_VALUE) {
                    bitmapBads.set(j, i);
                    ++totBP;
                }
            }
        }

        if (totBP > 0) {
            interpolateBadPixelsBayer(bitmapBads, rawData);
        }

    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:totBP) schedule(dynamic,16)
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                if (rawData[i][j] >= MAX_OUT_VALUE) {
                    bitmapBads.set(j, i);
                    totBP++;
                }
            }
        }

        if (totBP > 0) {
            interpolateBadPixelsXtrans(bitmapBads);
        }
    }

    t3.set();

    if (settings->verbose > 1) {
        printf("Bad pixels count: %d\n", totBP);
        printf("Bad pixels interpolation time us: %d\n", t3.etime(t2));
    }
}


//-----------------------------------------------------------------------------

// *** Non-raw processing ***


bool StdImageSource::getImageSpotValues(Coord2D spotCoord, int spotSize, int tran, const procparams::FilmNegativeParams &params, std::array<float, 3>& out)
{
    if (!imgCopy) {
        imgCopy = new Imagefloat(img->getWidth(), img->getHeight());
        img->getStdImage(wb, 0, imgCopy, PreviewProps(0,0,img->getWidth(), img->getHeight(), 1));
    }

    Coord spot;
    imgCopy->transformPixel(spotCoord.x, spotCoord.y, tran, spot.x, spot.y);

    if (settings->verbose > 1) {
        printf("Transformed coords: %d,%d\n", spot.x, spot.y);
    }

    if (spotSize < 4) {
        return false;
    }

    // Calculate averages of raw unscaled channels
    if (!channelsAvg(imgCopy, imgCopy->getWidth(), imgCopy->getHeight(), spot, spotSize, out)) {
        return false;
    }

    // normalize to 0-1
    // for (auto &v : out) {
    //     v /= 65535.f;
    // }

    if (settings->verbose > 1) {
        printf("Image spot values: R=%g, G=%g, B=%g\n", out[0], out[1], out[2]);
    }

    return true;
}


bool StdImageSource::getFilmNegativeExponents(Coord2D spotA, Coord2D spotB, int tran, const procparams::FilmNegativeParams &currentParams, std::array<float, 3>& newExps)
{
//    using rtengine::Imagefloat;

    if (!imgCopy) {
        imgCopy = new Imagefloat(img->getWidth(), img->getHeight());
        img->getStdImage(wb, 0, imgCopy, PreviewProps(0,0,img->getWidth(), img->getHeight(), 1));
    }


    newExps = {
        static_cast<float>(currentParams.redRatio * currentParams.greenExp),
        static_cast<float>(currentParams.greenExp),
        static_cast<float>(currentParams.blueRatio * currentParams.greenExp)
    };

    constexpr int spotSize = 32; // TODO: Make this configurable?

    Coord spot;
    std::array<float, 3> clearVals;
    std::array<float, 3> denseVals;

    // Sample first spot
    imgCopy->transformPixel(spotA.x, spotA.y, tran, spot.x, spot.y);

    if (!channelsAvg(imgCopy, imgCopy->getWidth(), imgCopy->getHeight(), spot, spotSize, clearVals)) {
        return false;
    }

    // Sample second spot
    imgCopy->transformPixel(spotB.x, spotB.y, tran, spot.x, spot.y);

    if (!channelsAvg(imgCopy, imgCopy->getWidth(), imgCopy->getHeight(), spot, spotSize, denseVals)) {
        return false;
    }

    // Detect which one is the dense spot, based on green channel
    if (clearVals[1] < denseVals[1]) {
        std::swap(clearVals, denseVals);
    }

    if (settings->verbose) {
        printf("Clear film values: R=%g G=%g B=%g\n", clearVals[0], clearVals[1], clearVals[2]);
        printf("Dense film values: R=%g G=%g B=%g\n", denseVals[0], denseVals[1], denseVals[2]);
    }

    const float denseGreenRatio = clearVals[1] / denseVals[1];

    // Calculate logarithms in arbitrary base
    const auto logBase =
        [](float base, float num) -> float
        {
            return std::log(num) / std::log(base);
        };

    // Calculate exponents for each channel, based on the ratio between the bright and dark values,
    // compared to the ratio in the reference channel (green)
    for (int ch = 0; ch < 3; ++ch) {
        if (ch == 1) {
            newExps[ch] = 1.f;  // Green is the reference channel
        } else {
            newExps[ch] = LIM(logBase(clearVals[ch] / denseVals[ch], denseGreenRatio), 0.3f, 4.f);
        }
    }

    if (settings->verbose) {
        printf("New exponents:  R=%g G=%g B=%g\n", newExps[0], newExps[1], newExps[2]);
    }

    return true;
}

void rtengine::StdImageSource::filmNegativeProcess(const procparams::FilmNegativeParams &params, std::array<float, 3>& filmBaseValues)
{
    if (!params.enabled) {
        // If filmneg is not enabled, restore the copy as main image.
        if (imgCopy) {
            if(img) {
                img->allocate(0,0);
                delete img;
            }
            img = imgCopy;
            imgCopy = nullptr;
        }

        return;
    }

    if (params.enabled && !imgCopy) {
        //printf("ONCE!!!!!!\n");
        imgCopy = new Imagefloat(img->getWidth(), img->getHeight());
        img->getStdImage(wb, 0, imgCopy, PreviewProps(0,0,img->getWidth(), img->getHeight(), 1));
    }

    // Destroy old buffer
    img->allocate(0,0);
    delete img;

    // Overwrite working buffer with a copy of the float version
    Imagefloat* posImg = imgCopy->copy();
    img = posImg;


    std::array<float, 3> exps = {
        float(-(params.greenExp * params.redRatio)), // 2.2f;
        float(-params.greenExp),  // 2.2f;
        float(-(params.greenExp * params.blueRatio)) // 2.2f;
    };

    std::array<float, 3> mults;  // Channel normalization multipliers
    constexpr float MAX_OUT_VALUE = 65000.f;

    const int W = imgCopy->getWidth();
    const int H = imgCopy->getHeight();

    // If film base values are set in params, use those
    if (filmBaseValues[0] <= 0.f) {
        // ...otherwise, the film negative tool might have just been enabled on this image,
        // whithout any previous setting. So, estimate film base values from channel medians

        std::array<float, 3> medians;

        // Cut 20% border from medians calculation. It will probably contain outlier values
        // from the film holder, which will bias the median result.
        const int bW = W * 20 / 100;
        const int bH = H * 20 / 100;
        calcMedians(imgCopy, bW, bH, W - bW, H - bH, 4, medians);

        for (int c = 0; c < 3; ++c) {

            // Estimate film base values, so that in the output data, each channel
            // median will correspond to 1/24th of MAX.
            filmBaseValues[c] = pow_F(24.f / 512.f, 1.f / exps[c]) * medians[c];
        }

        if (settings->verbose) {
            printf("Channel medians: R=%g, G=%g, B=%g\n", medians[0], medians[1], medians[2]);
            printf("Estimated base values: R=%g, G=%g, B=%g\n", filmBaseValues[0], filmBaseValues[1], filmBaseValues[2]);
        }
    }

    std::array<float, 3> fb;

    for (int c = 0; c < 3; ++c) {
        // Apply channel exponents, to obtain the corresponding base values in the output data
        fb[c] = pow_F(std::max(filmBaseValues[c], 1.f), exps[c]);

        // Determine the channel multiplier so that the film base value is 1/512th of max.
        mults[c] = (MAX_OUT_VALUE / 512.f) / fb[c];
    }

    if (settings->verbose) {
        printf("Output film base values: %g %g %g\n", static_cast<double>(fb[0]), static_cast<double>(fb[1]), static_cast<double>(fb[2]));
        printf("Computed multipliers: %g %g %g\n", static_cast<double>(mults[0]), static_cast<double>(mults[1]), static_cast<double>(mults[2]));
    }
    const float rmult = mults[0];
    const float gmult = mults[1];
    const float bmult = mults[2];

    const float rexp = exps[0];
    const float gexp = exps[1];
    const float bexp = exps[2];


#ifdef __SSE2__
    const vfloat clipv = F2V(MAXVALF);
    const vfloat rexpv = F2V(rexp);
    const vfloat gexpv = F2V(gexp);
    const vfloat bexpv = F2V(bexp);
    const vfloat rmultv = F2V(rmult);
    const vfloat gmultv = F2V(gmult);
    const vfloat bmultv = F2V(bmult);
#endif

    const int rheight = imgCopy->getHeight();
    const int rwidth = imgCopy->getWidth();

    for (int i = 0; i < rheight; i++) {
        float *rlinein = imgCopy->r(i);
        float *glinein = imgCopy->g(i);
        float *blinein = imgCopy->b(i);
        float *rlineout = posImg->r(i);
        float *glineout = posImg->g(i);
        float *blineout = posImg->b(i);
        int j = 0;
#ifdef __SSE2__

        for (; j < rwidth - 3; j += 4) {
            STVFU(rlineout[j], vminf(rmultv * pow_F(LVFU(rlinein[j]), rexpv), clipv));
            STVFU(glineout[j], vminf(gmultv * pow_F(LVFU(glinein[j]), gexpv), clipv));
            STVFU(blineout[j], vminf(bmultv * pow_F(LVFU(blinein[j]), bexpv), clipv));
        }

#endif

        for (; j < rwidth; ++j) {
            rlineout[j] = CLIP(rmult * pow_F(rlinein[j], rexp));
            glineout[j] = CLIP(gmult * pow_F(glinein[j], gexp));
            blineout[j] = CLIP(bmult * pow_F(blinein[j], bexp));
        }
    }
}

} // namespace rtengine
