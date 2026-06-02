//
//  DSTbackground.cxx
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution,
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#include <DSTools/DSTbackground.h>

#include <DSTfits/FITSimg.h>
#include <DSTfits/FITSexception.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DST { namespace Background {

// ---------------------------------------------------------------------------
//  Constructor
// ---------------------------------------------------------------------------
skyEstimator::skyEstimator() {}

// ---------------------------------------------------------------------------
//  Accessors
// ---------------------------------------------------------------------------
std::shared_ptr<FITScube> skyEstimator::getSkyMap() const
{
    if(!fSkyMap)
        throw std::runtime_error(
            "DST::Background::skyEstimator::getSkyMap: estimate() has not been called yet.");
    return fSkyMap;
}

std::shared_ptr<FITScube> skyEstimator::getRMSMap() const
{
    if(!fRMSMap)
        throw std::runtime_error(
            "DST::Background::skyEstimator::getRMSMap: estimate() has not been called yet.");
    return fRMSMap;
}

// ---------------------------------------------------------------------------
//  _sigmaClippedMedian
//  Iteratively rejects pixels beyond sigma * NMAD*1.4826 from the median.
//  The input vector is modified (reordered + shrunk) during clipping.
// ---------------------------------------------------------------------------
double skyEstimator::_sigmaClippedMedian(std::vector<double>& v,
                                          size_t               nIter,
                                          double               sigma,
                                          double&              outRms)
{
    if(v.empty()) { outRms = 0.0; return 0.0; }

    for(size_t iter = 0; iter < nIter; ++iter)
    {
        if(v.empty()) break;

        // Median via nth_element (partial sort)
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        double med = v[mid];
        if(v.size() % 2 == 0 && mid > 0)
        {
            // Lower median: max of the left half
            auto it = std::max_element(v.begin(), v.begin() + mid);
            med = (med + *it) * 0.5;
        }

        // NMAD: median of absolute deviations
        std::vector<double> absDev;
        absDev.reserve(v.size());
        for(double x : v) absDev.push_back(std::abs(x - med));
        const size_t mid2 = absDev.size() / 2;
        std::nth_element(absDev.begin(), absDev.begin() + mid2, absDev.end());
        const double nmad = absDev[mid2] * 1.4826;

        // Reject outliers
        const double lo = med - sigma * nmad;
        const double hi = med + sigma * nmad;
        v.erase(std::remove_if(v.begin(), v.end(),
                               [lo, hi](double x){ return x < lo || x > hi; }),
                v.end());
    }

    if(v.empty()) { outRms = 0.0; return 0.0; }

    // Final median
    const size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double med = v[mid];
    if(v.size() % 2 == 0 && mid > 0)
    {
        auto it = std::max_element(v.begin(), v.begin() + mid);
        med = (med + *it) * 0.5;
    }

    // Final NMAD * 1.4826
    std::vector<double> absDev;
    absDev.reserve(v.size());
    for(double x : v) absDev.push_back(std::abs(x - med));
    const size_t mid2 = absDev.size() / 2;
    std::nth_element(absDev.begin(), absDev.begin() + mid2, absDev.end());
    outRms = absDev[mid2] * 1.4826;

    return med;
}

// ---------------------------------------------------------------------------
//  estimate
// ---------------------------------------------------------------------------
void skyEstimator::estimate(const std::shared_ptr<FITScube>& img)
{
    if(!img)
        throw FITSexception(BAD_DIMEN, "DST::Background", "skyEstimator::estimate", "Input FITScube pointer is null.");
    
    if(img->GetDimension() != 2)
        throw FITSexception(BAD_DIMEN, "DST::Background", "skyEstimator::estimate", "Input FITScube must be 2-D.");

    const size_t W   = img->Size(1);
    const size_t H   = img->Size(2);
    const size_t nTX = (W + fTileSize - 1) / fTileSize;
    const size_t nTY = (H + fTileSize - 1) / fTileSize;

    fStats.clear();
    fStats.resize(nTX * nTY);

    for(size_t ty = 0; ty < nTY; ++ty)
    {
        for(size_t tx = 0; tx < nTX; ++tx)
        {
            const size_t x0 = tx * fTileSize;
            const size_t y0 = ty * fTileSize;
            const size_t x1 = std::min(x0 + fTileSize, W);
            const size_t y1 = std::min(y0 + fTileSize, H);

            TileStat& ts = fStats[ty * nTX + tx];
            ts.cx    = (static_cast<double>(x0) + static_cast<double>(x1 - 1)) * 0.5;
            ts.cy    = (static_cast<double>(y0) + static_cast<double>(y1 - 1)) * 0.5;
            ts.valid = false;

            // Collect unmasked pixel values
            std::vector<double> vals;
            vals.reserve((x1 - x0) * (y1 - y0));
            for(size_t y = y0; y < y1; ++y)
                for(size_t x = x0; x < x1; ++x)
                {
                    const size_t idx = y * W + x;
                    if(!img->Masked(idx))
                        vals.push_back(img->DoubleValueAtPixel(idx));
                }

            const size_t totalPix = (x1 - x0) * (y1 - y0);
            if(vals.empty() ||
               static_cast<double>(vals.size()) < fMinValidFraction * static_cast<double>(totalPix))
            {
                ts.sky = 0.0;
                ts.rms = 1.0; // guard: non-zero prevents division by zero downstream
                continue;
            }

            ts.valid = true;
            double rms  = 0.0;
            ts.sky = _sigmaClippedMedian(vals, fSigmaClipIter, fSigmaClipSigma, rms);
            ts.rms = (rms > 0.0) ? rms : 1.0;
        }
    }

    _interpolate(img, nTX, nTY);
}

// ---------------------------------------------------------------------------
//  _interpolate
//  Bilinear interpolation of per-tile sky/rms values onto every image pixel.
//  Uses clamp-to-edge for pixels whose fractional tile coordinate falls
//  outside the grid (borders of the image).
// ---------------------------------------------------------------------------
void skyEstimator::_interpolate(const std::shared_ptr<FITScube>& img,
                                 size_t nTX, size_t nTY)
{
    const size_t W = img->Size(1);
    const size_t H = img->Size(2);

    fSkyMap = std::make_shared<FITSimg<double>>(std::vector<size_t>{W, H});
    fRMSMap = std::make_shared<FITSimg<double>>(std::vector<size_t>{W, H});

    for(size_t y = 0; y < H; ++y)
    {
        for(size_t x = 0; x < W; ++x)
        {
            // Fractional tile index of this pixel (tile centres are at tile_mid)
            // ftx=0 means "at the centre of tile column 0"
            double ftx = (static_cast<double>(x) + 0.5) / static_cast<double>(fTileSize) - 0.5;
            double fty = (static_cast<double>(y) + 0.5) / static_cast<double>(fTileSize) - 0.5;

            // Clamp to valid tile range before splitting into integer + fraction.
            // Without this, negative ftx/fty (pixels in the first half of the edge
            // tile) would have floor(ftx)==-1 while tx0 is clamped to 0, causing
            // the fractional weight wx to be ~0.5 toward the *second* tile instead
            // of 0 — producing visible bands along the bottom and left borders.
            ftx = std::max(0.0, std::min(static_cast<double>(nTX - 1), ftx));
            fty = std::max(0.0, std::min(static_cast<double>(nTY - 1), fty));

            // Lower-left tile index
            const size_t tx0 = static_cast<size_t>(std::floor(ftx));
            const size_t ty0 = static_cast<size_t>(std::floor(fty));
            const size_t tx1 = std::min(tx0 + 1, nTX - 1);
            const size_t ty1 = std::min(ty0 + 1, nTY - 1);

            // Bilinear weights
            const double wx = ftx - std::floor(ftx);
            const double wy = fty - std::floor(fty);

            const TileStat& t00 = fStats[ty0 * nTX + tx0];
            const TileStat& t10 = fStats[ty0 * nTX + tx1];
            const TileStat& t01 = fStats[ty1 * nTX + tx0];
            const TileStat& t11 = fStats[ty1 * nTX + tx1];

            const double sky = (1.0 - wy) * ((1.0 - wx) * t00.sky + wx * t10.sky)
                             +         wy  * ((1.0 - wx) * t01.sky + wx * t11.sky);
            const double rms = (1.0 - wy) * ((1.0 - wx) * t00.rms + wx * t10.rms)
                             +         wy  * ((1.0 - wx) * t01.rms + wx * t11.rms);

            const size_t idx = y * W + x;
            fSkyMap->SetPixelValue(sky,                   idx);
            fRMSMap->SetPixelValue(rms > 0.0 ? rms : 1.0, idx);
        }
    }

    // Propagate mask from input image
    fSkyMap->MaskPixels(img->GetMask());
    fRMSMap->MaskPixels(img->GetMask());
}

}} // namespace DST::Background
