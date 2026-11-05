//
//  DSTbackground.h
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under PolyForm Noncommercial 1.0.0
//  You may use, share and adapt this code for noncommercial purposes.
//  Commercial use requires a separate licence from the author.
//  Licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/

#ifndef __DST_background__
#define __DST_background__

#include <memory>
#include <vector>
#include <cstddef>

#include <DSTfits/FITSdata.h>
#include <DSTfits/FITSimg.h>

using namespace DSL;

namespace DST { namespace Background {

/**
 * @class skyEstimator DSTools/DSTbackground.h
 * @brief Estimates the sky background and noise map of a 2-D astronomical image.
 *
 * The image is tessellated into square tiles.  For each tile the sigma-clipped
 * median (sky) and NMAD×1.4826 (noise RMS) are computed from the unmasked pixels.
 * The per-tile values are then bilinearly interpolated onto the full image grid to
 * produce smooth sky and RMS maps.
 *
 * All parameters have sensible defaults and are fully configurable through setters.
 * The class operates in double precision internally and accepts any FITScube pixel
 * type via DoubleValueAtPixel().
 *
 * @author GILLARD William
 */
class skyEstimator
{
public:

    /// Per-tile sky / noise statistics.
    struct TileStat
    {
        double sky   = 0.0;   ///< Sigma-clipped median sky value.
        double rms   = 0.0;   ///< Sigma-clipped NMAD * 1.4826 (noise RMS).
        double cx    = 0.0;   ///< Tile centre x [pixels].
        double cy    = 0.0;   ///< Tile centre y [pixels].
        bool   valid = false; ///< False when the tile has too few unmasked pixels.
    };

    skyEstimator();

#pragma region -- modifiers
    /// Side length of the tessellation tiles in pixels (default: 50).
    void setTileSize(size_t s)           { fTileSize = s; }
    /// Number of sigma-clipping iterations (default: 3).
    void setSigmaClipIter(size_t n)      { fSigmaClipIter = n; }
    /// Sigma-clipping rejection threshold (default: 3.0).
    void setSigmaClipThreshold(double s) { fSigmaClipSigma = s; }
    /// Minimum fraction of unmasked pixels required to mark a tile as valid (default: 0.2).
    void setMinValidFraction(double f)   { fMinValidFraction = f; }
#pragma endregion

    /**
     * @brief Run sky and noise estimation on the given 2-D image.
     * @param img  Input FITScube (must be 2-D).  Masked pixels are ignored.
     * @throw FITSexception if \a img is null or not 2-D.
     */
    void estimate(const std::shared_ptr<FITScube>& img);

#pragma region -- accessors
    /**
     * @brief Return the bilinearly interpolated sky map.
     * @throw std::runtime_error if estimate() has not been called yet.
     */
    std::shared_ptr<FITScube> getSkyMap() const;

    /**
     * @brief Return the bilinearly interpolated noise RMS map.
     * @throw std::runtime_error if estimate() has not been called yet.
     */
    std::shared_ptr<FITScube> getRMSMap() const;

    /// Per-tile statistics populated by estimate().
    const std::vector<TileStat>& tileStats() const { return fStats; }
#pragma endregion

private:
    size_t fTileSize         = 50;
    size_t fSigmaClipIter    = 3;
    double fSigmaClipSigma   = 3.0;
    double fMinValidFraction = 0.2;

    std::vector<TileStat>     fStats;
    std::shared_ptr<FITScube> fSkyMap;
    std::shared_ptr<FITScube> fRMSMap;

    /// Bilinear interpolation of per-tile sky/rms values onto the full image grid.
    void _interpolate(const std::shared_ptr<FITScube>& img,
                      size_t nTX, size_t nTY);

    /**
     * @brief Compute the sigma-clipped median and noise RMS of a pixel sample.
     * @param v      Pixel values (modified in place by the clipping).
     * @param nIter  Number of sigma-clipping passes.
     * @param sigma  Rejection threshold in units of NMAD*1.4826.
     * @param outRms Output: final NMAD*1.4826 after all clipping passes.
     * @return       Sigma-clipped median.
     */
    static double _sigmaClippedMedian(std::vector<double>& v,
                                      size_t nIter,
                                      double sigma,
                                      double& outRms);
};

}} // namespace DST::Background

#endif /* defined(__DST_background__) */
