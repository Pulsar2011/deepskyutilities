/**
 * @file gTest_stackingHelpers.h
 * @brief Test helpers for DST::stack unit tests.
 *
 * Provides:
 *   - WCS injection utilities (CDELT form and CD-matrix rotation form)
 *   - Dithered/rotated synthetic frame factory
 *   - Multi-extension FITS output for visual inspection
 *   - Image-quality metrics (NRMS, intensity-weighted centroid)
 *
 * Dithering model:
 *   Stars are at the same sky positions in every frame.  A dither of (dx, dy)
 *   pixels is implemented by shifting CRPIX by (+dx, +dy); consequently a star
 *   that was at reference pixel (xn, yn) appears at (xn+dx, yn+dy) in the
 *   dithered frame.  After WCS-based reprojection, the stacker places it back
 *   at (xn, yn) in the output.
 *
 *   For rotation θ (CCW) the star pixel in the dithered frame is:
 *     p1 = (50.5+dx) + cos(θ)·(xn−50.5) + sin(θ)·(yn−50.5)
 *     p2 = (50.5+dy) − sin(θ)·(xn−50.5) + cos(θ)·(yn−50.5)
 *
 *   CD-matrix convention (FITS WCS Paper II, PC matrix formulation):
 *     CD1_1 = −scale·cos θ,  CD1_2 = +scale·sin θ
 *     CD2_1 = +scale·sin θ,  CD2_2 = +scale·cos θ
 *   det(CD) = −scale² (constant, independent of θ). ✓
 */

#ifndef GTEST_STACKING_HELPERS_H
#define GTEST_STACKING_HELPERS_H

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>

#include <DSTfits/FITSimg.h>
#include <DSTfits/FITSmanager.h>
#include <DSTfits/FITSdata.h>
#include <DSTools/DSTstarFinder.h>

#include "test_utils.h"

namespace test_helpers
{

// ============================================================================
// WCS injection helpers
// ============================================================================

/**
 * @brief Inject a simple TAN WCS (CDELT form, no rotation) into an in-memory
 *        FITSimg.  Calls reLoadWCS() so the WCS is immediately active.
 *
 * @param img       Image to modify (modified in place).
 * @param crpix1    Reference pixel x (1-based FITS convention).
 * @param crpix2    Reference pixel y (1-based FITS convention).
 * @param crval1    RA at reference pixel [degrees].
 * @param crval2    Dec at reference pixel [degrees].
 * @param scale_deg Plate scale [degrees/pixel].
 *                  CDELT1 = −scale_deg (RA decreases with x),
 *                  CDELT2 = +scale_deg (Dec increases with y).
 */
template<typename T>
void injectWCS(DSL::FITSimg<T>& img,
               double crpix1,
               double crpix2,
               double crval1    = 180.0,
               double crval2    = 0.0,
               double scale_deg = 0.001)
{
    img.HDU().ValueForKey("CTYPE1", std::string("RA---TAN"), DSL::fChar, "WCS axis type");
    img.HDU().ValueForKey("CTYPE2", std::string("DEC--TAN"), DSL::fChar, "WCS axis type");
    img.HDU().ValueForKey("CRPIX1", crpix1,        "Reference pixel x");
    img.HDU().ValueForKey("CRPIX2", crpix2,        "Reference pixel y");
    img.HDU().ValueForKey("CRVAL1", crval1,        "RA at reference pixel [deg]");
    img.HDU().ValueForKey("CRVAL2", crval2,        "Dec at reference pixel [deg]");
    img.HDU().ValueForKey("CDELT1", -scale_deg,    "Pixel scale RA [deg/px]");
    img.HDU().ValueForKey("CDELT2",  scale_deg,    "Pixel scale Dec [deg/px]");
    img.reLoadWCS();
}

/**
 * @brief Inject a TAN WCS with field rotation via CD matrix.
 *
 * Calls reLoadWCS() so the WCS is immediately active.
 *
 * @param img       Image to modify (modified in place).
 * @param crpix1    Reference pixel x (1-based FITS convention).
 * @param crpix2    Reference pixel y (1-based FITS convention).
 * @param theta_rad Rotation angle [radians], CCW from North.
 * @param crval1    RA at reference pixel [degrees].
 * @param crval2    Dec at reference pixel [degrees].
 * @param scale_deg Plate scale [degrees/pixel].
 */
template<typename T>
void injectRotatedWCS(DSL::FITSimg<T>& img,
                      double crpix1,
                      double crpix2,
                      double theta_rad,
                      double crval1    = 180.0,
                      double crval2    = 0.0,
                      double scale_deg = 0.001)
{
    const double ct = std::cos(theta_rad);
    const double st = std::sin(theta_rad);

    // CD = PC * diag(CDELT), PC = [[cos θ, -sin θ],[sin θ, cos θ]]
    // CDELT_1 = -scale, CDELT_2 = +scale  →  det(CD) = -scale² (constant)
    const double cd11 = -scale_deg * ct;  // PC1_1 * CDELT_1
    const double cd12 =  scale_deg * st;  // PC1_2 * CDELT_1  (-(-st)*scale)
    const double cd21 =  scale_deg * st;  // PC2_1 * CDELT_2
    const double cd22 =  scale_deg * ct;  // PC2_2 * CDELT_2

    img.HDU().ValueForKey("CTYPE1", std::string("RA---TAN"), DSL::fChar, "WCS axis type");
    img.HDU().ValueForKey("CTYPE2", std::string("DEC--TAN"), DSL::fChar, "WCS axis type");
    img.HDU().ValueForKey("CRPIX1", crpix1,  "Reference pixel x");
    img.HDU().ValueForKey("CRPIX2", crpix2,  "Reference pixel y");
    img.HDU().ValueForKey("CRVAL1", crval1,  "RA at reference pixel [deg]");
    img.HDU().ValueForKey("CRVAL2", crval2,  "Dec at reference pixel [deg]");
    img.HDU().ValueForKey("CD1_1",  cd11,    "CD matrix element [deg/px]");
    img.HDU().ValueForKey("CD1_2",  cd12,    "CD matrix element [deg/px]");
    img.HDU().ValueForKey("CD2_1",  cd21,    "CD matrix element [deg/px]");
    img.HDU().ValueForKey("CD2_2",  cd22,    "CD matrix element [deg/px]");
    img.reLoadWCS();
}


// ============================================================================
// Dithered frame factory
// ============================================================================

/**
 * @brief Create one dithered (and optionally rotated) frame.
 *
 * All stars have the same sky coordinates across every frame; only their pixel
 * positions differ according to the dither and rotation of this frame.
 *
 * Pixel position of star n in this frame (TAN, small-angle, analytic):
 *   p1 = (50.5+dx) + cos(θ)·(xn−50.5) + sin(θ)·(yn−50.5)
 *   p2 = (50.5+dy) − sin(θ)·(xn−50.5) + cos(θ)·(yn−50.5)
 *
 * Stars falling within 5 pixels of any edge are silently skipped.
 *
 * @param catalog   Star catalog expressed in the reference frame
 *                  (CRPIX=(50.5,50.5)).
 * @param width     Image width in pixels.
 * @param height    Image height in pixels.
 * @param dx        Dither x [pixels] (positive → star moves right).
 * @param dy        Dither y [pixels] (positive → star moves up).
 * @param theta_rad Rotation angle [radians], CCW from North.
 * @param crval1    RA at reference pixel [degrees].
 * @param crval2    Dec at reference pixel [degrees].
 * @param scale_deg Plate scale [degrees/pixel].
 * @return Shared pointer to the generated FITSimg<double> with WCS set.
 */
inline std::shared_ptr<DSL::FITSimg<double>>
makeDitheredFrame(const std::vector<DST::Stars>& catalog,
                  size_t width,
                  size_t height,
                  double dx,
                  double dy,
                  double theta_rad  = 0.0,
                  double crval1     = 180.0,
                  double crval2     = 0.0,
                  double scale_deg  = 0.001,
                  double background = 0.0)
{
    auto frame = std::make_shared<DSL::FITSimg<double>>(
                     std::vector<size_t>{width, height});

    // Initialise with constant background
    for (size_t i = 0; i < frame->Nelements(); ++i)
        frame->SetPixelValue(background, i);

    constexpr double ref_crpix1 = 50.5;
    constexpr double ref_crpix2 = 50.5;
    const double ct = std::cos(theta_rad);
    const double st = std::sin(theta_rad);

    for (const auto& star : catalog) {
        const double dxn = star.position().X() - ref_crpix1;
        const double dyn = star.position().Y() - ref_crpix2;

        const double p1 = (ref_crpix1 + dx) + ct * dxn + st * dyn;
        const double p2 = (ref_crpix2 + dy) - st * dxn + ct * dyn;

        // Skip stars too close to the image border
        if (p1 < 5.0 || p1 > static_cast<double>(width)  - 5.0 ||
            p2 < 5.0 || p2 > static_cast<double>(height) - 5.0)
            continue;

        DST::Stars shifted = star;
        shifted.setPosition(DST::Math::point({p1, p2}));
        TestUtils::addGaussianStar(*frame, shifted);
    }

    // Inject WCS
    if (std::abs(theta_rad) < 1e-10)
        injectWCS(*frame, ref_crpix1 + dx, ref_crpix2 + dy,
                  crval1, crval2, scale_deg);
    else
        injectRotatedWCS(*frame, ref_crpix1 + dx, ref_crpix2 + dy,
                         theta_rad, crval1, crval2, scale_deg);

    return frame;
}


// ============================================================================
// Multi-extension FITS output for visual inspection
// ============================================================================

/**
 * @brief Write a 3-extension FITS file for visual inspection.
 *
 * Structure:
 *   Primary HDU   — reference image (@p refImg)
 *   Extension 1   — full stack cube (@p stackCube, one layer per frame)
 *   Extension 2   — combined output  (@p combined)
 *
 * The output directory is created if it does not exist.
 *
 * @param path      Absolute or relative path to the output FITS file.
 * @param refImg    Reference image (first layer / defines the output grid).
 * @param stackCube Internal stack cube returned by DST::stack::getStack().
 * @param combined  Combined output returned by DST::stack::combine().
 */
inline void saveInspectionFITS(const std::string& path,
                                DSL::FITScube&    refImg,
                                DSL::FITScube&    stackCube,
                                DSL::FITScube&    combined)
{
    std::filesystem::path dir = std::filesystem::path(path).parent_path();
    if (!dir.empty())
        std::filesystem::create_directories(dir);
    
    refImg.SetName("REF");
    stackCube.SetName("STACK");

    // Write primary HDU
    refImg.Write(path, /*replace=*/true);

    // Append additional extensions
    DSL::FITSmanager fm(path, /*readOnly=*/false);
    fm.AppendImage(stackCube);
    fm.AppendImage(combined);
}


// ============================================================================
// Image quality metrics
// ============================================================================

/**
 * @brief Compute the normalised RMS difference between two images.
 *
 * NRMS = sqrt( mean( (output − ref)² ) ) / max( |ref| )
 *
 * Masked pixels in either image are excluded from the calculation.
 *
 * @return NRMS value, or NaN if fewer than 2 unmasked pixels exist or the
 *         reference maximum is numerically zero.
 */
inline double computeNRMS(const DSL::FITScube& output,
                           const DSL::FITScube& ref)
{
    const size_t n = output.Nelements();
    if (n != ref.Nelements()) return std::numeric_limits<double>::quiet_NaN();

    const auto& out_mask = output.GetMask();
    const auto& ref_mask = ref.GetMask();

    double sum_sq = 0.0;
    double max_val = 0.0;
    size_t count   = 0;

    for (size_t i = 0; i < n; ++i) {
        if (out_mask.size() > 0 && out_mask[i]) continue;
        if (ref_mask.size() > 0 && ref_mask[i]) continue;

        const double rv = ref.DoubleValueAtPixel(i);
        const double ov = output.DoubleValueAtPixel(i);
        const double d  = ov - rv;
        sum_sq  += d * d;
        max_val  = std::max(max_val, std::abs(rv));
        ++count;
    }

    if (count < 2 || max_val < std::numeric_limits<double>::min())
        return std::numeric_limits<double>::quiet_NaN();

    return std::sqrt(sum_sq / static_cast<double>(count)) / max_val;
}

/**
 * @brief Compute the intensity-weighted centroid within a square window.
 *
 * @param img        2-D image to sample.
 * @param expected_x Expected centroid x in pixel coordinates.
 * @param expected_y Expected centroid y in pixel coordinates.
 * @param half_win   Half-width of the search window (full size = 2*half_win+1).
 * @return {cx, cy} centroid, or {-1.0, -1.0} if total weight is zero.
 */
inline std::pair<double, double>
computeCentroid(const DSL::FITScube& img,
                double expected_x, double expected_y,
                int half_win = 3)
{
    const size_t nx = img.Size(1);
    const size_t ny = img.Size(2);
    const auto& mask = img.GetMask();

    const int x0 = std::max(0, static_cast<int>(std::round(expected_x)) - half_win);
    const int x1 = std::min(static_cast<int>(nx) - 1,
                             static_cast<int>(std::round(expected_x)) + half_win);
    const int y0 = std::max(0, static_cast<int>(std::round(expected_y)) - half_win);
    const int y1 = std::min(static_cast<int>(ny) - 1,
                             static_cast<int>(std::round(expected_y)) + half_win);

    double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0;

    for (int iy = y0; iy <= y1; ++iy) {
        for (int ix = x0; ix <= x1; ++ix) {
            const size_t idx = static_cast<size_t>(ix) +
                               static_cast<size_t>(iy) * nx;
            if (mask.size() > 0 && mask[idx]) continue;

            const double v = img.DoubleValueAtPixel(idx);
            if (v <= 0.0) continue;

            sum_w  += v;
            sum_wx += v * static_cast<double>(ix);
            sum_wy += v * static_cast<double>(iy);
        }
    }

    if (sum_w < std::numeric_limits<double>::min())
        return {-1.0, -1.0};

    return {sum_wx / sum_w, sum_wy / sum_w};
}

} // namespace test_helpers

#endif // GTEST_STACKING_HELPERS_H
