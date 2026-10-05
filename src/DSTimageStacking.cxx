/**
 * @file DSTimageStacking.cxx
 * @brief Implementation of the DST::stack image-stacking class.
 *
 * @details
 * This file implements the full pixel-combination pipeline declared in
 * DSTimageStacking.h.  The two most complex pieces of work are:
 *
 * 1. **WCS-based area-weighted rebinning** (implemented in push_back()).
 *    Each input frame is reproejcted onto the reference pixel grid by computing
 *    the fractional overlap area (in sky/world-coordinate space) between every
 *    output pixel polygon and the candidate input pixel polygons.  This avoids
 *    any interpolation artefacts and correctly handles differences in pixel
 *    scale, orientation, and WCS projection between frames.
 *
 * 2. **Multi-threaded combination** (implemented in the four `_stackBy*`
 *    private helpers).  Each combination algorithm splits the pixel array into
 *    equal-sized chunks and processes them in parallel using `std::thread`,
 *    with one thread per available hardware CPU core.
 *
 * @author GILLARD William
 * @copyright Centre de Physique des Particules de Marseille
 * @par Licence
 * CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
 */

#include <DSTools/DSTimageStacking.h>
#include <DSTmath/DSTgeo.h>
#include <thread>
#include <algorithm>
#include <queue>           // std::queue — BFS frontier
#include <unordered_set>   // std::unordered_set — BFS visited tracking

namespace DST
{
#pragma region - stack class implementation

#pragma region -- private methods
    //double stack::pixelArea(const std::vector< DST::Math::point > &pts) const
    //{
    //    size_t n = pts.size();
    //    // Initialze area
    //    double area = 0.0;
    //
    //    // Calculate value of shoelace formula
    //    size_t j = n - 1;
    //    for (size_t i = 0; i < n; i++)
    //    {
    //        area += (pts[j].X() + pts[i].X()) * (pts[j].Y() - pts[i].Y());
    //        j = i;  // j is previous vertex to i
    //    }
    //
    //    // Return absolute value
    //    return std::abs(area / 2.);
    //}
#pragma endregion

#pragma region -- ctor/dtor
    /**
     * @details
     * Construction proceeds in three steps:
     *
     * 1. **Type coercion.**  If the reference image can be statically cast to
     *    `FITSimg<double>` directly (i.e. its underlying pixel type isn't already
     *    `double`), that pointer is reused.  Otherwise a new `FITSimg<double>`
     *    is allocated and the pixel values are copied one by one via
     *    `DoubleValueAtPixel()`, which performs any necessary type conversion.
     *
     * 2. **Header propagation.**  All FITS header keywords that are *not*
     *    related to the pixel representation (`BITPIX`, `BSCALE`, `BZERO`,
     *    `BLANK`) or to the data-array geometry (`NAXIS*`, `EXTEND`, `PCOUNT`,
     *    `GCOUNT`) are forwarded to the internal cube.  This preserves the WCS
     *    information and science metadata that later steps rely on.
     *
     * 3. **Weight initialisation.**  A single weight of 1.0 is pushed into
     *    `fWeights` for the reference layer.
     */
    stack::stack(const std::shared_ptr< const FITScube>& refImage, StackingMethod method):
        fStack(!refImage || std::dynamic_pointer_cast<const FITSimg<double>>(refImage)? nullptr : std::static_pointer_cast<FITSimg<double>>(refImage->cast<double>()))
        , fWeights({1})
        , fMethod(method)
    {
        if (!refImage)
            throw FITSexception(SHARED_NULPTR,"DST::stack","stack","Reference image pointer is null.");

        if (refImage->GetDimension() != 2)
            throw FITSexception(BAD_DIMEN,"DST::stack","stack","Only 2D images are supported for stacking.");

        if(fStack == nullptr)
        {
            // The reference image is not stored as double — allocate a fresh
            // FITSimg<double> and copy values pixel by pixel via the generic
            // DoubleValueAtPixel() accessor, which handles any internal scaling.
            fStack = std::make_shared< FITSimg<double> >(FITSimg<double>({refImage->Size(1), refImage->Size(2)}));
            for(size_t i = 0; i < refImage->Size(); ++i)
            {
                fStack->SetPixelValue(refImage->DoubleValueAtPixel(i),i);
            }

            // Copy WCS and relevant header keywords from reference image to stack.
            // Low-level type/size keywords are intentionally skipped because the
            // internal representation is always double and has its own NAXIS values.
            for(FITSDictionary::const_iterator it = refImage->HDU().begin(); it != refImage->HDU().end(); ++it)
            {
                if(it->first.find("BITPIX") != std::string::npos ||
                   it->first.find("BSCALE") != std::string::npos ||
                   it->first.find("BZERO" ) != std::string::npos ||
                   it->first.find("BLANK" ) != std::string::npos ||
                   it->first.find("NAXIS" ) != std::string::npos ||
                   it->first.find("EXTEND") != std::string::npos ||
                   it->first.find("PCOUNT") != std::string::npos ||
                   it->first.find("GCOUNT") != std::string::npos ||
                   it->first.find("NAXES" ) != std::string::npos )
                    continue; // skip type/size-specific keywords
                fStack->HDU().ValueForKey(it->first, it->second.value(), it->second.type(), it->second.comment());
            }
            fStack->reLoadWCS(); // rebuild FITSwcs from the copied header keywords
        }

    }
#pragma endregion

#pragma region -- Modifier
            void stack::switchMethod(const StackingMethod& method)
            {
                if (fMethod == method)
                    return; // no change

                fMethod = method;
            }

            /**
             * @details
             * Delegates to push_back(image, 1.0).  See the two-argument overload
             * for the full algorithm description.
             */
            void stack::push_back(const std::shared_ptr< const FITScube>& image)
            {
                push_back(image, 1.0); // default weight of 1.0 for unweighted stacking
            }

            /**
             * @details
             * ### WCS-based area-weighted rebinning
             *
             * The goal is to fill a new layer `iLayer` (same spatial dimensions as
             * the reference image) with flux values from `image` such that the total
             * flux is conserved and the WCS alignment is exact, regardless of
             * differences in pixel scale, orientation, or projection.
             *
             * The algorithm processes each output pixel @f$n@f$ independently:
             *
             * #### Step 1 — Find candidate input pixels
             *
             * The world coordinates (RA, Dec) of the output pixel centre are obtained
             * from the reference-image WCS via `fStack->WorldCoordinates(n, 0)`.
             * Those coordinates are projected back onto the input-image pixel grid to
             * find the seed input pixel.  A BFS flood-fill then expands 4-connected
             * (N/S/E/W) from the seed, gating each candidate on
             * `intersection_area > kOverlapThreshold * s_out`.  This correctly handles
             * arbitrary pixel-scale ratios and dither angles without fixed-radius
             * assumptions.
             *
             * #### Step 2 — Build the output-pixel polygon in sky coordinates
             *
             * The output pixel polygon `out_pixel` is constructed as an 8-vertex
             * polygon whose vertices are the **midpoints** (in sky coordinates)
             * between the output pixel centre and each of its 8 neighbours:
             * @f[
             *   V_k = \frac{\text{WCS}(\text{centre}) + \text{WCS}(\text{neighbour}_k)}{2},
             *   \quad k = 0, \ldots, 7
             * @f]
             * Using sky-coordinate midpoints rather than pixel-coordinate half-pixel
             * offsets correctly accounts for any WCS non-linearities (e.g. gnomonic
             * projection distortions near the image edges).
             *
             * The area @f$S_\text{out}@f$ of `out_pixel` (computed by the
             * `polygon2D::Area()` method via the shoelace formula in sky coordinates)
             * represents the solid angle of the output pixel.
             *
             * #### Step 3 — Build each input-pixel polygon in sky coordinates
             *
             * For each candidate input pixel @f$k@f$, a 4-vertex polygon `in_pixel`
             * is built from the **4 corner world coordinates** of that pixel:
             * @f[
             *   \text{corners} = \{(x \pm 0.5,\; y \pm 0.5)\}_\text{pixel} \to \text{WCS}
             * @f]
             * Working in world coordinates ensures that any distortion in the input
             * WCS is correctly captured.  The area @f$S_k@f$ of `in_pixel` represents
             * the solid angle of the input pixel.
             *
             * #### Step 4 — Intersection area and flux contribution
             *
             * The actual overlap polygon `intersect` is computed via
             * `out_pixel.GetIntersectionWithPolygons(in_pixel)` (Sutherland-Hodgman
             * or equivalent) and its area @f$s_k@f$ is obtained.
             *
             * The flux contribution of input pixel @f$k@f$ to the output pixel @f$n@f$
             * is:
             * @f[
             *   \Delta v_{n} \mathrel{+}= v_k \cdot \frac{s_k}{S_k}
             * @f]
             * where @f$v_k@f$ is the flux of input pixel @f$k@f$ and
             * @f$s_k / S_k@f$ is the fractional overlap (fraction of the input pixel
             * that falls inside the output pixel).  This ratio preserves surface
             * brightness (flux per unit solid angle).
             *
             * A numerical guard discards contributions with
             * @f$s_k < \varepsilon_\text{double}@f$ to prevent rounding noise
             * from accumulating.  If the total covered area
             * @f$\sum_k s_k / S_\text{out} > 1 + 10^{-3}@f$, the output pixel is
             * masked to flag a numerical anomaly.
             *
             * #### Final storage
             *
             * The completed layer `iLayer` is appended to the internal cube via
             * `fStack->AddLayer(iLayer)`, and `weight` is appended to `fWeights`.
             * Pixels that received no flux (either outside the input frame or fully
             * masked) are masked in `iLayer`.
             */
            void stack::push_back(const std::shared_ptr< const FITScube>& image, const double& weight)
            {
                // Point comparison precision: two world-coordinate points are considered
                // identical if they differ by less than 1e-10 (in the WCS unit, typically degrees).
                DST::Math::point::precision = 1e-10;

                if (!image)
                    throw FITSexception(SHARED_NULPTR,"DST::stack","push_back","Input image pointer is null.");

                if (image->GetDimension() != 2)
                    throw FITSexception(BAD_DIMEN,"DST::stack","push_back","Only 2D images are supported for stacking.");

                //if (image->Size(1) != fStack->Size(1) ||
                //    image->Size(2) != fStack->Size(2))
                //    throw FITSexception(BAD_DIMEN,"DST::stack","push_back","Input image dimensions do not match stack dimensions.");

                // Allocate the rebinned layer — same spatial size as the reference image.
                FITSimg<double> iLayer = FITSimg<double>({fStack->Size(1), fStack->Size(2)});

                size_t      pixel_index;
                pixelCoords dbl_pixel(2);
                worldCoords wcs;

                std::vector<size_t> int_pixel     = std::vector<size_t>(3);
                int_pixel[2] = 0; // layer index (0 for the new layer being built)

                // Cache of input pixels that genuinely overlap the current output pixel.
                // Built by BFS flood-fill; holds precomputed polygon geometry so the
                // flux accumulation loop does not need to redo any intersection work.
                struct OverlapCandidate {
                    size_t               pixel_index; // flat index into input image
                    DST::Math::polygon2D in_pixel;    // 8-vertex WCS polygon of input pixel
                    DST::Math::polygon2D intersect;   // out_pixel ∩ in_pixel polygon
                    double               s;           // intersect.Area()
                    double               S;           // in_pixel.Area()
                };

                // Input image dimensions needed for BFS bounds checks.
                const size_t img_width  = image->Size(1);
                const size_t img_height = image->Size(2);

                // BFS cardinal-direction offsets (4-connected, signed for bounds arithmetic).
                static const std::array<std::pair<int,int>, 4> kCardinal =
                    {{{0,+1},{0,-1},{+1,0},{-1,0}}};

                // Minimum overlap fraction (relative to out_pixel area) to include a
                // candidate in the BFS cache.  1e-4 suppresses edge/corner slivers
                // while keeping genuinely small partial overlaps.
                constexpr double kOverlapThreshold = 1e-4;

                size_t n_pixels = fStack->Size(1)*fStack->Size(2);

                for(size_t n = 0; n < n_pixels; n++ )
                {
                    // --- Progress bar (printed at every 5 % boundary) ---
                    double progress= static_cast<double>(n)/static_cast<double>(n_pixels)*100.;

                    if(!( static_cast<size_t>(progress+0.5)%5) && ( progress - static_cast<double>(static_cast<size_t>(progress+0.5)) ) < 1e-3 )
                    {
                        size_t barWidth = 70;

                        std::cout << "[";
                        size_t pos = static_cast<size_t>(barWidth * progress / 100.);
                        for (size_t i = 0; i < barWidth; ++i)
                        {
                            if (i < pos)
                                std::cout << "=";
                            else if (i == pos)
                                std::cout << ">";
                            else
                                std::cout << " ";
                        }
                        std::cout << "] " << size_t(progress+0.5) << " %\r";
                        std::cout.flush();
                    }

                    // -------------------------------------------------------
                    // Step 1a — Get world coordinates of the output pixel centre.
                    // -------------------------------------------------------
                    wcs = fStack->WorldCoordinates(n,0);

                    // Integer pixel coordinates of the output pixel centre (for neighbour enumeration).
                    const pixelCoords center_pix({static_cast<double>(fStack->PixelCoordinates(n)[0]),
                                                  static_cast<double>(fStack->PixelCoordinates(n)[1]),
                                                  0});

                    // `out_pixel` will accumulate the 8 midpoint vertices that define
                    // the output pixel polygon in sky coordinates (Step 2).
                    DST::Math::polygon2D out_pixel=DST::Math::polygon2D();

                    // -------------------------------------------------------
                    // Step 1b — Build out_pixel polygon from 8 WCS midpoints.
                    //   Each vertex is the sky coordinate of the midpoint between
                    //   this pixel centre and one of its 8 grid neighbours.
                    //   Mapping through WorldCoordinates is exact for any WCS
                    //   projection; averaging projected coordinates would not be.
                    // -------------------------------------------------------
                    static const std::array<std::pair<double,double>, 8> kNeighbors =
                        {{{-1.,-1.},{0.,-1.},{+1.,-1.},{-1.,0.},{+1.,0.},{-1.,+1.},{0.,+1.},{+1.,+1.}}};

                    for(auto const& [dx, dy] : kNeighbors)
                    {
                        const pixelCoords mid_pix({center_pix[0] + dx * 0.5,
                                                   center_pix[1] + dy * 0.5, 0.});
                        // WorldCoordinates on a 3-D fStack returns 3 elements
                        // {RA, Dec, layer_val}; strip the layer component.
                        const worldCoords mid_wcs = fStack->WorldCoordinates(mid_pix, 0);
                        const worldCoords krn_wcs({mid_wcs[0], mid_wcs[1]});
                        out_pixel.AddPoint(DST::Math::point(krn_wcs));
                    }

                    // -------------------------------------------------------
                    // Step 2 — Area of the output-pixel polygon (solid angle).
                    //   Computed here so the BFS gate can use it immediately.
                    // -------------------------------------------------------
                    double s_out = out_pixel.Area();

                    // Axis-aligned bounding box of out_pixel in WCS space.
                    // Used as a fast pre-filter before calling the expensive
                    // (and buggy-for-disjoint-pairs) GetIntersectionWithPolygons:
                    // if the AABBs don't overlap there is zero intersection and
                    // there is no point expanding BFS neighbours from that pixel.
                    double out_xmin =  std::numeric_limits<double>::max();
                    double out_xmax = -std::numeric_limits<double>::max();
                    double out_ymin =  std::numeric_limits<double>::max();
                    double out_ymax = -std::numeric_limits<double>::max();
                    for(size_t vi = 0; vi < out_pixel.size(); ++vi)
                    {
                        out_xmin = std::min(out_xmin, out_pixel[vi].X());
                        out_xmax = std::max(out_xmax, out_pixel[vi].X());
                        out_ymin = std::min(out_ymin, out_pixel[vi].Y());
                        out_ymax = std::max(out_ymax, out_pixel[vi].Y());
                    }

                    // -------------------------------------------------------
                    // Step 3 — BFS flood-fill to discover all input pixels that
                    //   genuinely overlap out_pixel.
                    //
                    //   Seed: the input pixel whose centre maps to the output
                    //   pixel centre.  Expand 4-connected (N/S/E/W); gate each
                    //   candidate first on AABB overlap, then on
                    //   intersection_area > kOverlapThreshold*s_out.
                    //   Cache the precomputed geometry so the flux loop below
                    //   does not redo any polygon work.
                    // -------------------------------------------------------
                    std::vector<OverlapCandidate> overlap_cache;
                    overlap_cache.reserve(16);

                    std::unordered_set<size_t> bfs_visited;
                    bfs_visited.reserve(16);

                    std::queue<size_t> bfs_queue;

                    // Seed from the 4 input pixels surrounding the projected
                    // output-pixel centre (floor/ceil in each axis).
                    // Using all 4 rather than just the rounded pixel handles
                    // the case where rounding steps outside the image boundary
                    // (e.g. projected x = 9.5 with img_width = 10 rounds to 10
                    // which is invalid, but floor gives 9 which overlaps).
                    {
                        dbl_pixel = image->World2Pixel({wcs[0], wcs[1]}, 0);
                        const long sx0 = static_cast<long>(std::floor(dbl_pixel[0]));
                        const long sy0 = static_cast<long>(std::floor(dbl_pixel[1]));
                        for(int sdx : {0, 1})
                        for(int sdy : {0, 1})
                        {
                            const long sx = sx0 + sdx;
                            const long sy = sy0 + sdy;
                            if(sx < 0 || sy < 0 ||
                               static_cast<size_t>(sx) >= img_width ||
                               static_cast<size_t>(sy) >= img_height)
                                continue;
                            std::vector<size_t> s_coords = {static_cast<size_t>(sx),
                                                             static_cast<size_t>(sy), 0};
                            const size_t s_idx = image->PixelIndex(s_coords);
                            if(s_idx >= image->Nelements() || image->isMasked(s_idx))
                                continue;
                            if(!bfs_visited.count(s_idx))
                            {
                                bfs_visited.insert(s_idx);
                                bfs_queue.push(s_idx);
                            }
                        }
                    }

                    while(!bfs_queue.empty())
                    {
                        const size_t cur_idx = bfs_queue.front();
                        bfs_queue.pop();

                        // Build in_pixel polygon for this candidate.
                        const size_t cx = image->PixelCoordinates(cur_idx)[0];
                        const size_t cy = image->PixelCoordinates(cur_idx)[1];

                        DST::Math::polygon2D in_pixel = DST::Math::polygon2D();
                        // Corner pixel coordinates (lower-left, lower-right,
                        // upper-right, upper-left) + 4 edge midpoints.
                        std::vector<double> llc({static_cast<double>(cx)-0.5, static_cast<double>(cy)-0.5});
                        std::vector<double> lrc({static_cast<double>(cx)+0.5, static_cast<double>(cy)-0.5});
                        std::vector<double> urc({static_cast<double>(cx)+0.5, static_cast<double>(cy)+0.5});
                        std::vector<double> ulc({static_cast<double>(cx)-0.5, static_cast<double>(cy)+0.5});
                        std::vector<double> U  ({static_cast<double>(cx)    , static_cast<double>(cy)+0.5});
                        std::vector<double> D  ({static_cast<double>(cx)    , static_cast<double>(cy)-0.5});
                        std::vector<double> R  ({static_cast<double>(cx)+0.5, static_cast<double>(cy)    });
                        std::vector<double> L  ({static_cast<double>(cx)-0.5, static_cast<double>(cy)    });
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(llc)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(lrc)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(urc)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(ulc)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(U)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(D)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(R)));
                        in_pixel.AddPoint(DST::Math::point(image->WorldCoordinates(L)));
                        
                        double S = in_pixel.Area();

                        // Fast AABB pre-filter: if the bounding boxes don't
                        // overlap there is definitely no intersection — skip
                        // the expensive (and unreliable for disjoint polygons)
                        // GetIntersectionWithPolygons call entirely.
                        double in_xmin =  std::numeric_limits<double>::max();
                        double in_xmax = -std::numeric_limits<double>::max();
                        double in_ymin =  std::numeric_limits<double>::max();
                        double in_ymax = -std::numeric_limits<double>::max();
                        for(size_t vi = 0; vi < in_pixel.size(); ++vi)
                        {
                            in_xmin = std::min(in_xmin, in_pixel[vi].X());
                            in_xmax = std::max(in_xmax, in_pixel[vi].X());
                            in_ymin = std::min(in_ymin, in_pixel[vi].Y());
                            in_ymax = std::max(in_ymax, in_pixel[vi].Y());
                        }
                        // Compute signed overlap width in each axis.
                        // Positive → genuine overlap; zero or negative → disjoint.
                        // A tolerance larger than WCS floating-point rounding (~1e-12 deg)
                        // but much smaller than any real pixel (~0.001 deg/px) ensures
                        // that edge-sharing pixels (in_xmax ≈ out_xmin up to ~1e-12)
                        // are correctly rejected without discarding real partial overlaps.
                        const double aabb_tol = 100. * DST::Math::point::precision; // ~1e-8 deg
                        const double x_overlap = std::min(in_xmax, out_xmax) - std::max(in_xmin, out_xmin);
                        const double y_overlap = std::min(in_ymax, out_ymax) - std::max(in_ymin, out_ymin);
                        if(x_overlap <= aabb_tol || y_overlap <= aabb_tol)
                            continue;  // disjoint or edge-only — don't call GetIntersectionWithPolygons

                        // Intersection gate.
                        DST::Math::polygon2D intersect = out_pixel.GetIntersectionWithPolygons(in_pixel);
                        double s = intersect.Area();

                        if(s > kOverlapThreshold * s_out)
                        {
                            // Cache precomputed geometry — flux loop reuses this directly.
                            overlap_cache.push_back({cur_idx,
                                                     std::move(in_pixel),
                                                     std::move(intersect),
                                                     s, S});

                            // Expand 4-connected (N/S/E/W) face-adjacent neighbours.
                            // Diagonal pixels with real overlap are reached transitively.
                            for(auto const& [ndx, ndy] : kCardinal)
                            {
                                const long nx = static_cast<long>(cx) + ndx;
                                const long ny = static_cast<long>(cy) + ndy;
                                if(nx < 0 || ny < 0 ||
                                   static_cast<size_t>(nx) >= img_width ||
                                   static_cast<size_t>(ny) >= img_height)
                                    continue;

                                std::vector<size_t> nb_coords = {static_cast<size_t>(nx),
                                                                  static_cast<size_t>(ny), 0};
                                const size_t nb_idx = image->PixelIndex(nb_coords);
                                if(nb_idx >= image->Nelements())
                                    continue;

                                if(image->isMasked(nb_idx))
                                    continue;
                                    
                                if(bfs_visited.count(nb_idx))
                                    continue;

                                bfs_visited.insert(nb_idx);
                                bfs_queue.push(nb_idx);
                            }
                        }
                        // else: no real overlap — do not cache, do not expand.
                    }

                    // No overlapping input pixels → skip this output pixel.
                    if(overlap_cache.empty())
                    {
                        iLayer.MaskPixel(n);
                        continue;
                    }

                    // -------------------------------------------------------
                    // Step 4 — Accumulate area-weighted flux from cached candidates.
                    //   All polygon geometry was computed during BFS; no recomputation
                    //   is needed here.
                    // -------------------------------------------------------
                    // Accumulated intersection area (sanity-check variable).
                    double s_this = 0;

                    for(const auto& cand : overlap_cache)
                    {
                        const double s = cand.s;
                        const double S = cand.S;

                        // Discard numerical noise slivers (different, tighter threshold
                        // than the BFS gate: applied per-candidate after the full BFS).
                        if(s < 1e-6 * s_out)
                            continue;

                        s_this += s;

                        // Fractional overlap: preserves surface brightness.
                        const double ratio = s / S;
#ifdef _DEBUG_
                        if((ratio - 1.) > 100 * DST::Math::point::precision)
                        {
                            std::cerr << "[DBG] n=" << n
                                      << " in_px=" << cand.pixel_index
                                      << " in_px_xy=("
                                      << image->PixelCoordinates(cand.pixel_index)[0] << ","
                                      << image->PixelCoordinates(cand.pixel_index)[1] << ")"
                                      << " S_in=" << std::setprecision(10) << S
                                      << " s_out=" << s_out
                                      << " s_intersect=" << s
                                      << " ratio=" << ratio
                                      << std::endl;
                            if(ratio > 1.01) {
                                auto polyRows = [](const DST::Math::polygon2D& p) {
                                    return p.size() > 0 ? p.size() + 1 : size_t(0);
                                };
                                const size_t nrows = std::max({polyRows(out_pixel),
                                                               polyRows(cand.in_pixel),
                                                               polyRows(cand.intersect)});

                                auto printCell = [](std::ostream& os,
                                                    const DST::Math::polygon2D& p,
                                                    size_t iv)
                                {
                                    if(iv < p.size())
                                        os << std::setw(14) << p[iv].X()
                                           << " | " << std::setw(14) << p[iv].Y();
                                    else if(iv == p.size() && p.size() > 0)
                                        os << std::setw(14) << p[0].X()
                                           << " | " << std::setw(14) << p[0].Y();
                                    else
                                        os << "             - |              -";
                                };

                                std::cerr << std::setprecision(12)
                                          << "[DBG]   out_pixel(" << out_pixel.size() << ")"
                                          << " | in_pixel(" << cand.in_pixel.size() << ")"
                                          << " | intersect(" << cand.intersect.size() << ")\n"
                                          << "[DBG]   x             | y             "
                                             "| x             | y             "
                                             "| x             | y\n";
                                for(size_t iv = 0; iv < nrows; ++iv) {
                                    std::cerr << "        ";
                                    printCell(std::cerr, out_pixel,       iv);
                                    std::cerr << " | ";
                                    printCell(std::cerr, cand.in_pixel,   iv);
                                    std::cerr << " | ";
                                    printCell(std::cerr, cand.intersect,  iv);
                                    std::cerr << "\n";
                                }
                            }
                        }
#endif

                        // Accumulate the flux contribution: Δv = v_input × (s / S_input)
                        double increment = image->DoubleValueAtPixel(cand.pixel_index) * ratio;
                        (*(iLayer.GetData<double>()))[n] += increment;

                    }

                    // Mask any output pixel whose accumulated flux is numerically zero.
                    // This happens when the output pixel is entirely outside the input frame
                    // or when all overlapping input pixels were masked.
                    if(std::abs((*(iLayer.GetData<double>()))[n]) < std::numeric_limits<double>::min())
                        iLayer.MaskPixel(n);
                }

            // Append the rebinned layer to the internal cube and record its weight.
            fStack->AddLayer(iLayer);
            fWeights.push_back(weight);

            // After AddLayer the cube may have become 3-D for the first time.
            // Inject trivial axis-3 WCS keywords so WCSLIB can parse a valid
            // NAXIS=3 header (plain linear axis: layer index with unit scale).
            // This is only strictly necessary on the 2D→3D transition, but
            // writing them every time is harmless and keeps the header consistent.
            fStack->HDU().ValueForKey("CTYPE3", std::string("FRAME"),   DSL::fChar, "Layer index axis");
            fStack->HDU().ValueForKey("CRPIX3", 1.0,               "Reference pixel layer axis");   // FITS: layer 0 is pixel 1, world = layer index
            fStack->HDU().ValueForKey("CRVAL3", 0.0,               "Reference value layer axis");
            fStack->HDU().ValueForKey("CDELT3", 1.0,               "Pixel scale layer axis");

            // Rebuild the live FITSwcs from the updated header. Without this the
            // WCS object remains 2-D and WorldCoordinates() returns (0,0) for
            // every pixel on the second and subsequent push_back calls.
            fStack->reLoadWCS();

            // Print 100 % completion bar.
            double progress= 100.;

            if(!(static_cast<size_t>(progress)%5) && ( progress - static_cast<double>(static_cast<size_t>(progress)) ) < 1e-2 )
            {
                int barWidth = 70;

                std::cout << "[";
                int pos = barWidth * progress/100.;
                for (int i = 0; i < barWidth; ++i)
                {
                    if (i < pos)
                        std::cout << "=";
                    else if (i == pos)
                        std::cout << ">";
                    else
                        std::cout << " ";
                }
                std::cout << "] " << int(progress) << " %"<<std::endl;
            }

            wcs.clear();
            int_pixel.clear();
        }
#pragma endregion

#pragma region -- Accessor
            /**
             * @details
             * Uses a simple switch over the `StackingMethod` enumerator.
             * Returns `"undefined"` as a safe fallback for values that are added to
             * the enum in the future but not yet handled here.
             */
            const std::string stack::methodName(const stack::StackingMethod& m)
            {
                switch (m)
                {
                    case StackingMethod::MEAN:
                        return "mean";
                    case StackingMethod::MEDIAN:
                        return "median";
                    case StackingMethod::WEIGHTED_MEAN:
                        return "weighted mean";
                    case StackingMethod::LSF:
                        return "linear least squares fit";
                    default:
                        break;
                }

                return "undefined";
            }

            const stack::StackingMethod stack::method() const
            {
                return fMethod;
            }

            std::shared_ptr< const FITScube > stack::getStack() const
            {
                return fStack;
            }

            std::shared_ptr< const FITScube > stack::getStackAt(const size_t& iL) const
            {
                if(!fStack)
                    return nullptr;

                return fStack->Layer(iL);
            }

            const size_t stack::size() const
            {
                if(!fStack)
                    return 0;

                return fStack->GetDimension() == 3 ? fStack->Size(3) : 1;
            }
#pragma endregion

#pragma region -- Processing function
            /**
             * @details
             * Allocates a fresh `FITSimg<double>` of the same spatial dimensions as
             * the reference image, copies non-structural header keywords from the
             * internal cube, then dispatches to the appropriate private worker.
             *
             * The internal cube `fStack` is read-only during combination; no layer is
             * modified.  The caller therefore retains the full unmodified cube and may
             * call combine() again with a different method.
             */
            std::shared_ptr< FITSimg<double> > stack::combine()
            {
                if(!fStack || fStack->Size(3) < 1)
                    throw FITSexception(BAD_DIMEN,"DST::stack","combine","No images in stack to combine.");

                std::shared_ptr< FITSimg<double> > result = std::make_shared< FITSimg<double> >(FITSimg<double>({fStack->Size(1), fStack->Size(2)}));

                // Propagate header keywords to the output image, skipping low-level
                // type/geometry keywords that are specific to the cube representation.
                for(FITSDictionary::const_iterator iK = fStack->HDU().begin(); iK != fStack->HDU().end(); iK++)
                {
                    if(iK->first.find("BITPIX") != std::string::npos ||
                       iK->first.find("BSCALE") != std::string::npos ||
                       iK->first.find("BZERO" ) != std::string::npos ||
                       iK->first.find("BLANK" ) != std::string::npos ||
                       iK->first.find("NAXIS" ) != std::string::npos ||
                       iK->first.find("NAXES" ) != std::string::npos ||
                       iK->first.find("EXTEND") != std::string::npos ||
                       iK->first.find("PCOUNT") != std::string::npos ||
                       iK->first.find("GCOUNT") != std::string::npos ||
                       iK->first.find("CDELT3") != std::string::npos ||
                       iK->first.find("CRPIX3") != std::string::npos ||
                       iK->first.find("CRVAL3") != std::string::npos ||
                       iK->first.find("CTYPE3") != std::string::npos)
                        continue; // skip type/size-specific keywords

                    result->HDU().ValueForKey(iK->first, iK->second.value(), iK->second.type(), iK->second.comment());
                }

                switch (fMethod)
                {
                    case StackingMethod::MEAN:
                        _stackByMean(result);
                        break;
                    case StackingMethod::MEDIAN:
                        _stackByMedian(result);
                        break;
                    case StackingMethod::WEIGHTED_MEAN:
                        _stackByWeightedMean(result);
                        break;
                    case StackingMethod::LSF:
                        _stackByLSF(result);
                        break;
                    default:
                        throw FITSexception(VALUE_UNDEFINED,"DST::stack","combine","Unknown stacking method.");
                }

                return result;
            }

#pragma region -- private processing methods
            /**
             * @details
             * ### Parallelisation strategy
             *
             * The pixel array (length = `stride = Nx × Ny`) is divided into
             * `n_threads` contiguous chunks, where `n_threads` is the number of
             * hardware CPU cores reported by `std::thread::hardware_concurrency()`.
             * Each chunk is processed independently by one `std::thread`; no shared
             * mutable state is accessed (each thread writes to a disjoint set of output
             * pixel indices), so no synchronisation primitives are needed.
             *
             * ### Data layout
             *
             * The internal cube stores all layers in a flat `std::valarray<double>`
             * with stride `stride = Nx × Ny`.  To extract all layer values for pixel
             * @f$n@f$, a `std::slice(n, nlayer, stride)` is used:
             * @f[
             *   \text{indices} = \{n,\; n+\text{stride},\; n+2\cdot\text{stride},\; \ldots\}
             * @f]
             * This pattern is identical across all four combination workers.
             *
             * ### Mean computation
             *
             * For pixel @f$n@f$:
             * @f[
             *   \text{result}[n] = \frac{\sum_{l \,:\, \text{mask}[l]=0} v_{l,n}}
             *                           {N_\text{valid}}
             * @f]
             * where the mask convention is `false` = valid, `true` = masked.
             */
            void stack::_stackByMean(std::shared_ptr< FITSimg<double> >& outImg)         const
            {
                const std::valarray<double>& all_data = *(fStack->GetData<double>());
                const std::valarray<bool>&   all_mask =   fStack->GetMask();
                const size_t stride   = fStack->Size(1) * fStack->Size(2);
                const size_t nlayer   = fStack->Size(3);
                const size_t n_pixels = stride;

                const size_t n_threads = std::max(1u, std::thread::hardware_concurrency());
                const size_t chunk     = (n_pixels + n_threads - 1) / n_threads;

                auto worker = [&](size_t begin, size_t end)
                {
                    for(size_t n = begin; n < end; ++n)
                    {
                        // Extract the depth-wise slice for pixel n across all layers.
                        std::valarray<double> pixel_values = all_data[std::slice(n, nlayer, stride)];
                        std::valarray<bool>   pixel_masks  = all_mask[std::slice(n, nlayer, stride)];

                        // Select only unmasked values (!mask → true where valid).
                        double sum   = ((std::valarray<double>) pixel_values[!pixel_masks]).sum();
                        size_t count = ((std::valarray<bool>)   pixel_masks [!pixel_masks]).size();

                        if(count > 0)
                            outImg->SetPixelValue(sum / static_cast<double>(count), n);
                        else
                        {
                            // No valid layer for this pixel: output zero and mask it.
                            outImg->SetPixelValue(0., n);
                            outImg->MaskPixel(n);
                        }
                    }
                };

                std::vector<std::thread> threads;
                threads.reserve(n_threads);

                for(size_t t = 0; t < n_threads; ++t)
                    threads.emplace_back(worker, t * chunk, std::min((t + 1) * chunk, n_pixels));
                for(auto& th : threads)
                    th.join();
            }

            /**
             * @details
             * ### Median computation
             *
             * For pixel @f$n@f$, the unmasked layer values are collected into a
             * temporary `std::vector<double>`.  `std::nth_element` is used to find
             * the median in O(N) average time without fully sorting the vector:
             *
             * - **Odd count N:** the element at position @f$\lfloor N/2 \rfloor@f$
             *   after `nth_element` is the median.
             * - **Even count N:** the median is the average of the element at position
             *   @f$N/2@f$ and the maximum of the lower half
             *   @f$\{v_0, \ldots, v_{N/2-1}\}@f$:
             *   @f[
             *     \text{median} = \frac{v_{N/2} + \max(v_0,\ldots,v_{N/2-1})}{2}
             *   @f]
             *   (After `nth_element`, the element at `mid` is in its sorted position
             *   and all elements before it are ≤ it, so `std::max_element` over
             *   `[begin, begin+mid)` gives the lower-half maximum.)
             *
             * Parallelisation follows the same thread-per-chunk scheme as
             * `_stackByMean`.
             */
            void stack::_stackByMedian(std::shared_ptr< FITSimg<double> >& outImg)       const
            {
                const std::valarray<double>& all_data = *(fStack->GetData<double>());
                const std::valarray<bool>&   all_mask =   fStack->GetMask();
                const size_t stride   = fStack->Size(1) * fStack->Size(2);
                const size_t nlayer   = fStack->Size(3);
                const size_t n_pixels = stride;

                const size_t n_threads = std::max(1u, std::thread::hardware_concurrency());
                const size_t chunk     = (n_pixels + n_threads - 1) / n_threads;

                auto worker = [&](size_t begin, size_t end)
                {
                    for(size_t n = begin; n < end; ++n)
                    {
                        std::valarray<double> pixel_values = all_data[std::slice(n, nlayer, stride)];
                        std::valarray<bool>   pixel_masks  = all_mask[std::slice(n, nlayer, stride)];

                        // Collect the unmasked layer values into a contiguous vector.
                        std::valarray<double> unmasked = pixel_values[!pixel_masks];
                        const size_t count = unmasked.size();

                        if(count > 0)
                        {
                            std::vector<double> v(std::begin(unmasked), std::end(unmasked));
                            const size_t mid = count / 2;

                            // Partially sort so that v[mid] is in its correct sorted position.
                            std::nth_element(v.begin(), v.begin() + mid, v.end());

                            // For even counts, average v[mid] with the maximum of the lower half.
                            double median = (count % 2 != 0) ? v[mid] : (v[mid] + *std::max_element(v.begin(), v.begin() + mid)) / 2.;
                            outImg->SetPixelValue(median, n);
                        }
                        else
                        {
                            outImg->SetPixelValue(0., n);
                            outImg->MaskPixel(n);
                        }
                    }
                };

                std::vector<std::thread> threads;
                threads.reserve(n_threads);

                for(size_t t = 0; t < n_threads; ++t)
                    threads.emplace_back(worker, t * chunk, std::min((t + 1) * chunk, n_pixels));
                for(auto& th : threads)
                    th.join();
            }

            /**
             * @details
             * ### Weighted mean computation
             *
             * For pixel @f$n@f$, the per-layer weights @f$w_l@f$ stored in `fWeights`
             * are applied to the layer values before summation:
             * @f[
             *   \text{result}[n] =
             *     \frac{\sum_{l \,:\, \text{unmasked}} w_l \, v_{l,n}}
             *          {\sum_{l \,:\, \text{unmasked}} w_l}
             * @f]
             *
             * The weight vector is pre-converted to a `std::valarray<double>` once
             * before the thread loop so that `valarray` arithmetic (element-wise
             * multiplication, masked sum) can be used efficiently inside the worker
             * lambda.
             *
             * Pixels whose total unmasked weight is zero (or numerically below
             * `double` minimum) are masked in the output.
             *
             * Parallelisation follows the same thread-per-chunk scheme as
             * `_stackByMean`.
             */
            void stack::_stackByWeightedMean(std::shared_ptr< FITSimg<double> >& outImg) const
            {
                if(fWeights.size() != fStack->Size(3))
                    throw FITSexception(BAD_DIMEN,"DST::stack","_stackByWeightedMean","Weights size does not match number of layers in stack.");

                const std::valarray<double>& all_data = *(fStack->GetData<double>());
                const std::valarray<bool>&   all_mask =   fStack->GetMask();
                const size_t stride   = fStack->Size(1) * fStack->Size(2);
                const size_t nlayer   = fStack->Size(3);
                const size_t n_pixels = stride;

                // Convert fWeights to a valarray once so we can use valarray masked operations.
                const std::valarray<double> weighted_values(fWeights.data(), fWeights.size());

                const size_t n_threads = std::max(1u, std::thread::hardware_concurrency());
                const size_t chunk     = (n_pixels + n_threads - 1) / n_threads;

                auto worker = [&](size_t begin, size_t end)
                {
                    for(size_t n = begin; n < end; ++n)
                    {
                        std::valarray<double> pixel_values = all_data[std::slice(n, nlayer, stride)];
                        std::valarray<bool>   pixel_masks  = all_mask[std::slice(n, nlayer, stride)];

                        // Element-wise multiply: each layer value is scaled by its weight.
                        pixel_values *= weighted_values; // apply weights to pixel values

                        // Sum the weighted values over unmasked layers.
                        double sum   = ((std::valarray<double>) pixel_values    [!pixel_masks]).sum();
                        // Normalise by the sum of unmasked weights (not the count) so that
                        // frames with higher weights contribute proportionally more.
                        double count = ((std::valarray<double>) weighted_values [!pixel_masks]).sum();

                        if(count > 0)
                            outImg->SetPixelValue(sum / count, n);
                        else
                        {
                            outImg->SetPixelValue(0., n);
                            outImg->MaskPixel(n);
                        }
                    }
                };

                std::vector<std::thread> threads;
                threads.reserve(n_threads);

                for(size_t t = 0; t < n_threads; ++t)
                    threads.emplace_back(worker, t * chunk, std::min((t + 1) * chunk, n_pixels));
                for(auto& th : threads)
                    th.join();


            }

            /**
             * @details
             * ### Linear least-squares fit to the cumulative flux
             *
             * This method provides a robust estimate of the per-layer mean flux by
             * fitting a line @f$C_l = a \cdot l + b@f$ to the running cumulative
             * flux @f$C_l@f$ as a function of layer index @f$l@f$.  The slope
             * @f$a@f$ is the mean flux increment per layer and is written to the
             * output image.
             *
             * **Why cumulative flux?**
             * Fitting to the cumulative sum rather than to the raw flux values makes
             * the estimator equivalent to computing the slope of the integrated signal,
             * which suppresses the contribution of individual outlier frames while
             * retaining sensitivity to the underlying persistent signal.
             *
             * #### Per-pixel algorithm
             *
             * 1. **Cumulative sum (with masking):**
             *    @f[
             *      C_l = \sum_{k=0}^{l} v_{k,n} \cdot \mathbb{1}[\text{unmasked}(k,n)]
             *    @f]
             *    Masked layers contribute 0 to the running total so that the
             *    cumulative curve remains monotonically related to the number of
             *    unmasked frames.
             *
             * 2. **Least-squares sums over unmasked layers @f$l@f$:**
             *    @f[
             *      S_x = \sum l, \quad S_y = \sum C_l, \quad
             *      S_{xx} = \sum l^2, \quad S_{xy} = \sum l \, C_l, \quad N
             *    @f]
             *    Only layers for which the mask is `false` contribute to these sums;
             *    this prevents masked frames from biasing the fit.
             *
             * 3. **Slope (mean flux per layer):**
             *    @f[
             *      a = \frac{N \, S_{xy} - S_x \, S_y}{N \, S_{xx} - S_x^2}
             *    @f]
             *    The denominator @f$D = N \, S_{xx} - S_x^2@f$ is the total sum of
             *    squared deviations of the layer indices from their mean.  If
             *    @f$|D| < \varepsilon_\text{double}@f$ (all unmasked layers have
             *    the same index, which can only happen when @f$N = 1@f$), or if
             *    @f$N < 2@f$, the pixel is masked in the output.
             *
             * Parallelisation follows the same thread-per-chunk scheme as
             * `_stackByMean`.
             */
            void stack::_stackByLSF(std::shared_ptr< FITSimg<double> >& outImg) const
            {
                const std::valarray<double>& all_data = *(fStack->GetData<double>());
                const std::valarray<bool>&   all_mask =   fStack->GetMask();
                const size_t stride   = fStack->Size(1) * fStack->Size(2);
                const size_t nlayer   = fStack->Size(3);
                const size_t n_pixels = stride;

                const size_t n_threads = std::max(1u, std::thread::hardware_concurrency());
                const size_t chunk     = (n_pixels + n_threads - 1) / n_threads;

                auto worker = [&](size_t begin, size_t end)
                {
                    for(size_t n = begin; n < end; ++n)
                    {
                        std::valarray<double> pixel_values = all_data[std::slice(n, nlayer, stride)];
                        std::valarray<bool>   pixel_masks  = all_mask[std::slice(n, nlayer, stride)];

                        // Step 1 — Build cumulative flux; masked layers contribute 0 to the increment.
                        std::vector<double> cum(nlayer, 0.);
                        cum[0] = pixel_masks[0] ? 0. : pixel_values[0];
                        for(size_t iL = 1; iL < nlayer; ++iL)
                            cum[iL] = cum[iL-1] + (pixel_masks[iL] ? 0. : pixel_values[iL]);

                        // Step 2 — Accumulate least-squares sums over unmasked layers only.
                        //   x = layer index (0-based), y = cumulative flux at that index.
                        double sx = 0., sy = 0., sxx = 0., sxy = 0.;
                        size_t N = 0;
                        for(size_t iL = 0; iL < nlayer; ++iL)
                        {
                            if(pixel_masks[iL]) continue; // skip masked layers
                            const double x = static_cast<double>(iL);
                            const double y = cum[iL];
                            sx  += x;
                            sy  += y;
                            sxx += x * x;
                            sxy += x * y;
                            ++N;
                        }

                        // Need at least 2 unmasked layers to fit a line.
                        if(N < 2)
                        {
                            outImg->SetPixelValue(0., n);
                            outImg->MaskPixel(n);
                            continue;
                        }

                        // Step 3 — Compute the slope a = (N·Sxy − Sx·Sy) / (N·Sxx − Sx²).
                        const double dN    = static_cast<double>(N);
                        const double denom = dN * sxx - sx * sx;

                        // Guard against a degenerate denominator (all unmasked layers at the same index).
                        if(std::abs(denom) < std::numeric_limits<double>::min())
                        {
                            outImg->SetPixelValue(0., n);
                            outImg->MaskPixel(n);
                            continue;
                        }

                        // The slope is the estimated mean flux per additional layer.
                        outImg->SetPixelValue((dN * sxy - sx * sy) / denom, n);
                    }
                };

                std::vector<std::thread> threads;
                threads.reserve(n_threads);
                for(size_t t = 0; t < n_threads; ++t)
                    threads.emplace_back(worker, t * chunk, std::min((t + 1) * chunk, n_pixels));
                for(auto& th : threads)
                    th.join();
            }

#pragma endregion

#pragma endregion

#pragma region -- dump
            /**
             * @details
             * Prints a human-readable summary of the stack's contents and configuration
             * to the provided output stream.  The format is not guaranteed to be stable
             * across versions and should not be parsed by external code; it is intended
             * for quick debugging and inspection by developers.
             */
            void stack::PrintInfo(std::ostream& os) const
            {
                os << "DST::stack dump:" << std::endl;
                os << "  Stacking method: " << methodName(fMethod) << std::endl;
                os << "  Number of layers: " << fStack->Size(3) << std::endl;
                os << "  Image dimensions: " << fStack->Size(1) << " x " << fStack->Size(2) << std::endl;
                os << "  Weights: [";
                for(size_t i = 0; i < fWeights.size(); ++i)
                {
                    os << fWeights[i];
                    if(i < fWeights.size() - 1)
                        os << ", ";
                }
                os << "]" << std::endl;
            }

#pragma endregion
}
