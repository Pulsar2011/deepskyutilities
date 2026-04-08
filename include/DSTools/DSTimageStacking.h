/**
 * @file DSTimageStacking.h
 * @brief WCS-aware image stacking for 2-D FITS frames.
 *
 * @details
 * Declares the `DST::stack` class, which accumulates a set of 2-D FITS images
 * into an internal cube and then combines them into a single output image using
 * one of four statistical estimators:
 *
 * | Method              | Description                                         |
 * |---------------------|-----------------------------------------------------|
 * | `MEAN`              | Arithmetic mean over unmasked layers                |
 * | `MEDIAN`            | Median over unmasked layers                         |
 * | `WEIGHTED_MEAN`     | Weighted arithmetic mean (user-supplied weights)    |
 * | `LSF`               | Slope of a linear least-squares fit to the running  |
 * |                     | cumulative flux (robust mean-flux estimator)        |
 *
 * Before co-adding, every input frame is **rebinned** onto the pixel grid of the
 * reference image using the WCS metadata embedded in the FITS headers.  The
 * rebinning is performed pixel-by-pixel by computing the fractional overlap area
 * (in world-coordinate space) between each output pixel polygon and the
 * candidate input pixel polygons, so it is exact to machine precision and
 * correctly handles images with different pixel scales, orientations, or
 * projections.
 *
 * @author GILLARD William
 * @copyright Centre de Physique des Particules de Marseille
 * @par Licence
 * CC BY-NC 4.0 — https://creativecommons.org/licenses/by-nc/4.0/
 *
 * @see DSTimageStacking.cxx for full algorithm details.
 */

#ifndef __DST_background__
#define __DST_background__

#include <memory>
#include <vector>
#include <cstddef>

#include <DSTfits/FITSdata.h>
#include <DSTfits/FITSimg.h>
#include <DSTmath/DSTmath.h>
#include <DSTmath/DSTgeo.h>

using namespace DSL;

namespace DST
{
#pragma region - stack class definition
    /**
     * @class stack
     * @brief Combines multiple pre-aligned 2-D FITS images into a single stacked image.
     *
     * @details
     * `stack` performs pixel-level co-addition of a set of 2-D astronomical images
     * stored as FITS files.  The class supports four combination strategies
     * (mean, median, weighted mean, linear-least-squares fit) and produces a single
     * output image of type `FITSimg<double>`.
     *
     * ### Assumptions and pre-conditions
     * - **Alignment is the caller's responsibility.**  `stack` does not perform any
     *   geometric registration or PSF matching between the input images.  All images
     *   passed to push_back() must already be spatially aligned (i.e. share a
     *   coherent WCS) before being added to the stack.
     * - **WCS-based rebinning.**  Before co-adding, each input image is rebinned onto
     *   the pixel grid of the reference image using the World Coordinate System (WCS)
     *   metadata embedded in the FITS headers.  This ensures that images acquired with
     *   different pixel scales or camera orientations are projected onto a common grid
     *   prior to combination.
     * - **Reference image.**  The reference image is supplied to the constructor and
     *   defines the output pixel grid, pixel scale, and WCS frame.  It becomes the
     *   first layer of the internal data cube.
     *
     * ### Typical workflow
     * @code{.cpp}
     * // 1. Load the reference image (defines the output grid)
     * std::shared_ptr<FITScube> ref ; //Use DST::FITSmanager to open fits file and load reference image
     *
     * // 2. Construct the stacker (default method: MEAN)
     * DST::stack stacker(ref, DST::stack::StackingMethod::MEDIAN);
     *
     * // 3. Push additional frames (same sky region, any pixel scale)
     * for (auto& img : frames)
     *     stacker.push_back(img);          // or push_back(img, weight)
     *
     * // 4. Combine and retrieve the result
     * auto result = stacker.combine();
     * @endcode
     *
     * @note Only 2-D images are supported; 3-D cubes will throw `FITSexception`.
     *
     * @see FITSimg, FITScube, StackingMethod
     */
    class stack
    {
        public:
            // ---------------------------------------------------------------
            /**
             * @enum StackingMethod
             * @brief Selects the pixel-combination algorithm used by combine().
             *
             * @details
             * The chosen method is applied independently to each pixel position
             * across all layers of the internal data cube.  Masked pixels are
             * excluded from every computation.
             */
            enum class StackingMethod
            {
                MEAN,           ///< Arithmetic mean of unmasked pixel values across all layers.
                MEDIAN,         ///< Median of unmasked pixel values across all layers.
                WEIGHTED_MEAN,  ///< Weighted arithmetic mean; weights are set via push_back(image, weight).
                LSF             /**< Slope of a linear least-squares fit to the running cumulative
                                 *   flux as a function of layer index.  Equivalent to a
                                 *   robust estimator of the per-layer mean flux.
                                 *   Requires at least 2 unmasked layers per pixel. */
            };

        private:
            /// @name Internal state
            /// @{

            /**
             * @brief Internal data cube storing all rebinned input layers.
             *
             * @details
             * The cube has dimensions @f$(N_x \times N_y \times N_{\text{layers}})@f$.
             * Layer 0 is the reference image (cast to `double` at construction time);
             * subsequent layers are added by push_back().  All values are stored as
             * `double` regardless of the original FITS pixel type.
             */
            std::shared_ptr< FITSimg<double> > fStack;

            /**
             * @brief Per-layer combination weights for `WEIGHTED_MEAN` stacking.
             *
             * @details
             * `fWeights[0]` is always 1.0 (the reference image).  Each call to
             * push_back(image, weight) appends one entry.  If push_back(image) is
             * called without an explicit weight, a default weight of 1.0 is used.
             * The vector must have the same number of entries as layers in `fStack`
             * when combine() is called with `StackingMethod::WEIGHTED_MEAN`.
             */
            std::vector<double> fWeights;

            /**
             * @brief Currently active stacking method.
             *
             * @details Set at construction time and changeable at any point before
             * combine() is called via switchMethod().
             */
            StackingMethod fMethod;

            /// @}

            // ---------------------------------------------------------------
            /// @name Private combination workers
            /// @{

            /**
             * @brief Combine all layers by arithmetic mean.
             *
             * @details
             * For each pixel position @f$n@f$, iterates over the @f$N_l@f$ layers,
             * ignoring masked pixels, and writes
             * @f[
             *   \text{result}[n] = \frac{\sum_{l \,:\, \text{unmasked}} v_{l,n}}{N_{\text{valid}}}
             * @f]
             * where @f$v_{l,n}@f$ is the value of layer @f$l@f$ at pixel @f$n@f$.
             * Pixels with no valid layer are masked in the output.
             * Computation is parallelised across all available hardware threads.
             *
             * @param[out] outImg  Pre-allocated output image of matching dimensions.
             */
            void _stackByMean        ( std::shared_ptr< FITSimg<double> >& outImg) const;

            /**
             * @brief Combine all layers by median.
             *
             * @details
             * For each pixel position @f$n@f$, collects all unmasked layer values
             * into a temporary vector, partially sorts it using `std::nth_element`
             * (O(N) average), and writes:
             * - **Odd count:** the middle element.
             * - **Even count:** the average of the two middle elements.
             *   @f[
             *     \text{result}[n] =
             *     \begin{cases}
             *       v_{\lfloor N/2 \rfloor} & \text{if } N \text{ is odd}\\
             *       \dfrac{v_{N/2-1} + v_{N/2}}{2} & \text{if } N \text{ is even}
             *     \end{cases}
             *   @f]
             * Pixels with no valid layer are masked in the output.
             * Computation is parallelised across all available hardware threads.
             *
             * @param[out] outImg  Pre-allocated output image of matching dimensions.
             */
            void _stackByMedian      ( std::shared_ptr< FITSimg<double> >& outImg) const;

            /**
             * @brief Combine all layers by weighted arithmetic mean.
             *
             * @details
             * For each pixel position @f$n@f$, applies per-layer weights
             * @f$w_l@f$ stored in `fWeights` and writes:
             * @f[
             *   \text{result}[n] =
             *     \frac{\sum_{l \,:\, \text{unmasked}} w_l \, v_{l,n}}
             *          {\sum_{l \,:\, \text{unmasked}} w_l}
             * @f]
             * Pixels whose total unmasked weight is zero are masked in the output.
             * Computation is parallelised across all available hardware threads.
             *
             * @param[out] outImg  Pre-allocated output image of matching dimensions.
             * @throw FITSexception if `fWeights.size() != fStack->Size(3)`.
             */
            void _stackByWeightedMean( std::shared_ptr< FITSimg<double> >& outImg) const;

            /**
             * @brief Combine all layers via linear least-squares fit to cumulative flux.
             *
             * @details
             * This method fits a linear model to the **running cumulative flux** as a
             * function of layer index, which provides a robust estimator of the
             * per-layer mean flux even in the presence of outlier frames.
             *
             * **Algorithm (per pixel @f$n@f$):**
             *
             * 1. Build the cumulative sum over layers, skipping masked ones:
             *    @f[
             *      C_l = \sum_{k=0}^{l} v_{k,n} \cdot \mathbb{1}[\text{unmasked}(k,n)]
             *    @f]
             *
             * 2. Fit @f$C_l = a \cdot l + b@f$ by ordinary least squares over all
             *    unmasked layers @f$l@f$:
             *    @f[
             *      a = \frac{N \sum x_l C_l - \sum x_l \sum C_l}
             *               {N \sum x_l^2 - \left(\sum x_l\right)^2}
             *    @f]
             *    where @f$x_l = l@f$ (layer index, 0-based) and the sums run over
             *    unmasked layers only.
             *
             * 3. The slope @f$a@f$ is the estimated mean flux per layer and is stored
             *    as the output pixel value.
             *
             * Pixels with fewer than 2 unmasked layers, or for which the denominator
             * is numerically zero, are masked in the output.
             * Computation is parallelised across all available hardware threads.
             *
             * @param[out] outImg  Pre-allocated output image of matching dimensions.
             */
            void _stackByLSF         ( std::shared_ptr< FITSimg<double> >& outImg) const;

            /// @}

        public:

            // ---------------------------------------------------------------
            /// @name Construction / destruction
            /// @{

            /**
             * @brief Default constructor.  Produces an empty, unusable stack.
             *
             * @note A default-constructed `stack` cannot be used for stacking until it
             * is assigned from a properly constructed instance.
             */
            stack() = default;

            /**
             * @brief Construct a stack and register the reference image.
             *
             * @details
             * The reference image defines the output pixel grid, WCS frame, and pixel
             * scale to which all subsequent images (added via push_back()) will be
             * rebinned before co-addition.
             *
             * The image is immediately cast to `double` (if not already) and stored as
             * the first layer of the internal data cube `fStack`.  All WCS-related and
             * scientifically meaningful FITS header keywords are copied to the output
             * header; low-level type/size keywords (`BITPIX`, `BSCALE`, `BZERO`,
             * `BLANK`, `NAXIS*`, `EXTEND`, `PCOUNT`, `GCOUNT`) are intentionally
             * omitted because the internal representation always uses `double`.
             *
             * @param refImage  Shared pointer to the reference 2-D FITS image.
             *                  Must be non-null and 2-D (`GetDimension() == 2`).
             * @param method    Combination algorithm to use when combine() is called.
             *                  Defaults to `StackingMethod::MEAN`.
             *
             * @throw FITSexception (`SHARED_NULPTR`) if @p refImage is null.
             * @throw FITSexception (`BAD_DIMEN`)     if @p refImage is not 2-D.
             */
            stack(const std::shared_ptr< const FITScube>& refImage,
                        StackingMethod method = StackingMethod::MEAN);

            /// @brief Default destructor.
            ~stack() = default;

            /// @}

            // ---------------------------------------------------------------
            /// @name Modifiers
            /// @{

            /**
             * @brief Change the combination algorithm for future calls to combine().
             *
             * @details
             * Has no effect on data already stored in the internal cube.
             * If @p method equals the current method, the call is a no-op.
             *
             * @param method  New stacking method.
             */
            void switchMethod(const StackingMethod& method);

            /**
             * @brief Rebin and append an input image with unit weight.
             *
             * @details
             * Equivalent to `push_back(image, 1.0)`.
             * Use this overload when all frames should contribute equally to the stack
             * (i.e. for `MEAN` or `MEDIAN` stacking).
             *
             * @param image  Shared pointer to the 2-D input FITS image.
             *               Must share a coherent WCS with the reference image.
             *
             * @throw FITSexception (`SHARED_NULPTR`) if @p image is null.
             * @throw FITSexception (`BAD_DIMEN`)     if @p image is not 2-D or its
             *                                         dimensions differ from those of
             *                                         the reference image.
             */
            void push_back(const std::shared_ptr< const FITScube>& image);

            /**
             * @brief Rebin and append an input image with an explicit weight.
             *
             * @details
             * The rebinning step reprojects @p image onto the output pixel grid
             * by computing, for each output pixel, the fractional overlap area
             * between the output pixel polygon and the candidate input pixel
             * polygons (see the implementation in DSTimageStacking.cxx for the full
             * algorithm description).  The result is stored as a new layer in the
             * internal cube `fStack`.
             *
             * @p weight is stored in `fWeights` and is only used by the
             * `WEIGHTED_MEAN` stacking method.  For other methods the weight is
             * ignored at combination time, but it is still recorded for
             * informational purposes.
             *
             * @param image   Shared pointer to the 2-D input FITS image.
             * @param weight  Relative weight assigned to this frame
             *                (must be positive; a value of 1.0 is the neutral weight).
             *
             * @throw FITSexception (`SHARED_NULPTR`) if @p image is null.
             * @throw FITSexception (`BAD_DIMEN`)     if @p image is not 2-D or its
             *                                         dimensions differ from those of
             *                                         the reference image.
             */
            void push_back(const std::shared_ptr< const FITScube>& image, const double& weight);

            /// @}

            // ---------------------------------------------------------------
            /// @name Accessors
            /// @{

            /**
             * @brief Return a human-readable name for a stacking method.
             *
             * @param m  A `StackingMethod` enumerator value.
             * @return   A lowercase string such as `"mean"`, `"median"`,
             *           `"weighted mean"`, or `"linear least squares fit"`.
             *           Returns `"undefined"` for unrecognised values.
             */
            static const std::string methodName(const StackingMethod& m);

            /**
             * @brief Return the currently active stacking method.
             * @return The `StackingMethod` enumerator currently stored in `fMethod`.
             */
            const StackingMethod method() const;

            /**
             * @brief Return a read-only pointer to the full internal data cube.
             *
             * @details
             * The returned cube contains one layer per image registered so far
             * (reference image + all push_back() calls).  Each layer has been
             * rebinned to the reference pixel grid and is stored as `double`.
             *
             * @return Shared pointer to the `FITScube` internal data structure,
             *         or `nullptr` if the object was default-constructed.
             */
            std::shared_ptr< const FITScube > getStack() const;

            /**
             * @brief Return a single layer of the internal data cube.
             *
             * @param iL  Zero-based layer index (0 = reference image).
             * @return    Shared pointer to the requested `FITScube` layer,
             *            or `nullptr` if the object has no data.
             *
             * @see FITScube::Layer()
             */
            std::shared_ptr< const FITScube > getStackAt(const size_t& iL) const;

            /**
             * @brief Return the number of layers currently in the internal cube.
             *
             * @details
             * Equals 1 after construction (reference image only) and increments
             * by 1 with every successful push_back() call.
             *
             * @return Number of layers, or 0 if the object is empty.
             */
            const size_t size() const;

            /// @}

            // ---------------------------------------------------------------
            /// @name Processing
            /// @{

            /**
             * @brief Combine all registered layers into a single output image.
             *
             * @details
             * Applies the algorithm selected by `fMethod` across all layers of
             * the internal cube independently for each pixel position, then
             * returns the result.  WCS and non-structural FITS header keywords
             * are propagated from the internal cube to the output image.
             *
             * The internal cube is **not** modified by this call; combine() may
             * be called multiple times (e.g. to compare methods by calling
             * switchMethod() between invocations).
             *
             * @return A newly allocated `FITSimg<double>` with the same spatial
             *         dimensions as the reference image, containing the combined
             *         pixel values.  Pixels that have no valid (unmasked) input
             *         contribution in any layer are masked in the output.
             *
             * @throw FITSexception (`BAD_DIMEN`)      if the internal cube is empty.
             * @throw FITSexception (`VALUE_UNDEFINED`) if `fMethod` holds an
             *                                          unrecognised enumerator value.
             */
            std::shared_ptr< FITSimg<double> > combine();

            /// @}

            // ---------------------------------------------------------------
            /// @name Diagnostics
            /// @{

            /**
             * @brief Print a human-readable summary of the stack to an output stream.
             *
             * @param os  Destination output stream (e.g. `std::cout`).
             */
            void PrintInfo(std::ostream& os) const;

            /// @}
    };
#pragma endregion
}

#endif
