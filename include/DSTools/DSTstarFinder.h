/**
 * @file DSTstarFinder.h
 * @brief Star detection and PSF fitting interface for FITS images.
 * @author GILLARD William — Centre de Physique des Particules de Marseille
 *
 * @copyright Licensed under PolyForm Noncommercial 1.0.0.
 * You may use, share and adapt this code for noncommercial purposes.
 * Commercial use requires a separate licence from the author.
 * Full licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/
 *
 * @details
 * This header exposes three cooperating classes:
 *
 * - **DST::Stars** — a plain-data object that stores the centroid position,
 *   pixel membership, photometric properties (flux, SNR, FWHM, …) and PSF fit
 *   parameters of a single detected star.
 *
 * - **DST::algorithms::starFitter** — a Minuit2 cost function that fits a
 *   parametric PSF profile (Gaussian, Moffat, Lorentzian, Voigt, Sérsic, King,
 *   Airy, or a user-defined model) to a small image stamp centred on a star.
 *
 * - **DST::algorithms::starFinder** — the top-level detection engine.  It runs
 *   one of five detection pipelines (MORPHO, THRESHOLD, SNR, NOISECHISEL,
 *   CUSTOM) on a 2D or 3D `FITScube`, then evaluates basic stellar properties
 *   for each blob using moment analysis or a PSF model fit.  An optional
 *   GPU-accelerated morphology backend is selected automatically at construction
 *   time via `DST::GPU::Backend::create()`.
 *
 * ### Typical usage
 * @code
 * auto sky = manager.GetPrimary();                        // shared_ptr<FITScube>
 * DST::algorithms::starFinder finder(sky);
 * finder.findStars();                                     // default: MORPHO
 *
 * // Refine with Gaussian PSF fit
 * for (auto& [id, star] : finder.detectedStars()) {
 *     finder.EvalStellarProperties(star, DST::Stars::evalType::GAUSS2D);
 * }
 *
 * // Export for SCAMP
 * finder.toFITS_LDAC("catalogue.ldac");
 * @endcode
 */

#ifndef __DST_starFinder__
#define __DST_starFinder__

#include <stdexcept>
#include <vector>
#include <map>
#include <memory>
#include <valarray>
#include <cmath>
#include <algorithm>
#include <thread>
#include <future>
#include <limits>
#include <type_traits>
#include <filesystem>

/**
 * @def DST_HAS_SPAN
 * @brief 1 when std::span is available, 0 otherwise.
 *
 * @details starFitter::Eval() takes its parameters as a view where the standard library
 * provides one, so that a caller holding a longer vector can pass a sub-range without
 * copying. That is a C++20 facility, and this header is public: including <span>
 * unconditionally would make C++20 a hard floor for every consumer of the library, while
 * the build still degrades to C++17 and below.
 *
 * The C++17 fallback takes a const std::vector<double>& instead. Callers that pass a whole
 * vector are unaffected either way -- a vector binds implicitly to std::span<const double>
 * -- so only the sub-range case differs, and there it costs a copy rather than failing to
 * compile.
 */
#if __cplusplus >= 202002L
#  include <span>
#  define DST_HAS_SPAN 1
#else
#  define DST_HAS_SPAN 0
#endif

#include <DSTfits/FITSdata.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITStable.h>
#include <DSTmath/DSM_version.h>
#include <DSTmath/DSTgeo.h>
#include <DSTmath/DSTmath.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>



#include "DSTmorpho.h"
#include "DSTmorph_gpu.h"
#include "DSTbackground.h"

#ifdef __WITH_QT__
#include <DSTdisplay/QtDisplay.h>
#include <DSTdisplay/QCustomContours.h>
#include <DSTdisplay/QtFitsPloting.h>
#endif

using namespace DSL;
using namespace DST::Math;

namespace DST
{

#pragma region - DST::Stars class definition
    /**
     * @brief Plain-data container for a single detected stellar object.
     *
     * @details
     * A `Stars` instance is produced by `starFinder` for each blob that passes
     * the detection criteria.  It bundles:
     * - the list of image pixels belonging to the star (`fpixels`),
     * - derived photometric quantities (centroid, peak, FWHM, total flux, SNR,
     *   local background), and
     * - optionally a full PSF fit parameter vector (`fPSFparameters`) together
     *   with the profile model that was used (`fevalMethod`).
     *
     * All physical quantities are in **pixel** and **image-unit** space.  If
     * WCS information is needed the caller should convert centroid positions
     * using `FITScube::WorldCoordinates()` after detection.
     *
     * Instances are normally created by `starFinder`; direct construction is
     * also supported for unit tests or custom pipelines.
     */
    class Stars
        {
            public:
                /**
                 * @brief PSF profile model used to evaluate stellar properties.
                 *
                 * @details
                 * Passed to `starFinder::EvalStellarProperties()` and stored in
                 * each `Stars` instance via `setEvalMethod()`.  The `MOMENT`
                 * model uses image moments only (no fit); all other values trigger
                 * a Minuit2 minimisation via `starFitter`.  `UNDEFINED` is a
                 * sentinel meaning "no evaluation requested yet".
                 */
                enum class evalType
                {
                    MOMENT,       ///< Image-moment centroid and flux (no PSF fit).
                    GAUSS,        ///< Symmetric 2-D Gaussian + constant background.
                    GAUSS2D,      ///< Asymmetric (elliptical) 2-D Gaussian + background.
                    iGAUSS,       ///< Integrated symmetric Gaussian (pixel-convolved).
                    iGAUSS2D,     ///< Integrated asymmetric Gaussian (pixel-convolved).
                    MOFFAT,       ///< Symmetric Moffat profile + background.
                    MOFFAT2D,     ///< Asymmetric Moffat profile + background.
                    LORENTZIAN,   ///< Symmetric Lorentzian (Cauchy) profile + background.
                    LORENTZIAN2D, ///< Asymmetric Lorentzian profile + background.
                    VOIGT,        ///< Symmetric pseudo-Voigt profile + background.
                    VOIGT2D,      ///< Asymmetric pseudo-Voigt profile + background.
                    SERSIC,       ///< Symmetric Sérsic profile + background.
                    SERSIC2D,     ///< Asymmetric Sérsic profile + background.
                    KING,         ///< Symmetric King profile + background.
                    KING2D,       ///< Asymmetric King profile + background.
                    AIRY,         ///< Symmetric Airy disk profile + background.
                    CUSTOM,       ///< User-defined method implemented in `starFinder::_evalStellarPropertiesCustom()`.
                    CUSTOM_MODEL, ///< User-defined method with custom PSF model implemented in `starFinder::_evalStellarPropertiesCustomModel()`.
                    UNDEFINED     ///< Sentinel — no evaluation performed yet.
                };

                static ROOT::Minuit2::MnUserParameters generateLabelParameters(const evalType& method);

            private:
                const size_t fID;                       ///< Unique star identifier, assigned at construction.
                evalType fevalMethod;                   ///< Profile model used in the last `EvalStellarProperties()` call.
                std::vector<point> fpixels;             ///< Pixels classified as belonging to this star.
                point fposition;                        ///< Centroid position [pixels].
                std::vector<double> fPSFparameters;     ///< Best-fit PSF parameter vector; layout depends on `fevalMethod`.

                double fMaxima;          ///< Peak pixel value in the stellar region [image units].
                point fMaximaPosition;   ///< Position of the peak pixel [pixels].
                point fFWHM;             ///< FWHM along each axis (x, y) [pixels].
                double fangle;           ///< PSF orientation angle [radians].
                double fTotalFlux;       ///< Integrated flux over the stellar region [image units].
                double fSNR;             ///< Signal-to-noise ratio (peak / RMS background).
                double fBackground;      ///< Local background level [image units].
                std::string finURL;      ///< Source FITS file path (set by `starFinder`).

            public:
                /// @brief Default constructor — creates a star with ID 0 and NaN photometry.
                Stars();
                /**
                 * @brief Construct a star with a given identifier.
                 * @param id  Unique star identifier (assigned by `starFinder`).
                 */
                Stars(const size_t& id);
                /**
                 * @brief Construct a star at a given pixel position (vector form).
                 * @param id   Unique star identifier.
                 * @param pos  Centroid pixel coordinates; must have at least 2 elements (x, y).
                 * @throws std::invalid_argument if `pos.size() < 2`.
                 */
                Stars(const size_t& id, const pixelCoords& pos);
                /**
                 * @brief Construct a star at a given pixel position (point form).
                 * @param id   Unique star identifier.
                 * @param pos  Centroid position as a `point`; must have at least 2 dimensions.
                 * @throws std::invalid_argument if `pos.size() < 2`.
                 */
                Stars(const size_t& id, const point& pos);
                /// @brief Destructor — clears pixel list and PSF parameter vector.
                ~Stars();

#pragma region -- Modifiers
                /// @brief Append a single pixel to the star's pixel membership list.
                /// @param pix  Pixel coordinates (must have at least 2 elements).
                /// @throws std::invalid_argument if `pix.size() < 2`.
                void addPixel(const pixelCoords& );
                /// @brief Append multiple pixels from a `pixelVectors` container.
                void addPixel(const pixelVectors& );
                /// @brief Append multiple pixels from a `std::vector<point>` container.
                void addPixel(const std::vector<point>& );
                /// @brief Append a single pixel to the star's pixel membership list.
                /// @param pix  Pixel position as a `point` (must have at least 2 dimensions).
                /// @throws std::invalid_argument if `pix.size() < 2`.
                void addPixel(const point&);

                /// @brief Set the centroid position from a pixel-coordinate vector.
                void setPosition(const pixelCoords& );
                /// @brief Set the centroid position from a `point`.
                void setPosition(const point& );

                /// @brief Set the peak pixel value in the stellar region [image units].
                void setMaxima(const double& );
                /// @brief Set the position of the peak pixel (vector form).
                void setMaximaPosition(const pixelCoords& );
                /// @brief Set the position of the peak pixel (point form).
                void setMaximaPosition(const point& );
                /// @brief Set the FWHM along both axes separately [pixels].
                void setFWHM(const double& fwhmX, const double& fwhmY);
                /// @brief Set the FWHM from a `point` holding (FWHM_x, FWHM_y) [pixels].
                void setFWHM(const point& fwhm);
                /// @brief Set the PSF orientation angle [radians].
                void setAngle(const double& angle);
                /// @brief Set the integrated flux over the stellar region [image units].
                void setTotalFlux(const double& flux);
                /// @brief Set the signal-to-noise ratio.
                void setSNR(const double& snr);
                /// @brief Set the local background level [image units].
                void setBackground(const double& background);
                /// @brief Record which PSF model was used to derive the stored properties.
                void setEvalMethod(const evalType& method);
                /// @brief Store the best-fit PSF parameter vector.
                /// @param params  Parameter vector whose layout matches `evalMethod()`.
                void setPSFparameters(const std::vector<double>& params);
#pragma endregion
#pragma region -- Accessors
                inline const size_t&             ID()             const {return fID;}             ///< Unique star identifier.
                inline const std::vector<point>& FilledPixels()   const {return fpixels;}         ///< Pixels belonging to this star.
                inline const point&              position()       const {return fposition;}        ///< Centroid position [pixels].
                inline const double&             maxima()         const {return fMaxima;}          ///< Peak pixel value [image units].
                inline const point&              maximaPosition() const {return fMaximaPosition;}  ///< Position of the peak pixel.
                inline const point&              FWHM()           const {return fFWHM;}            ///< FWHM (x, y) [pixels].
                inline const double&             angle()          const {return fangle;}           ///< PSF orientation [radians].
                inline const double&             totalFlux()      const {return fTotalFlux;}       ///< Integrated flux [image units].
                inline const double&             SNR()            const {return fSNR;}             ///< Signal-to-noise ratio.
                inline const double&             background()     const {return fBackground;}      ///< Local background [image units].
                inline const evalType&           evalMethod()     const {return fevalMethod;}      ///< PSF model used in the last evaluation.
                inline const std::vector<double>& PSFparameters() const {return fPSFparameters;}  ///< Best-fit PSF parameters.

                /// @brief Return the active `evalType` as a human-readable string (e.g. `"GAUSS-FIT"`).
                const std::string          evalMethodString() const;
                /// @brief Return the string name for an arbitrary `evalType` value.
                /// @param m  Profile model to convert.
                /// @return   Uppercase string such as `"MOFFAT-FIT"` or `"UNKNOWN"`.
                static const std::string   ModelString(const evalType& m);
                /// @brief Return the source FITS file path associated with this star.
                inline const std::string&  srcFile() const {return finURL;}
#pragma endregion
#pragma region -- Dumper
                /**
                 * @brief Print stellar properties to an output stream.
                 *
                 * @param os       Destination stream (default: `std::cout`).
                 * @param _header_ If `true`, emit a `#`-prefixed column header line first.
                 */
                void dump(std::ostream& os = std::cout, bool _header_=true) const;
#pragma endregion
        };
#pragma endregion

    namespace algorithms
    {

#pragma region - DST::starFitter class definition
        /**
         * @brief Minuit2 cost function for fitting a parametric PSF to a stellar image stamp.
         *
         * @details
         * `starFitter` implements `ROOT::Minuit2::FCNBase` so it can be passed directly
         * to any Minuit2 minimiser (e.g. `ROOT::Minuit2::MnMigrad`).  The cost function
         * is a weighted χ² sum over all unmasked pixels of the input stamp:
         * @f[
         *   \chi^2 = \sum_i w_i \bigl(d_i - m_i(\boldsymbol{p})\bigr)^2
         * @f]
         * where @f$d_i@f$ is the observed pixel value, @f$m_i(\boldsymbol{p})@f$ the model
         * value from `Eval()`, and @f$w_i@f$ is either 1 (no weight map) or @f$1/\sigma_i^2@f$
         * (error map) or the user-supplied weight.
         *
         * ### Supported PSF models
         * The profile model is selected at construction via `Stars::evalType`.  All models
         * add a constant background term as the last parameter.  The parameter vector
         * convention for each model is documented in `Eval()`.
         *
         * ### Extending with a custom profile
         * Derive from `starFitter` and override `user_eval(x, y, params)`.  That method
         * is called automatically by `Eval()` when `fevalMethod == Stars::evalType::UNDEFINED`.
         *
         * @note All fitting is done in **pixel** space.  WCS conversion (if needed) is the
         *       responsibility of the caller.
         *
         * @see Stars::evalType
         * @see DST::algorithms::starFinder::EvalStellarProperties()
         */
        class starFitter : public ROOT::Minuit2::FCNBase
        {
            private:
                std::shared_ptr<const FITSimg<double>> fsky;     ///< Image stamp to fit (always double precision).
                std::shared_ptr<const FITSimg<double>> fweight;  ///< Per-pixel weight map; null → uniform weight 1.
                const DST::Stars::evalType fevalMethod;          ///< PSF profile selected at construction.

                /**
                 * @brief Build a per-pixel weight map from an error or weight image.
                 * @param err        Error or weight `FITSimg<double>`.
                 * @param isErrorMap `true` → `weight = 1/err²`; `false` → use @p err as weights directly.
                 * @return           Const shared pointer to the weight image.
                 */
                static std::shared_ptr<const FITSimg<double>> _makeWeight(const std::shared_ptr<FITScube>& err, bool isErrorMap);
                /// @copydoc _makeWeight(const std::shared_ptr<FITScube>&, bool)
                static std::shared_ptr<const FITSimg<double>> _makeWeight(const std::shared_ptr<FITSimg<double>>& err, bool isErrorMap);

            protected:
                /**
                 * @brief Override point for user-defined PSF models.
                 *
                 * @details
                 * Called by `Eval()` when `fevalMethod == Stars::evalType::UNDEFINED`.
                 * The base implementation prints a warning and throws
                 * `std::invalid_argument`.  Override in a derived class to implement
                 * a custom profile without touching `operator()()`.
                 *
                 * @param x       Pixel x-coordinate (column, can be fractional).
                 * @param y       Pixel y-coordinate (row, can be fractional).
                 * @param params  Parameter vector passed through from `operator()()`.
                 * @return        Model value at (x, y).
                 * @throws std::invalid_argument in the base implementation.
                 */
                virtual double user_eval(double x, double y, const std::vector<double>& params) const;

                /**
                 * @brief Override point for the parameter layout of user-defined PSF models.
                 *
                 * @details
                 * Called by `labelParameters()` when `fevalMethod` designates a model whose
                 * layout is not known to the library (`MOMENT`, `CUSTOM`, `CUSTOM_MODEL`,
                 * `UNDEFINED`).  The base implementation returns an **empty**
                 * `MnUserParameters`; override it in a derived class to declare the
                 * parameters expected by `user_eval()` (names, start values, steps and
                 * optional limits), in the very same order.
                 *
                 * @return Parameter set matching the custom model; empty by default.
                 */
                virtual ROOT::Minuit2::MnUserParameters user_parameters() const;

            public:

                /**
                 * @brief Build the Minuit2 parameter set (names, seeds, steps, limits) of a PSF model.
                 *
                 * @details
                 * The returned layout matches exactly the parameter vector convention of
                 * `Eval()`: `psfFlux, psfX0, psfY0, <profile parameters...>, psfBKG`.  Values are
                 * generic seeds (`psfFlux = 1`, centroid at the origin); the caller is expected
                 * to re-seed them with the actual star measurements, e.g.
                 * @code
                 * auto upar = starFitter::generateLabelParameters(method);
                 * upar.SetValue("psfFlux", star.maxima());
                 * upar.SetValue("psfX0",   star.position().X() - xMin);
                 * upar.SetValue("psfY0",   star.position().Y() - yMin);
                 * @endcode
                 *
                 * @param method PSF profile model.
                 * @return Parameter set for @p method; an **empty** set for `MOMENT`,
                 *         `CUSTOM`, `CUSTOM_MODEL`, `UNDEFINED` and any unknown value —
                 *         those layouts are supplied by `user_parameters()`.
                 */
                static ROOT::Minuit2::MnUserParameters generateLabelParameters(const DST::Stars::evalType& method);

                /**
                 * @brief Parameter set for the model this fitter was constructed with.
                 *
                 * @details
                 * Dispatches to `generateLabelParameters(fevalMethod)` for the built-in
                 * profiles, and to the virtual `user_parameters()` for custom/unknown ones,
                 * so a derived class only has to override `user_eval()` and
                 * `user_parameters()` to be fully usable.
                 */
                ROOT::Minuit2::MnUserParameters labelParameters() const;

#pragma region -- ctor/dtor
                /**
                 * @brief Construct from a `FITScube` stamp with uniform pixel weights.
                 * @param sky     2-D image stamp; cast to double internally.
                 * @param method  PSF profile model to use.
                 * @throws FITSexception if @p sky is null or not 2-D.
                 */
                starFitter(const std::shared_ptr<FITScube>& sky, const DST::Stars::evalType& method);
                /**
                 * @brief Construct from a `FITSimg<double>` stamp with uniform pixel weights.
                 * @param sky     2-D double-precision image stamp.
                 * @param method  PSF profile model.
                 * @throws FITSexception if @p sky is null or not 2-D.
                 */
                starFitter(const std::shared_ptr<FITSimg<double>>& sky, const DST::Stars::evalType& method);
                /**
                 * @brief Construct from a `FITScube` stamp and a `FITScube` error/weight image.
                 * @param sky        Image stamp (cast to double).
                 * @param err        Error or weight image (same dimensions as @p sky).
                 * @param method     PSF profile model.
                 * @param isErrorMap `true` (default) → weight = 1/err²; `false` → use @p err as weights.
                 * @throws FITSexception if any pointer is null, or dimensions mismatch.
                 */
                starFitter(const std::shared_ptr<FITScube>& sky, const std::shared_ptr<FITScube>& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                /// @copydoc starFitter(const std::shared_ptr<FITScube>&, const std::shared_ptr<FITScube>&, const DST::Stars::evalType&, bool)
                starFitter(const std::shared_ptr<FITScube>& sky, const std::shared_ptr<FITSimg<double>>& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                /// @copydoc starFitter(const std::shared_ptr<FITScube>&, const std::shared_ptr<FITScube>&, const DST::Stars::evalType&, bool)
                starFitter(const std::shared_ptr<FITSimg<double>>& sky, const std::shared_ptr<FITScube>& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                /// @copydoc starFitter(const std::shared_ptr<FITScube>&, const std::shared_ptr<FITScube>&, const DST::Stars::evalType&, bool)
                starFitter(const std::shared_ptr<FITSimg<double>>& sky, const std::shared_ptr<FITSimg<double>>& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                ~starFitter();
#pragma endregion
#pragma region -- operator
                /**
                 * @brief Minuit2 cost function — returns the weighted χ² for parameter vector @p params.
                 *
                 * @details
                 * Iterates over all unmasked pixels of `fsky`, evaluates the PSF model
                 * via `Eval()`, and accumulates:
                 * @f$ \chi^2 = \sum_i w_i (d_i - m_i)^2 @f$.
                 * If `Eval()` throws (invalid parameters), returns `1e30` to steer Minuit2
                 * away from that region of parameter space.
                 *
                 * @param params  Parameter vector; layout depends on `fevalMethod` (see `Eval()`).
                 * @return        Weighted χ² value.
                 */
                double operator()(const std::vector<double>&) const override;
#pragma endregion

                /**
                 * @brief Minuit2 error definition — returns 1.0 (χ² formalism).
                 *
                 * @details
                 * Returning 1.0 means Minuit2 defines the 1-σ parameter uncertainty as the
                 * change in parameters that increases χ² by 1.  Override in a derived class
                 * to switch to likelihood (return 0.5) or a custom scale.
                 */
                double Up() const override { return 1.0; }

                /// @brief Evaluate the PSF model at a single pixel (x, y).
                /// @param x       Pixel x-coordinate (column, can be fractional).
                /// @param y       Pixel y-coordinate (row, can be fractional).
                /// @param params  Parameter vector (same convention as operator()).
                /// @return Model value (signal + background) at (x, y).
                /// @note For UNDEFINED evalType calls user_eval(); override that in a derived class.
                /// @note Under C++20 this takes a view, so a `std::vector<double>` binds to it
                ///       directly and a caller holding a longer vector can pass a sub-range
                ///       without copying. Under C++17 it takes the vector itself and the
                ///       sub-range case costs a copy. See DST_HAS_SPAN.
#if DST_HAS_SPAN
                double Eval(double x, double y, std::span<const double> params) const;
#else
                double Eval(double x, double y, const std::vector<double>& params) const;
#endif

                /**
                 * @brief Evaluate the PSF model from a single coordinates-first vector.
                 *
                 * @details
                 * Convenience overload for callers that carry the evaluation point inside the
                 * parameter vector rather than beside it:
                 * @code
                 * // { x, y, psfFlux, psfX0, psfY0, <profile...>, psfBKG }
                 * double v = fitter.Eval(params);
                 * @endcode
                 * Equivalent to `Eval(params[0], params[1], <the rest>)`. Under C++20 the tail
                 * is passed as a view and nothing is copied; under C++17 the tail is copied
                 * into a temporary vector, which is the one place the fallback costs anything.
                 *
                 * @param params Evaluation point followed by the model parameters.
                 * @return Model value (signal + background) at (params[0], params[1]).
                 * @throws std::invalid_argument if @p params holds fewer than three entries.
                 */
                double Eval(const std::vector<double>& params) const;

#pragma region -- ploting
#ifdef __WITH_QT__
                /// @brief Display fsky heatmap + centroid cross + model contours on axisRect.
                /// @param axisRect   Target QCPAxisRect (must belong to a QCustomPlot).
                /// @param params     Fitted parameter vector (params[1]=x0, params[2]=y0).
                /// @param oversample Per-dimension model oversampling factor (total pixels =
                ///                   oversample² × fsky pixels). Clamped to [1,100]. Default 10.
                /// @param saveFile   If non-empty, save the plot to file
                ///                   (.png by default; .jpg and .pdf are also supported).
                void Plot(QCPAxisRect* axisRect,
                          const std::vector<double>& params,
                          size_t oversample = 10,
                          const std::string& saveFile = "") const;

                /// @brief Convenience overload accepting ROOT::Minuit2::MnUserParameters.
                void Plot(QCPAxisRect* axisRect,
                          const ROOT::Minuit2::MnUserParameters& params,
                          size_t oversample = 10,
                          const std::string& saveFile = "") const;
#endif
#pragma endregion
        };
#pragma endregion

#pragma region - DST::starFinder class definition

        /**
         * @brief Multi-method star detection engine for 2-D and 3-D FITS images.
         *
         * @details
         * `starFinder` wraps five detection pipelines selectable at runtime via
         * `setDetectionMethod()` or the two-argument `findStars()` overload:
         *
         * | Method | Description |
         * |--------|-------------|
         * | `MORPHO` (default) | Morphological opening/erosion/dilation to isolate stellar blobs, followed by watershed flooding and moment analysis. |
         * | `THRESHOLD` | Pixels above a user-defined threshold are grown into blobs. |
         * | `SNR` | S/N map thresholding — blobs whose peak S/N exceeds the configured value survive. |
         * | `NOISECHISEL` | Gnuastro-inspired tile-by-tile quantile thresholding with sky subtraction and re-detection. |
         * | `CUSTOM` | Hook for user-defined pipelines via `_customPreprocessing()`. |
         *
         * After detection, moment-based stellar properties (centroid, FWHM, flux, SNR) are
         * computed automatically.  Refined PSF fits can be requested via
         * `EvalStellarProperties()`.  The result can be exported as a `FITStable` via
         * `detectedStars()`, or as a SCAMP-compatible FITS-LDAC file via `toFITS_LDAC()`.
         *
         * 3-D data cubes are processed layer by layer with up to 4 concurrent threads.
         *
         * @see Stars
         * @see starFitter
         */
        class starFinder
        {
            public:
                /**
                 * @brief Detection algorithm to run on the input image.
                 * @see findStars()
                 */
                enum class detectionMethod
                {
                    MORPHO,     ///< Morphological pipeline (opening + watershed): robust default.
                    THRESHOLD,  ///< Simple pixel-value threshold.
                    SNR,        ///< Signal-to-noise ratio map threshold.
                    NOISECHISEL,///< Tile-based quantile thresholding with sky subtraction.
                    CUSTOM      ///< User-defined pipeline via `_customPreprocessing()`.
                };

                typedef std::map<const std::string, DST::Morphology::kernel<double> > kernelMap; ///< Type alias for a map of named morphology kernels.
                typedef std::map<const std::string, DST::Morphology::kernel<double> >::iterator kIterator;
                typedef std::map<const std::string, DST::Morphology::kernel<double> >::const_iterator const_kIterator;

            protected:
                std::map<size_t, Stars> fdetectedStars;         ///< All detected stars, keyed by ID.

            private:
                const std::shared_ptr<FITScube> fsky;          ///< Input image (never null after construction).
                std::shared_ptr<FITScube> flabel;               ///< Integer label image produced during detection.
                detectionMethod fmethod = detectionMethod::MORPHO; ///< Active detection pipeline.
                double fthresholdValue;                         ///< Pixel threshold for `THRESHOLD` mode (default: 97th percentile).

                kernelMap fKernels;

                std::unique_ptr<DST::GPU::Backend> fbackend_; ///< GPU backend; always set, falls back to CPU if no GPU is present.

                /**
                 * @brief Tunable parameters for the NOISECHISEL detection pipeline.
                 *
                 * @details
                 * Control every stage of the tile-based sky estimation and blob-detection
                 * pipeline inspired by Gnuastro NoiseChisel.  All fields have conservative
                 * defaults suitable for typical astronomical images.
                 */
                struct NoiseChiselParams
                {
                    double convSigma         = 1.5;  ///< Gaussian pre-smoothing kernel σ [pixels].
                    double detectionQuantile = 0.70; ///< Per-tile quantile that defines the detection threshold (0–1).
                    size_t tileSize          = 50;   ///< Tessellation tile side length [pixels].
                    double snrThreshold      = 3.0;  ///< Minimum blob S/N to survive the S/N filter.
                    size_t skyIterations     = 3;    ///< Number of σ-clipping iterations for sky estimation.
                    double skyClipSigma      = 3.0;  ///< Rejection threshold for sky σ-clipping.
                    bool   reDetect          = true; ///< If `true`, re-run detection on the sky-subtracted image.
                } fNC;

                void _init_NCgaussKernel(const size_t &);
                void _init_NCerodeKernel(const size_t &);
                void _init_NCdilateKernel(const size_t &);

                /// @brief Run the morphological detection pipeline on a single 2-D layer.
                void _findStarsMorpho(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);

                /// @brief Run the threshold-based detection pipeline on a single 2-D layer.
                void _findStarsThreshold(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);

                /// @brief Run the S/N-based detection pipeline on a single 2-D layer.
                void _findStarsSNR(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);

                /// @brief Run the NoiseChisel detection pipeline on a single 2-D layer.
                void _findStarsNoiseChisel(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);

                /// @brief Run the custom detection pipeline on a single 2-D layer.
                void _findStarsCustom(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);

                /// @brief Grow labelled blobs via priority-queue watershed flooding.
                void _watershedFlooding(const std::shared_ptr<FITScube>&, const double&, const double&, const size_t& iLayer=0);

                /// @brief Compute image-moment stellar properties (centroid, flux, FWHM).
                void _evalStellarPropertiesMoment(Stars& star);

                /// @brief Fit a parametric PSF model to the star stamp via Minuit2.
                void _evalStellarPropertiesModel(Stars& star, const DST::Stars::evalType& how);

                /// @brief Return the pixel indices of the 4- or 8-connected neighbours of @p pixIdx.
                std::vector<size_t> _getNeighborIndices(const size_t& pixIdx, const size_t& w, const size_t& h) const;

                /// @brief Flood-fill a flat-topped plateau in the priority-queue watershed.
                std::vector<point> _floodingPlateau(const size_t&, std::valarray<size_t>&, const std::shared_ptr<FITScube>, const double&, const size_t&);

            protected:
                std::string foutURL; ///< Optional output directory for intermediate debug products.
                bool is3D;           ///< `true` when the input image has 3 FITS axes.

                /// @brief Compute stellar properties using user-defined method.
                virtual void _evalStellarPropertiesCustom(Stars& star) {_evalStellarPropertiesMoment(star);} ///< Default: fall back to moment analysis; override for custom behaviour.

                /// @brief Fit a parametric PSF model to the star stamp via Minuit2.
                virtual void _evalStellarPropertiesCustomModel(Stars& star, const DST::Stars::evalType& how) {_evalStellarPropertiesModel(star, how);} ///< Default: fall back to built-in models; override for custom behaviour.

                /// @brief Compute aperture photometry for a star using the label mask.
                virtual void _evalApperturePhotometry(Stars& star);

                /**
                 * @brief Override point for user-defined image pre-processing.
                 *
                 * @details
                 * Called by `_findStarsCustom()` before blob detection.  The base
                 * implementation returns the input unchanged.  Override in a derived class
                 * to apply deconvolution, cosmic-ray rejection, or any other custom
                 * processing before star detection.
                 *
                 * @param img    Input 2-D layer.
                 * @param iLayer Layer index (for 3-D cubes).
                 * @return       Pre-processed image to use for detection.
                 */
                virtual std::shared_ptr<FITScube> _customPreprocessing(const std::shared_ptr<FITScube>&, const size_t& iLayer=0);

                std::shared_ptr<FITScube> _morphologicalPreprocessing(const std::shared_ptr<FITScube>&,const size_t& iLayer=0) const;
                std::shared_ptr<FITScube> _SNRPreprocessing(const std::shared_ptr<FITScube>&,const size_t& iLayer=0) const;
                std::shared_ptr<FITScube> _ThresholdPreprocessing(const std::shared_ptr<FITScube>&,const size_t& iLayer=0) const;
                std::shared_ptr<FITScube> _NoiseChiselPreprocessing(const std::shared_ptr<FITScube>&,const size_t& iLayer=0) const;

            public:
#pragma region -- initilizers
                void defaultMorphologicalKernels(); ///< @brief Set the default morphological kernels.
                void defaultThresholdKernels();
                void defaultSNRKernels();
                void defaultNoiseChiselKernels();
#pragma endregion

#pragma region -- constructors / destructor
                /**
                 * @brief Construct a `starFinder` bound to a sky image with a specific detection method.
                 *
                 * @details
                 * Initialises the morphological kernels, selects the GPU backend via
                 * `DST::GPU::Backend::create()` (falls back to CPU if no GPU is present),
                 * and sets the default threshold to the 97th image percentile.
                 *
                 * @param skyImage  Input 2-D or 3-D FITS image (shared ownership).
                 * @param method    Detection algorithm to use (default overload uses `MORPHO`).
                 * @throws FITSexception if @p skyImage is null or its dimension is not 2 or 3.
                 */
                starFinder(const std::shared_ptr<FITScube>& skyImage, const detectionMethod& method);
                /**
                 * @brief Construct a `starFinder` using the default `MORPHO` detection method.
                 * @param skyImage  Input 2-D or 3-D FITS image (shared ownership).
                 * @throws FITSexception if @p skyImage is null or its dimension is not 2 or 3.
                 */
                starFinder(const std::shared_ptr<FITScube>& skyImage);
                ~starFinder();
#pragma endregion
#pragma region -- modifiers

#pragma region kernel modifiers
                
                inline const bool kernelExists(const std::string& name) const {return fKernels.find(name) != fKernels.end();} ///< @brief Check if a kernel with the given name exists.

                inline kIterator findKernel(const std::string& name) {return fKernels.find(name);} ///< @brief Retrieve a kernel by name;
                inline const_kIterator findKernel(const std::string& name) const {return fKernels.find(name);} ///< @brief Retrieve a kernel by name (const version).

                inline kIterator Kbegin() {return fKernels.begin();} ///< @brief Get an iterator to the first kernel in the map.
                inline kIterator Kend() {return fKernels.end();} ///< @brief Get an iterator to one past the last kernel in the map.
                
                inline const_kIterator Kcbegin() const {return fKernels.cbegin();} ///< @brief Get a const iterator to the first kernel in the map.
                inline const_kIterator Kcend() const {return fKernels.cend();} ///< @brief Get a const iterator to one past the last kernel in the map.

                /// @brief Add a custom morphological kernel with a given name.
                void insertKernel(const std::string& name, const DST::Morphology::kernel<double>& k);

                void updateKernel(const std::string& name, const DST::Morphology::kernel<double>& k); ///< @brief Update an existing kernel by name; throws if the kernel does not exist.
                
                /// @brief Remove a custom kernel by name.
                void removeKernel(const std::string& name);
                
#pragma endregion


                /// @brief Select the detection algorithm for the next `findStars()` call.
                inline void setDetectionMethod(const detectionMethod& m) {fmethod = m;}

                /// @brief Remove all previously detected stars from the internal map.
                inline void clearDetectedStars() {fdetectedStars.clear();}

                /// @brief Set the pixel-value threshold used by the `THRESHOLD` pipeline.
                inline void setThresholdValue(const double& v) {fthresholdValue = v;}

                /**
                 * @brief Set the output directory for debug / intermediate products.
                 *
                 * @details
                 * When set, the label image and other intermediate results are written to
                 * this directory after each `findStars()` call.  The directory is created
                 * if it does not already exist.
                 *
                 * @param url  Absolute or relative path to the output directory (length > 3).
                 * @throws std::invalid_argument if @p url is shorter than 4 characters.
                 */
                void setOutputUrl(const std::string& url);

                /// @brief Change the operation type of a kernel.
                inline void switchkernelType(const std::string& name, const DST::Morphology::ktype& t)
                {
                    kIterator it = findKernel(name);
                    if(it == Kend())
                        throw std::invalid_argument("Kernel with name '" + name + "' does not exist.");

                    it->second.setType(t);
                }

                /// @name NoiseChisel parameter setters
                /// @{
                void setNC_convSigma(double v)         { fNC.convSigma         = v; } ///< Set Gaussian pre-smoothing σ [pixels].
                void setNC_convRadius(size_t v)        { _init_NCgaussKernel(v);} ///< Set Gaussian kernel half-size [pixels].
                void setNC_detectionQuantile(double v) { fNC.detectionQuantile = v; } ///< Set per-tile quantile threshold (0-1).
                void setNC_tileSize(size_t v)          { fNC.tileSize          = v; } ///< Set tessellation tile side [pixels].
                void setNC_erosionRadius(size_t v)     { _init_NCerodeKernel(v); } ///< Set erosion element half-size [pixels].
                void setNC_snrThreshold(double v)      { fNC.snrThreshold      = v; } ///< Set minimum blob S/N to survive the S/N filter.
                void setNC_dilationRadius(size_t v)    { _init_NCdilateKernel(v); } ///< Set dilation half-size for mask growth [pixels].
                void setNC_skyIterations(size_t v)     { fNC.skyIterations     = v; } ///< Set number of sky σ-clipping iterations.
                void setNC_skyClipSigma(double v)      { fNC.skyClipSigma      = v; } ///< Set sky σ-clipping rejection threshold.
                void setNC_reDetect(bool v)            { fNC.reDetect          = v; } ///< Enable/disable re-detection on sky-subtracted image.
                /// @}

#pragma endregion

#pragma region -- accessors
                /**
                 * @brief Build and return a `FITStable` with the given name containing all detected stars.
                 *
                 * @details
                 * Columns: `ID`, `coord_px`, `maxcoord_px`, optional WCS columns (read from
                 * the image header), `fits` (PSF parameters), `maxima`, `FWHM_X`, `FWHM_Y`,
                 * `angle`, `TotalFlux`, `SNR`, `Background`, `EvalMethod`, `n_pixels`, `source`.
                 *
                 * @param tblName  Extension name for the binary table HDU.
                 * @return         `FITStable` ready for writing with `FITSmanager`.
                 */
                FITStable detectedStars(const std::string& tblName) const;
                /// @brief Build a `FITStable` using the default extension name `"DETECTED.STARS"`.
                /// @return `FITStable` containing all detected stars.
                FITStable detectedStars() const;

                /// @return the threshold value used by the `THRESHOLD` detection pipeline.
                inline double thresholdValue() const { return fthresholdValue; }

                /**
                 * @brief Export the star catalogue to a FITS-LDAC file for use with SCAMP.
                 *
                 * @details
                 * Produces a two-extension FITS file:
                 * - `LDAC_IMHEAD` — one row containing the full source image header.
                 * - `LDAC_OBJECTS` — one row per star with SExtractor-compatible columns
                 *   (`X_IMAGE`, `Y_IMAGE`, `FLUX_AUTO`, `MAG_AUTO`, `FWHM_IMAGE`, …).
                 *
                 * If @p method is not `UNDEFINED`, stars whose `evalMethod()` differs are
                 * re-evaluated in parallel before writing.
                 *
                 * @param filename  Output FITS file path (created or overwritten).
                 * @param method    PSF model to use for re-evaluation; `UNDEFINED` skips it.
                 */
                void toFITS_LDAC(const std::string& filename,
                                 Stars::evalType method = Stars::evalType::UNDEFINED);

                /// @brief Return the active detection method.
                inline const detectionMethod& getDetectionMethod() const {return fmethod;}
                /// @brief Return the active detection method as a string (e.g. `"MORPHO"`).
                std::string detectionMethodString() const;

                /// @brief Return the opening (BS) kernel (const reference).
                const DST::Morphology::kernel<double>& getKernel(const std::string& name) const
                {
                    const_kIterator it = findKernel(name);
                    if (it == Kcend())
                    {
                        throw std::invalid_argument("Kernel with name '" + name + "' does not exist.");
                    }
                    return it->second;
                }

                /// @brief Return `true` when a real GPU backend (Vulkan/CUDA/HIP) is active.
                /// When `false`, morphological operations run on CPU regardless of image type.
                bool isGPUEnabled() const { return fbackend_ && fbackend_->available(); }
#pragma endregion
#pragma region -- processing
                /**
                 * @brief Detect stars in the bound image using the currently selected method.
                 *
                 * @details
                 * For 3-D cubes each layer is processed in a separate `std::async` task
                 * (up to 4 concurrent threads).  After detection, moment-based stellar
                 * properties are computed automatically for all blobs.
                 *
                 * @param quantizationStep  Intensity quantization step used by the watershed
                 *                          flooding queue (default: 1000).  Smaller values give
                 *                          finer separation of close objects at the cost of speed.
                 * @throws FITSexception if the input image pointer is null or has wrong dimensions.
                 */
                void findStars(double quantizationStep=1e3);
                /**
                 * @brief Detect stars using @p method, then delegate to `findStars(quantizationStep)`.
                 * @param method           Detection pipeline to activate.
                 * @param quantizationStep Watershed intensity step (default: 1000).
                 */
                void findStars(const enum detectionMethod& method, double quantizationStep=1e3);

                inline void EvalStellarProperties(const enum Stars::evalType& method)
                {
                    for (auto& [id, star] : fdetectedStars)
                    {
                        EvalStellarProperties(star, method);
                    }
                }

                /**
                 * @brief Evaluate stellar properties for the star identified by @p id.
                 * @param id      Identifier of the star in `fdetectedStars`.
                 * @param method  PSF profile model to fit.
                 * @throws std::out_of_range if @p id is not found.
                 */
                void EvalStellarProperties(const size_t&, const enum Stars::evalType& method);
                /**
                 * @brief Evaluate stellar properties for a `Stars` object in place.
                 * @param star    Reference to the `Stars` object to update.
                 * @param method  PSF profile model to fit.
                 */
                void EvalStellarProperties(Stars&, const enum Stars::evalType& method);
#pragma endregion
#pragma region -- ploting

#ifdef __WITH_QT__
                /**
                 * @brief Overlay detected star markers on a Qt axis rect.
                 * @param axisRect  Target `QCPAxisRect` (must belong to a live `QCustomPlot`).
                 * @param iL        Layer index to display (0 for 2-D images).
                 * @param starPen   Pen used to draw the star marker ellipses (default: red, 2 px).
                 */
                void Plot(QCPAxisRect* axisRect, size_t iL=0, const QPen& starPen=QPen(Qt::red, 2)) const;
#endif
#pragma endregion
        };
#pragma endregion
    }
}

#endif