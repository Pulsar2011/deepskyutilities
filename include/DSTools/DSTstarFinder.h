//
//  DSTstarFinder.h
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

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
    class Stars
        {
            public:
                enum class evalType
                {
                    MOMENT,
                    GAUSS,      //< SYMMETRIC GAUSSIAN FIT
                    GAUSS2D,    //< ASYMMETRIC GAUSSIAN FIT
                    iGAUSS,     //< INTEGRATED SYMMETRIC GAUSSIAN FIT
                    iGAUSS2D,   //< INTEGRATED ASYMMETRIC GAUSSIAN FIT
                    MOFFAT,      //< SYMMETRIC MOFFAT FIT
                    MOFFAT2D,   //< ASYMMETRIC MOFFAT FIT
                    LORENTZIAN,  //< SYMMETRIC LORENTZIAN FIT
                    LORENTZIAN2D,//< ASYMMETRIC LORENTZIAN FIT
                    VOIGT,       //< SYMMETRIC PSEUDO-VOIGT FIT
                    VOIGT2D,     //< ASYMMETRIC PSEUDO-VOIGT FIT
                    SERSIC,      //< SYMMETRIC SERSIC FIT
                    SERSIC2D,    //< ASYMMETRIC SERSIC FIT
                    KING,        //< SYMMETRIC KING FIT
                    KING2D,      //< ASYMMETRIC KING FIT
                    AIRY,        //< SYMMETRIC AIRY FIT
                    UNDEFINED
                };

            private:
                const size_t fID;
                evalType fevalMethod;       //< Method used to evaluate star parameters
                std::vector<point> fpixels; //< List of pixel detected to belong to the star
                point fposition;            //< Centroid position of the star [pixels]
                std::vector<double> fPSFparameters; //< Parameters of the PSF fitted to the star. Depends on the fEvalMethods.

                double fMaxima;             //< Maximum intensity of the star                 [UNIT of image pixel values]
                point fMaximaPosition;      //< Position of the maximum intensity of the star [pixels]
                point fFWHM;                //< FWHM of the star (x and y)                    [pixels]
                double fangle;              //< Angle of the star (orientation)               [red]
                double fTotalFlux;          //< Total flux of the star                        [UNTI of image pixel values]
                double fSNR;                //< Signal to noise ratio of the star
                double fBackground;         //< Background level around the star              [UNIT of image pixel values]
                std::string finURL;

            public:
                Stars();
                Stars(const size_t& id);
                Stars(const size_t& id, const pixelCoords& pos);
                Stars(const size_t& id, const point& pos);
                ~Stars();

#pragma region -- Modifiers
                void addPixel(const pixelCoords& );
                void addPixel(const pixelVectors& );
                void addPixel(const std::vector<point>& );
                void addPixel(const point&);

                void setPosition(const pixelCoords& );
                void setPosition(const point& );

                void setMaxima(const double& );
                void setMaximaPosition(const pixelCoords& );
                void setMaximaPosition(const point& );
                void setFWHM(const double& fwhmX, const double& fwhmY);
                void setFWHM(const point& fwhm);
                void setAngle(const double& angle);
                void setTotalFlux(const double& flux);
                void setSNR(const double& snr);
                void setBackground(const double& background);
                void setEvalMethod(const evalType& method);
                void setPSFparameters(const std::vector<double>& params);
#pragma endregion
#pragma region -- Accessors
                inline const size_t&                            ID() const {return fID;}
                inline const std::vector<point>&                FilledPixels() const {return fpixels;}
                inline const point&                             position() const {return fposition;}
                inline const double&                            maxima() const {return fMaxima;}
                inline const point&                             maximaPosition() const {return fMaximaPosition;}
                inline const point&                             FWHM() const {return fFWHM;}
                inline const double&                            angle() const {return fangle;}
                inline const double&                            totalFlux() const {return fTotalFlux;}
                inline const double&                            SNR() const {return fSNR;}
                inline const double&                            background() const {return fBackground;}
                inline const evalType&                          evalMethod() const {return fevalMethod;}
                inline const std::vector<double>&               PSFparameters() const {return fPSFparameters;}

                const std::string evalMethodString() const;
                static const std::string ModelString(const evalType& m);
                inline const std::string&                       srcFile() const {return finURL;}
#pragma endregion
#pragma region -- Dumper
                void dump(std::ostream& os = std::cout, bool _header_=true) const;
#pragma endregion
        };
#pragma endregion

    namespace algorithms
    {

#pragma region - DST::starFitter class definition
        // Helper classe to perform Star fitting profile.
        class starFitter : public ROOT::Minuit2::FCNBase 
        {
            /**
            * @class starFitter DSTools/DSTstarFinder.h
            * @details Helper class used to fit stellar profile and evaluate stars properties. The fit is performed in the pixel space. If the input image has WCS defined, the conversion can be achieved afterward using WCS information and DST::FITScube conversion tools. the \c fevalMethod member is used to select the fitting function to be use. It is the responsibility of the user to define correct parameter vector according to the methods. The fitting is performed using the ROOT::Minuit2 library. To use PSF profile not difined in the current implementation, the user can inherit from this class and define its own \c user_operator() method to perform the fit. This function will be automatically called by the \c operator() if the \c fevalMethod is undefined. The \c Up() method is used to define the step size for the fit parameters. By default, it is set to 1.0, but it can be overridden in a derived class if needed.
            * @author GILLARD William
            */
            private:
                std::shared_ptr<const FITSimg<double> > fsky;
                std::shared_ptr<const FITSimg<double> > fweight;
                const DST::Stars::evalType fevalMethod;

            static std::shared_ptr< const FITSimg<double> > _makeWeight(const std::shared_ptr<FITScube>& err, bool isErrorMap);
            static std::shared_ptr< const FITSimg<double> > _makeWeight(const std::shared_ptr< FITSimg<double> >& err, bool isErrorMap);

            protected:
                virtual double user_eval(double x, double y, const std::vector<double>& params) const;

            public:
#pragma region -- ctor/dtor
                starFitter(const std::shared_ptr<FITScube>& sky, const DST::Stars::evalType& method);
                starFitter(const std::shared_ptr< FITSimg<double> >& sky, const DST::Stars::evalType& method);
                starFitter(const std::shared_ptr<FITScube>& sky, const std::shared_ptr<FITScube>& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                starFitter(const std::shared_ptr<FITScube>& sky, const std::shared_ptr< FITSimg<double> >& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                starFitter(const std::shared_ptr< FITSimg<double> >& sky, const std::shared_ptr<FITScube>& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                starFitter(const std::shared_ptr< FITSimg<double> >& sky, const std::shared_ptr< FITSimg<double> >& err, const DST::Stars::evalType& method, bool isErrorMap = true);
                ~starFitter();
#pragma endregion
#pragma region -- operator
                double operator()(const std::vector<double>&) const override;
#pragma endregion
                
                double Up() const override { return 1.0; }

                /// @brief Evaluate the PSF model at a single pixel (x, y).
                /// @param x       Pixel x-coordinate (column, can be fractional).
                /// @param y       Pixel y-coordinate (row, can be fractional).
                /// @param params  Parameter vector (same convention as operator()).
                /// @return Model value (signal + background) at (x, y).
                /// @note For UNDEFINED evalType calls user_eval(); override that in a derived class.
                double Eval(double x, double y, const std::vector<double>& params) const;

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
        
        class starFinder
        {
            /**
            * @class starFinder DSTools/DSTstarFinder.h
            * @details This class implement algorithm to find stars on an images. Algorithm are based on morphological operators to isolate and locate stars. Parameters can be tuned to adapt the detection to different kind of images (pixel scale, PSF size, background noise, etc...).
            * @author GILLARD William
            */
            public:

                enum class detectionMethod
                {
                    MORPHO,
                    THRESHOLD,
                    SNR,
                    CUSTOM
                };

            private:
                const std::shared_ptr<FITScube> fsky; //< Pointer to the input image.
                std::shared_ptr<FITScube> flabel; //< Pointer to the label image.
                detectionMethod fmethod = detectionMethod::MORPHO; //< Detection method to use.
                std::map<size_t, Stars> fdetectedStars; //< Map of detected stars with their ID and position.
                double fthresholdValue; //< Threshold value for THRESHOLD detection method.

                //-- Kernel for the morphological detection
                DST::Morphology::kernel<double> fBS;
                DST::Morphology::kernel<double> fBE;
                DST::Morphology::kernel<double> fBM;

                //-- GPU backend (always created; falls back to CPUBackend when no GPU hardware)
                std::unique_ptr<DST::GPU::Backend> fbackend_;

                void _findStarsMorpho(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _findStarsThreshold(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _findStarsSNR(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _findStarsCustom(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _watershedFlooding(const std::shared_ptr<FITScube>&, const double&, const double&, const size_t& iLayer=0);
                void _evalStellarPropertiesMoment  (Stars& star);
                void _evalStellarPropertiesModel   (Stars& star, const DST::Stars::evalType& how);
                void _evalApperturePhotometry      (Stars& star);
                std::vector<size_t> _getNeighborIndices(const size_t& pixIdx, const size_t& w, const size_t& h) const;
                std::vector<point>  _floodingPlateau(const size_t&, std::valarray<size_t>&, const std::shared_ptr<FITScube>, const double&, const size_t&);

            protected:
                std::string foutURL;
                bool is3D;

                virtual std::shared_ptr<FITScube> _customPreprocessing(const std::shared_ptr<FITScube>&,const size_t& iLayer=0);

            public:

#pragma region -- constructors / destructor
                starFinder(const std::shared_ptr<FITScube>& skyImage);
                ~starFinder();
#pragma endregion
#pragma region -- modifiers
                inline void setDetectionMethod(const detectionMethod& m) {fmethod = m;}
                inline void clearDetectedStars() {fdetectedStars.clear();}
                inline void setThresholdValue(const double& v) {fthresholdValue = v;}
                void setOutputUrl(const std::string& url);

                void switchfBSkernel(DST::Morphology::kernel<double>& k) {fBS = k;}
                void switchfBSkernelType(const DST::Morphology::ktype& t) {fBS.setType(t);}
                void switchfBEkernel(DST::Morphology::kernel<double>& k) {fBE = k;}
                void switchfBEkernelType(const DST::Morphology::ktype& t) {fBE.setType(t);}
                void switchfBMkernel(DST::Morphology::kernel<double>& k) {fBM = k;}
                void switchfBMkernelType(const DST::Morphology::ktype& t) {fBM.setType(t);}
#pragma endregion
#pragma region -- accessors
                FITStable detectedStars(const std::string& tblName) const;
                FITStable detectedStars() const;
                //inline const std::map<size_t, Stars>& detectedStars() const {return fdetectedStars;}
                inline const detectionMethod& detectionMethod() const {return fmethod;}
                std::string detectionMethodString() const;

                const DST::Morphology::kernel<double>& getfBSkernel() const {return fBS;}
                const DST::Morphology::kernel<double>& getfBEkernel() const {return fBE;}
                const DST::Morphology::kernel<double>& getfBMkernel() const {return fBM;}

                /// Returns true when a real GPU backend (Vulkan/CUDA/HIP) is active.
                /// When false, morphological operations run on CPU regardless of image type.
                bool isGPUEnabled() const { return fbackend_ && fbackend_->available(); }
#pragma endregion
#pragma region -- processing
                void findStars(double quantizationStep=1e3);
                void findStars(const enum detectionMethod& method, double quantizationStep=1e3);
                void measureStars();
                void measureStar(const size_t& starID);
                void EvalStellarProperties(const size_t&, const enum Stars::evalType& method);
                void EvalStellarProperties(Stars&, const enum Stars::evalType& method);
#pragma endregion
#pragma region -- ploting

#ifdef __WITH_QT__
                void Plot(QCPAxisRect* axisRect, size_t iL=0, const QPen& starPen=QPen(Qt::red, 2)) const;
#endif
#pragma endregion
        };
#pragma endregion
    }
}

#endif