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
#include <DSTmath/DSM_version.h>
#include <DSTmath/DSTgeo.h>

#include "DSTmorpho.h"


using namespace DSL;
using namespace DST::Math;

namespace DST
{

#pragma region - DST::starFinder class definition
    class Stars
        {
            public:
                enum class evalType
                {
                    MOMENT,
                    GAUSS,
                    GAUSS2D,
                    MOFFAT,
                    MOFFAT2D
                };

            private:
                const size_t fID;
                evalType fevalMethod;       //< Method used to evaluate star parameters
                std::vector<point> fpixels; //< List of pixel detected to belong to the star
                point fposition;            //< Centroid position of the star   

                double fMaxima;             //< Maximum intensity of the star
                point fMaximaPosition;     //< Position of the maximum intensity of the star
                point fFWHM;                //< FWHM of the star (x and y)    
                double fangle;              //< Angle of the star (orientation)
                double fTotalFlux;          //< Total flux of the star
                double fSNR;                //< Signal to noise ratio of the star
                double fBackground;         //< Background level around the star

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
#pragma endregion
#pragma region -- Accessors
                inline const size_t& ID()                       const {return fID;}
                inline const std::vector<point>& FilledPixels() const {return fpixels;}
                inline const point& position()                  const {return fposition;}
                inline const double& maxima()                   const {return fMaxima;}
                inline const point& maximaPosition()            const {return fMaximaPosition;}
                inline const point&                             FWHM() const {return fFWHM;}
                inline const double&                            angle() const {return fangle;}
                inline const double&                            totalFlux() const {return fTotalFlux;}
                inline const double&                            SNR() const {return fSNR;}
                inline const double&                            background() const {return fBackground;}
        };
#pragma endregion

    namespace algorithms
    {

        /**
         * @class starFinder DSTools/DSTstarFinder.h
         * @details This class implement algorithm to find stars on an images. Algorithm are based on morphological operators to isolate and locate stars. Parameters can be tuned to adapt the detection to different kind of images (pixel scale, PSF size, background noise, etc...).
         * 
         */
        class starFinder
        {
            public:

                enum class detectionMethod
                {
                    MORPHO,
                    THRESHOLD,
                    SNR
                };

            private:
                const std::shared_ptr<FITScube> fsky; //< Pointer to the input image.
                detectionMethod fmethod = detectionMethod::MORPHO; //< Detection method to use.
                std::map<size_t, Stars> fdetectedStars; //< Map of detected stars with their ID and position.
                double fthresholdValue; //< Threshold value for THRESHOLD detection method.

                //-- Kernel for the morphological detection
                DST::Morphology::kernel<double> fBS;
                DST::Morphology::kernel<double> fBE;
                DST::Morphology::kernel<double> fBM;

                std::string foutURL;
                bool is3D;

                void _findStarsMorpho(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _findStarsThreshold(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _findStarsSNR(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _watershedFlooding(const std::shared_ptr<FITScube>&, const double&, const size_t& iLayer=0);
                void _evalStellarPropertiesMoment(Stars& star);
                std::vector<size_t> _getNeighborIndices(const size_t& pixIdx, const size_t& w, const size_t& h) const;
                std::vector<point>  _floodingPlateau(const size_t&, std::valarray<size_t>&, const std::shared_ptr<FITScube>, const double&, const size_t&);

            public:

#pragma region • constructors / destructor
                starFinder(const std::shared_ptr<FITScube>& skyImage);
                ~starFinder();
#pragma endregion
#pragma region • modifiers
                inline void setDetectionMethod(const detectionMethod& m) {fmethod = m;}
                inline void clearDetectedStars() {fdetectedStars.clear();}
                inline void setThresholdValue(const double& v) {fthresholdValue = v;}
                void setOutputUrl(const std::string& url);
#pragma endregion
#pragma region • accessors
                inline const std::map<size_t, Stars>& detectedStars() const {return fdetectedStars;}
                inline const detectionMethod& detectionMethod() const {return fmethod;}
                std::string detectionMethodString() const;
#pragma endregion
#pragma region • processing
                void findStars(double quantizationStep=1e3);
                void findStars(const enum detectionMethod& method, double quantizationStep=1e3);
                void measureStars();
                void measureStar(const size_t& starID);
                void EvalStellarProperties(const size_t&, const enum Stars::evalType& method);
                void EvalStellarProperties(Stars&, const enum Stars::evalType& method);
#pragma endregion
        };
#pragma endregion
    }
}

#endif