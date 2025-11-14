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
#include <memory>
#include <valarray>
#include <cmath>
#include <algorithm>
#include <thread>
#include <future>
#include <limits>
#include <type_traits>

#include <DSTfits/FITSdata.h>
#include <DSTfits/FITSimg.h>
#include <DSTmath/DSM_version.h>
#include <DSTmath/DSTgeo.h>

#include "DST_version.h"
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
                size_t fID;
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
                void addPixel(const point&);
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
                    MAXIMA,
                    MORPHO,
                    THRESHOLD,
                    SNR
                };

            private:
                const std::shared_ptr<FITScube> fsky; //< Pointer to the input image.
                detectionMethod fmethod = detectionMethod::MORPHO; //< Detection method to use.
                std::map<size_t, DST::point> fdetectedStars; //< Map of detected stars with their ID and position.

            public:

#pragma region • constructors / destructor
                starFinder(const std::shared_ptr<FITScube>& skyImage): fsky(skyImage)
                {
                    if(fsky == nullptr)
                        throw FITSexception(BAD_DIMEN,"DST","starFinder::starFinder","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");
                }

                ~starFinder() {fsky.reset();};
#pragma endregion
#pragma region • modifiers
                void setDetectionMethod(const detectionMethod& m) {fmethod = m;}
#pragma endregion
#pragma region • accessors
                detectionMethod getDetectionMethod() const {return fmethod;}
                std::string getDetectionMethodString() const 
                {
                    switch(fmethod)
                    {
                        case detectionMethod::MAXIMA:
                            return "MAXIMA";
                        case detectionMethod::MORPHO:
                            return "MORPHO";
                        case detectionMethod::THRESHOLD:
                            return "THRESHOLD";
                        case detectionMethod::SNR:
                            return "SNR";
                        default:
                            return "UNKNOWN";
                    }
                }
#pragma endregion
        };
    }
}