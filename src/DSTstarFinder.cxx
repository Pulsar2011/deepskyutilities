//
//  DSTstarFinder.cxx
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#include "DSTstarFinder.h"

namespace DST
{
#pragma region - DST::starFinder class implementation

    /** 
     * @brief Default constructor
     */
    Stars::Stars(): Stars(0)
    { }

    /** 
     * @brief Constructor with star ID
     * @param id Star identifier
     */
    Stars::Stars(const size_t& id): fID(id), fevalMethod(evalType::MOMENT), fposition(), fMaxima(std::numeric_limits<double>::quiet_NaN), fMaximaPosition(), fFWHM(), fangle(std::numeric_limits<double>::quiet_NaN), fTotalFlux(std::numeric_limits<double>::quiet_NaN), fSNR(std::numeric_limits<double>::quiet_NaN), fBackground(std::numeric_limits<double>::quiet_NaN)
    { }

    Stars::Stars(const size_t& id, const pixelCoords& pos): Stars(id)
    {
        if(pos.size() < 2)
            throw std::invalid_argument("DST::Stars::Stars: Pixel coordinates size must be at least 2 (x,y).["+std::to_string(__LINE__)+"]");
        
        fposition.x = pos[0];
        fposition.y = pos[1];
    }

    Stars::Stars(const size_t& id, const point& pos): Stars(id)
    {
        if(pos.size() < 2)
            throw std::invalid_argument("DST::Stars::Stars: Pixel coordinates size must be at least 2 (x,y).["+std::to_string(__LINE__)+"]");
        
        fposition.x = pos[0];
        fposition.y = pos[1];
    }

    Stars::~Stars()
    {
        fpixels.clear();
    }
#pragma endregion
    namespace algorithms
    {

    }
    
} // namespace name
