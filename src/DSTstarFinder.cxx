//
//  DSTstarFinder.cxx
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#include <DSTools/DSTstarFinder.h>

#include <queue>
#include <set>
#include <mutex>
#include <cstdlib>

#ifdef __WITH_QT__

#include <DSTdisplay/QtFITSviewer.h>

#include <QtConcurrent/QtConcurrent>
#include <QMetaObject>
#include <QMainWindow>
#include <QStatusBar>
#include <QSlider>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QProxyStyle>
#include <QPainter>
#include <QStyleOption>
#include <QPolygon>
#include <QToolTip>
#include <QCursor>

#endif

#include <cmath>
#include <chrono>
#include <ctime>


std::mutex flabelMutex;
/// Serializes concurrent GPU folding calls — the Vulkan backend is not thread-safe.
std::mutex fgpuMutex;

// -----------------------------------------------------------------------
//  File-static GPU helpers
// -----------------------------------------------------------------------

/// Convert a kernel<From> to kernel<To> element-wise (for GPU dispatch).
template<typename To, typename From>
static DST::Morphology::kernel<To>
convertKernel(const DST::Morphology::kernel<From>& k)
{
    std::valarray<To> w(k.size());
    for(size_t i = 0; i < k.size(); ++i)
        w[i] = static_cast<To>(k.data()[i]);
    return DST::Morphology::kernel<To>(w, k.sizes(), k.type());
}

/// Returns true when the image's pixel storage type is GPU-eligible.
/// Uses GetEqBitPerPixel() which reflects the effective type after
/// BSCALE/BZERO (via cfitsio fits_get_img_equivtype).
/// double (DOUBLE_IMG=-64), int64 (LONGLONG_IMG=64), uint64 (ULONGLONG_IMG=80)
/// are NOT GPU-eligible; all others (int8/16/32 + float) are.
static bool isGPUEligible(const std::shared_ptr<DSL::FITScube>& img)
{
    int bpp = img->GetEqBitPerPixel();
    return bpp != DOUBLE_IMG && bpp != LONGLONG_IMG && bpp != ULONGLONG_IMG;
}

namespace DST
{
#pragma region - DST::Stars class implementation
#pragma region -- ctor/dtor
    /** 
     * @brief Default constructor
     */
    Stars::Stars(): Stars(0)
    { }

    /** 
     * @brief Constructor with star ID
     * @param id Star identifier
     */
    Stars::Stars(const size_t& id):
        fID(id),
        fevalMethod(evalType::MOMENT),
        fposition(),
        fPSFparameters(),
        fMaxima(std::numeric_limits<double>::quiet_NaN()),
        fMaximaPosition(),
        fFWHM(),
        fangle(std::numeric_limits<double>::quiet_NaN()),
        fTotalFlux(std::numeric_limits<double>::quiet_NaN()),
        fSNR(std::numeric_limits<double>::quiet_NaN()),
        fBackground(std::numeric_limits<double>::quiet_NaN()),
        finURL()
    { }

    Stars::Stars(const size_t& id, const pixelCoords& pos): Stars(id)
    {
        if(pos.size() < 2)
            throw std::invalid_argument("DST::Stars::Stars: Pixel coordinates size must be at least 2 (x,y).["+std::to_string(__LINE__)+"]");
        
        point tmp = point::pointFromVector<double>(pos);
        fposition.SetPoints(tmp);
    }

    Stars::Stars(const size_t& id, const point& pos): Stars(id)
    {
        if(pos.size() < 2)
            throw std::invalid_argument("DST::Stars::Stars: Pixel coordinates size must be at least 2 (x,y).["+std::to_string(__LINE__)+"]");
        
        fposition.SetPoints(pos);
    }

    Stars::~Stars()
    {
        fpixels.clear();
        fPSFparameters.clear();
        finURL.clear();
    }

#pragma endregion
#pragma region -- Modifiers

    /**
     * @brief Insert pixel to the star pixel list
     * 
     * @param pix pixel coordinate vector
     */
    void Stars::addPixel(const pixelCoords& pix)
    {
        if(pix.size() < 2)
            throw std::invalid_argument("DST::Stars::addPixel: Pixel coordinates size must be at least 2 (x,y).["+std::to_string(__LINE__)+"]");
        
        point p;
        p.SetPoints(pix);

        fpixels.push_back(p);
    }

    /**
     * @brief Insert multiple pixels to the star pixel list
     * 
     * @param pixs array of pixel coordinate vectors
     */
    void Stars::addPixel(const pixelVectors& pixs)
    {
        for(pixelVectors::const_iterator pix=pixs.cbegin(); pix!=pixs.cend(); ++pix)
        {
            addPixel(*pix);
        }
    }

    
    /**
     * @brief Insect pixel to the star pixel list
     * 
     * @param pix Coordinates of the pixel to add
     */
    void Stars::addPixel(const point& pix)
    {
        if(pix.size() < 2)
            throw std::invalid_argument("DST::Stars::addPixel: Pixel coordinates size must be at least 2 (x,y).["+std::to_string(__LINE__)+"]");
        
        fpixels.push_back(point(pix));
    }

    /**
     * @brief Insert pixel to the star pixel list
     * 
     * @param pix vector of coordinates of pixels to be added
     */
    void Stars::addPixel(const std::vector<point>& pix)
    {
        for(std::vector<point>::const_iterator p=pix.cbegin(); p!=pix.cend(); ++p)
            fpixels.push_back(*p);
    }

    /**
     * @brief Set centroid position of the star
     * 
     * @param p vector ofpixel coordinates (x,y). 
     * @note Pixel coordinates ar pixel along image axis in pixel coordinates. This isn't WCS coordinat (still will works if given WCS but algorithm will have undefined behaviour)
     */
    void Stars::setPosition(const pixelCoords& p)
    {
        fposition.SetPoints(p);
    }

    /**
     * @brief Set centroid position of the star
     * 
     * @param p Pixel coordinates of the star.
     */
    void Stars::setPosition(const point& p)
    {
        fposition.SetPoints(p);
    }

    /**
     * @brief Set the maximum value of the star
     * 
     * @param max maximum vlaue in the stellar area
     */
    void Stars::setMaxima(const double& max)
    {
        fMaxima = max;
    }

    /**
     * @brief Set the position of the maximum value of the star
     * 
     */
    void Stars::setMaximaPosition(const pixelCoords& maxPos)
    {
        fMaximaPosition.SetPoints(maxPos);
    }
    /**
     * @brief Set the position of the maximum value of the star
     * 
     */
    void Stars::setMaximaPosition(const point& maxPos)
    {
        fMaximaPosition.SetPoints(maxPos);
    }

    /**
     * @brief Set the FWHM of the star
     * 
     * @param fwhmX FWHM along X axis
     * @param fwhmY FWHM along Y axis
     */
    void Stars::setFWHM(const double& fwhmX, const double& fwhmY)
    {
        fFWHM.SetPoints({fwhmX,fwhmY});
    }

    /**
     * @brief Set the FWHM of the star
     * 
     * @param fwhm Vector of FWHM values (FWHM_x,FWHM_y)
     */
    void Stars::setFWHM(const point& fwhm)
    {
        fFWHM.SetPoints(fwhm);
    }

    /**
     * @brief Set the orientation of the star
     * 
     * @param angle Stellar PSF orientation in radiant
     */
    void Stars::setAngle(const double& angle)
    {
        fangle = angle;
    }

    /**
     * @brief Set the total flux of the star
     * 
     * @param flux 
     */
    void Stars::setTotalFlux(const double& flux)
    {
        fTotalFlux = flux;
    }

    /**
     * @brief Set the SNR of the star
     * 
     */
    void Stars::setSNR(const double& snr)
    {
        fSNR = snr;
    }

    /**
     * @brief Set the background value used for star evaluation
     * 
     */
    void Stars::setBackground(const double& background)
    {
        fBackground = background;
    }

    /**
     * @brief Define which evaluation method to use when computing star parameters
     * 
     * @param method 
     */
    void Stars::setEvalMethod(const evalType& method)
    {
        fevalMethod = method;
    }

    void Stars::setPSFparameters(const std::vector<double>& params)
    {
        fPSFparameters = params;
    }

#pragma endregion

#pragma region -- Accessor

    const std::string Stars::ModelString(const evalType& m)
    {
        switch (m)
        {
            case evalType::MOMENT:
                return "MOMENT";
            case evalType::GAUSS:
                return "GAUSS-FIT";
            case evalType::iGAUSS:
                return "INT-GAUSS-FIT";
            case evalType::GAUSS2D:
                return "2D-GAUSS-FIT";
            case evalType::iGAUSS2D:
                return "2D-INT-GAUSS-FIT";
            case evalType::MOFFAT:
                return "MOFFAT-FIT";
            case evalType::MOFFAT2D:
                return "2D-MOFFAT-FIT";
            case evalType::LORENTZIAN:
                return "LORENTZIAN-FIT";
            case evalType::LORENTZIAN2D:
                return "2D-LORENTZIAN-FIT";
            case evalType::VOIGT:
                return "VOIGT-FIT";
            case evalType::VOIGT2D:
                return "2D-VOIGT-FIT";
            case evalType::SERSIC:
                return "SERSIC-FIT";
            case evalType::SERSIC2D:
                return "2D-SERSIC-FIT";
            case evalType::KING:
                return "KING-FIT";
            case evalType::KING2D:
                return "2D-KING-FIT";
            case evalType::AIRY:
                return "AIRY-FIT";
            default:
                return "UNKNOWN";
        }
    }

    /** 
     * @details Get the evaluation method as string
     * @return \c std::string Evaluation method string
     */
    const std::string Stars::evalMethodString() const
    {
        return ModelString(fevalMethod);
    }

#pragma endregion

#pragma region -- Dumper
    /**
     * @brief Printout stellar properties
     * @details Print to the given output stream the stellar properties
     * 
     * @param os 
     */
    void Stars::dump(std::ostream& os, bool _header_) const
    {
        bool is3D = false;
        if(fposition.size() >=3)
            is3D = true;

        if(_header_)
            os << "#ID ; X ; Y ; "<<((is3D)?"Z ; ":"")<<"Maxima ; Maxima_X ; Maxima_Y ; FWHM_X ; FWHM_Y ; Angle(rad) ; TotalFlux ; SNR ; Background ; EvalMethods" << std::endl;

        os << fID << " ; "
              << fposition.X() << " ; "
              << fposition.Y() << " ; "
              << ((is3D)? std::to_string(static_cast<size_t>(fposition.Z()))+" ; " : "")
              << fMaxima << " ; "
              << fMaximaPosition.X() << " ; "
              << fMaximaPosition.Y() << " ; "
              << ( (fFWHM.size()>=2) ? fFWHM.X() : std::numeric_limits<double>::quiet_NaN() ) << " ; "
              << ( (fFWHM.size()>=2) ? fFWHM.Y() : std::numeric_limits<double>::quiet_NaN() ) << " ; "
              << fangle << " ; "
              << fTotalFlux << " ; "
              << fSNR << " ; "
              << fBackground<<" ; "
              << evalMethodString()<<" ; "
              << finURL
              << std::endl;
    }

#pragma endregion

#pragma endregion

#pragma region - Algorithms 
    namespace algorithms
    {

#pragma region - starFitter class implementation

#pragma region -- private member function

        std::shared_ptr< const FITSimg<double> > starFitter::_makeWeight(const std::shared_ptr< FITSimg<double> >& err, bool isErrorMap)
        {       
            if (!isErrorMap)
                return std::static_pointer_cast< const FITSimg<double> >( err );  // pass-through: user gave a weight map directly

            // build weight = 1/err^2 without touching the user's array
            std::shared_ptr< FITSimg<double> > dblErr   = std::static_pointer_cast<FITSimg<double> >(err->cast<double>());
            std::shared_ptr< FITSimg<double> > weightMap= std::make_shared<FITSimg<double> >(FITSimg<double>({dblErr->Size(1), dblErr->Size(2)}));
            
            *(weightMap->GetData<double>()) += 1;

            (*weightMap) /= (*dblErr);
            (*weightMap) /= (*dblErr);

            dblErr.reset();  // free the intermediate double-cast error map
            return weightMap;
        }

        std::shared_ptr< const FITSimg<double> > starFitter::_makeWeight(const std::shared_ptr<FITScube>& err, bool isErrorMap)
        {       
            if (!isErrorMap)
                return std::static_pointer_cast< const FITSimg<double> >( err->cast<double>() );  // pass-through: user gave a weight map directly

            // build weight = 1/err^2 without touching the user's array
            std::shared_ptr< FITSimg<double> > dblErr   = std::static_pointer_cast<FITSimg<double> >(err->cast<double>());
            std::shared_ptr< FITSimg<double> > weightMap= std::make_shared<FITSimg<double> >(FITSimg<double>({dblErr->Size(1), dblErr->Size(2)}));
            
            *(weightMap->GetData<double>()) += 1;

            (*weightMap) /= (*dblErr);
            (*weightMap) /= (*dblErr);

            dblErr.reset();  // free the intermediate double-cast error map
            return weightMap;
        }
#pragma endregion


#pragma region -- ctor/dtor
        starFitter::starFitter(const std::shared_ptr<FITScube>& sky, const DST::Stars::evalType& method):
                    ROOT::Minuit2::FCNBase(),
                    fsky(std::static_pointer_cast<const FITSimg<double> >(sky ? sky->cast<double>() : nullptr)),
                    fweight(nullptr),
                    fevalMethod(method)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube pointer is null. ["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube must be 2D or 3D.["+std::to_string(__LINE__)+"]");
        }

        starFitter::starFitter(const std::shared_ptr< FITSimg<double> >& sky, const DST::Stars::evalType& method):
                    ROOT::Minuit2::FCNBase(),
                    fsky(std::static_pointer_cast<const FITSimg<double> >(sky)),
                    fweight(nullptr),
                    fevalMethod(method)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube pointer is null. ["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube must be 2D or 3D.["+std::to_string(__LINE__)+"]");
        }

        starFitter::starFitter(const std::shared_ptr<FITScube>& sky, const std::shared_ptr<FITScube>& err, const DST::Stars::evalType& method, bool isErrorMap):
                    ROOT::Minuit2::FCNBase(),
                    fsky(std::static_pointer_cast<const FITSimg<double> >(sky ? sky->cast<double>() : nullptr)),
                    fweight(_makeWeight(err,isErrorMap)),
                    fevalMethod(method)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube pointer is null. ["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube must be 2D or 3D. ["+std::to_string(__LINE__)+"]");

            if(fweight == nullptr || fweight->GetDimension() != fsky->GetDimension() || fweight->Size(1) != fsky->Size(1) || fweight->Size(2) != fsky->Size(2))
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input error FITScube pointer is null. ["+std::to_string(__LINE__)+"]");
        }

        starFitter::starFitter(const std::shared_ptr< FITSimg<double> >& sky, const std::shared_ptr< FITSimg<double> >& err, const DST::Stars::evalType& method, bool isErrorMap):
                    ROOT::Minuit2::FCNBase(),
                    fsky(std::static_pointer_cast<const FITSimg<double> >(sky)),
                    fweight(_makeWeight(err,isErrorMap)),
                    fevalMethod(method)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube pointer is null. ["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube must be 2D or 3D. ["+std::to_string(__LINE__)+"]");

            if(fweight == nullptr || fweight->GetDimension() != fsky->GetDimension() || fweight->Size(1) != fsky->Size(1) || fweight->Size(2) != fsky->Size(2))
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input error FITScube pointer is null. ["+std::to_string(__LINE__)+"]");
        }

        starFitter::starFitter(const std::shared_ptr< FITSimg<double> >& sky, const std::shared_ptr< FITScube >& err, const DST::Stars::evalType& method, bool isErrorMap):
                    ROOT::Minuit2::FCNBase(),
                    fsky(std::static_pointer_cast<const FITSimg<double> >(sky)),
                    fweight(_makeWeight(err,isErrorMap)),
                    fevalMethod(method)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube pointer is null. ["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube must be 2D or 3D. ["+std::to_string(__LINE__)+"]");

            if(fweight == nullptr || fweight->GetDimension() != fsky->GetDimension() || fweight->Size(1) != fsky->Size(1) || fweight->Size(2) != fsky->Size(2))
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input error FITScube pointer is null. ["+std::to_string(__LINE__)+"]");
        }

        starFitter::starFitter(const std::shared_ptr< FITScube >& sky, const std::shared_ptr< FITSimg<double> >& err, const DST::Stars::evalType& method, bool isErrorMap):
                    ROOT::Minuit2::FCNBase(),
                    fsky(std::static_pointer_cast<const FITSimg<double> >(sky ? sky->cast<double>() : nullptr)),
                    fweight(_makeWeight(err,isErrorMap)),
                    fevalMethod(method)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube pointer is null. ["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input FITScube must be 2D or 3D. ["+std::to_string(__LINE__)+"]");

            if(fweight == nullptr || fweight->GetDimension() != fsky->GetDimension() || fweight->Size(1) != fsky->Size(1) || fweight->Size(2) != fsky->Size(2))
                throw FITSexception(BAD_DIMEN,"DST::starFitter","starFitter","Input error FITScube pointer is null. ["+std::to_string(__LINE__)+"]");
        }

        starFitter::~starFitter() = default;

#pragma endregion

        double starFitter::operator()(const std::vector<double>& params) const
        {
            double chi2 = 0.0;

            for(size_t i = 0 ; i < fsky->Size(); i++)
            {
                if(fsky->isMasked(i))
                    continue;

                const auto   pix    = fsky->PixelCoordinates(i);
                const double x      = static_cast<double>(pix[0]);
                const double y      = static_cast<double>(pix[1]);
                const double val    = (*fsky)[i];
                const double weight = (fweight) ? (*fweight)[i] : 1.0;

                double model_val;
                try {
                    model_val = Eval(x, y, params);
                } catch (...) {
                    return 1e30; // invalid parameters — steer Minuit away
                }
                chi2 += std::pow(val - model_val, 2) * weight;
            }

            return chi2;
        }

        double starFitter::user_eval(double /*x*/, double /*y*/, const std::vector<double>& /*p*/) const
        {
            std::cerr << "\033[1;35;47m[DST::starFitter::user_eval] Warning:\033[0m No user_eval defined for UNDEFINED evaluation method.\n"
                         "                               To use this method, create a child class and override \033[32mstarFitter::user_eval\033[0m.\n";
            throw std::invalid_argument("\033[1;35;47m[DST::starFitter::user_eval] :\033[0m Unknown evaluation method. ["+std::to_string(__LINE__)+"]");
            
            return 0.0;
        }

        double starFitter::Eval(double x, double y, const std::vector<double>& params) const
        {
            switch(fevalMethod)
            {
                case DST::Stars::evalType::GAUSS:
                    //                                  x  y  A          x0         y0         sigma
                    return DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3])
                           + params[4];

                case DST::Stars::evalType::iGAUSS:
                    //                                   x  y  A          x0         y0         sigma      dx  dy
                    return DST::Math::function::igauss2D(x, y, params[0], params[1], params[2], params[3], 1., 1.)
                           + params[params.size()-1];

                case DST::Stars::evalType::GAUSS2D:
                    //                                  x  y  A          x0         y0         sigma_x    sigma_y    theta
                    return DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5])
                           + params[params.size()-1];

                case DST::Stars::evalType::iGAUSS2D:
                    //                                  x  y  A          x0         y0         sigma_x    sigma_y    theta      dx  dy
                    return DST::Math::function::igauss2D(x, y, params[0], params[1], params[2],params[3], params[4], params[5], 1., 1.)
                           + params[params.size()-1];

                case DST::Stars::evalType::MOFFAT:
                    //                                   x  y  A          x0         y0        alpha      beta
                    return DST::Math::function::moffat2D(x, y, params[0], params[1], params[2],params[3], params[4])
                           + params[params.size()-1];

                case DST::Stars::evalType::MOFFAT2D:
                    //                                   x  y  A          x0         y0        alpha      beta       gamma      theta
                    return DST::Math::function::moffat2D(x, y, params[0], params[1], params[2],params[3], params[4], params[5], params[6])
                           + params[params.size()-1];

                case DST::Stars::evalType::LORENTZIAN:
                    // params: A, x0, y0, gamma, bg
                    return DST::Math::function::lorentzian2D(x, y, params[0], params[1], params[2], params[3])
                           + params[4];

                case DST::Stars::evalType::LORENTZIAN2D:
                    // params: A, x0, y0, gamma, q, theta, bg
                    return DST::Math::function::lorentzian2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5])
                           + params[params.size()-1];

                case DST::Stars::evalType::VOIGT:
                    // params: A, x0, y0, fG, fL, bg
                    return DST::Math::function::pseudo_voigt2D(x, y, params[0], params[1], params[2], params[3], params[4])
                           + params[5];

                case DST::Stars::evalType::VOIGT2D:
                    // params: A, x0, y0, fG, fL, q, theta, bg
                    return DST::Math::function::pseudo_voigt2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5], params[6])
                           + params[params.size()-1];

                case DST::Stars::evalType::SERSIC:
                    // params: A, x0, y0, re, n, bg
                    return DST::Math::function::sersic2D(x, y, params[0], params[1], params[2], params[3], params[4])
                           + params[5];

                case DST::Stars::evalType::SERSIC2D:
                    // params: A, x0, y0, re, n, q, theta, bg
                    return DST::Math::function::sersic2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5], params[6])
                           + params[params.size()-1];

                case DST::Stars::evalType::KING:
                    // params: A, x0, y0, rc, rt, bg
                    return DST::Math::function::king2D(x, y, params[0], params[1], params[2], params[3], params[4])
                           + params[5];

                case DST::Stars::evalType::KING2D:
                    // params: A, x0, y0, rc, rt, q, theta, bg
                    return DST::Math::function::king2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5], params[6])
                           + params[params.size()-1];

                case DST::Stars::evalType::AIRY:
                    // params: A, x0, y0, a, bg
                    return DST::Math::function::airy2D(x, y, params[0], params[1], params[2], params[3])
                           + params[4];

                default:
                    return user_eval(x, y, params);
            }
        }

#pragma endregion

#pragma region - starFinder class implementation
#pragma region -- ctor/dtor

        /**
         * @brief Default constructor
         */
        starFinder::starFinder(const std::shared_ptr<FITScube>& skyImage):
                                    fsky(skyImage),
                                    flabel(nullptr),
                                    fmethod(detectionMethod::MORPHO),fdetectedStars(),
                                    fthresholdValue(std::numeric_limits<double>::min()),
                                    fBS(Morphology::kernel<double>::one({3,3})),
                                    fBE(Morphology::kernel<double>::one({9,9})),
                                    fBM(Morphology::kernel<double>::one({13,13})),
                                    fbackend_(DST::GPU::Backend::create()),
                                    foutURL(""),
                                    is3D(false)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::starFinder","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() < 2 || fsky->GetDimension() > 3)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::starFinder","Input FITScube must be 2D or 3D.["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() == 3)
                is3D = true;

            // Set default threshold.
            fthresholdValue=fsky->Getpercentil(0.97);

            std::valarray<double> bmData(1,13*13);
            bmData[std::gslice(2+2*13,{13-4,13-4},{1,13})] = 0;
            fBM *= bmData;

            fBS.setType(Morphology::ktype::OPEN);
            fBE.setType(Morphology::ktype::ERODE);
            fBM.setType(Morphology::ktype::DILATE);
        }

        starFinder::~starFinder()
        {
            fdetectedStars.clear();
            flabel.reset();
        }
#pragma endregion
#pragma region -- modifiers
        /**
         * @brief Set the output URL for star finder results
         * @param url Output URL (directory path)
         * @throw std::invalid_argument if URL is too short
         * @note The output url is mostly used for debug and saving intermediate product of the starfinder. This function creates the output directory if it does not exist
         */
        void starFinder::setOutputUrl(const std::string& url)
        {
            // Ensure output directory exists
            if(url.size() > 3)
            {
                std::filesystem::path outDir(url);
                std::error_code ec;
                std::filesystem::create_directories(outDir, ec);
            }
            else
            {
                throw std::invalid_argument("DST::starFinder::setOutputUrl: Output URL is too short.["+std::to_string(__LINE__)+"]");
            }

            foutURL = url;
        }

        /**
         * @brief Get the detection method as string
         * 
         * @return std::string Detection method string
         */
        std::string starFinder::detectionMethodString() const
        {
            switch (fmethod)
            {
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

        /** 
         * @details Get the list of identified stars
         */
        DSL::FITStable starFinder::detectedStars() const
        {
            return detectedStars("DETECTED.STARS");
        }

        /** 
         * @details Get the list of identified stars
         * @param tblName Name of the output FITStable
         * @return \c DSL::FITStable Table containing the list of detected stars
         */
        DSL::FITStable starFinder::detectedStars(const std::string& tblName) const
        {

            std::string pixel_units = "ADU";
            if(fsky->HDU().Exists("BUNIT"))
            {
                pixel_units = fsky->HDU().GetValueForKey("BUNIT");
            }

            std::vector<std::string> wcsCoord;
            std::vector<std::string> wcsUnits;
            size_t k=1;
            while(fsky->HDU().Exists("CTYPE"+std::to_string(k)))
            {
                wcsCoord.push_back(fsky->HDU().GetValueForKey("CTYPE"+std::to_string(k)));

                if(fsky->HDU().Exists("CUNIT"+std::to_string(k)))
                    wcsUnits.push_back(fsky->HDU().GetValueForKey("CUNIT"+std::to_string(k)));
                else
                    wcsUnits.push_back("");

                    k++;
            }

            DSL::FITStable table(tblName);
            size_t ncol = 0;
            FITScolumn<uint64_t>  ID      ("ID", tulonglong, "", ncol); ncol++;
            FITScolumn<std::vector<double>> coo    ("coord_px", tdouble, "pixel", ncol); ncol++;
            FITScolumn<std::vector<double>> coo_max("maxcoord_px", tdouble, "pixel", ncol); ncol++;

            std::vector< FITScolumn<double> > wcsColumns;
            for(size_t i=0; i<wcsCoord.size(); i++)
            {
                wcsColumns.push_back( FITScolumn<double>( wcsCoord[i], tdouble, wcsUnits[i], ncol)); ncol++;
            }

            FITScolumn<std::vector<double>>  params("fits", tdouble, "", ncol); ncol++;

            FITScolumn<double>       maxima        ("maxima", tdouble, pixel_units, ncol); ncol++;
            FITScolumn<double>       FWHM_X        ("FWHM_X", tdouble, "pixel", ncol); ncol++;
            FITScolumn<double>       FWHM_Y        ("FWHM_Y", tdouble, "pixel", ncol); ncol++;
            FITScolumn<double>       angle         ("angle", tdouble, "radian", ncol); ncol++;
            FITScolumn<double>       TotalFlux     ("TotalFlux", tdouble, pixel_units, ncol); ncol++;
            FITScolumn<double>       SNR           ("SNR", tdouble, "", ncol); ncol++;
            FITScolumn<double>       Background    ("Background", tdouble, pixel_units, ncol); ncol++;
            FITScolumn<std::string>  EvalMethod    ("EvalMethod", tstring, "", ncol); ncol++;
            FITScolumn<uint32_t>     n_pixels      ("n_pixels", tuint, "", ncol); ncol++;
            FITScolumn<std::string>  source        ("source", tstring, "", ncol); ncol++;

            table.InsertColumn(std::make_shared< FITScolumn<uint64_t> >(ID));
            table.InsertColumn(std::make_shared< FITScolumn<std::vector<double>> >(coo));
            table.InsertColumn(std::make_shared< FITScolumn<std::vector<double>> >(coo_max));

            for(size_t i=0; i<wcsColumns.size(); i++)
                table.InsertColumn(std::make_shared< FITScolumn<double> >(wcsColumns[i]));

            table.InsertColumn(std::make_shared< FITScolumn<std::vector<double>> >(params));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(maxima));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(FWHM_X));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(FWHM_Y));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(angle));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(TotalFlux));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(SNR));
            table.InsertColumn(std::make_shared< FITScolumn<double> >(Background));
            table.InsertColumn(std::make_shared< FITScolumn<std::string> >(EvalMethod));
            table.InsertColumn(std::make_shared< FITScolumn<uint32_t> >(n_pixels));
            table.InsertColumn(std::make_shared< FITScolumn<std::string> >(source));

            for(const auto& [id, star] : fdetectedStars)
            {
                size_t newRow = table.push_back("ID", static_cast<uint64_t>(star.ID()));

                DSL::FITSform::doubleVector pos;
                for(size_t dim = 0; dim < star.position().size(); dim++)
                    pos.push_back(star.position()[dim]);
                table.updateAt("coord_px", newRow, pos);

                DSL::FITSform::doubleVector maxpos;
                for(size_t dim = 0; dim < star.maximaPosition().size(); dim++)
                    maxpos.push_back(star.maximaPosition()[dim]);
                table.updateAt("maxcoord_px", newRow, maxpos);

                // WCS coordinates
                if(wcsColumns.size() > 0)
                {
                    worldCoords wcsPos = fsky->WorldCoordinates(star.position().GetPoint());
                    for(size_t i=0; i<wcsColumns.size(); i++)
                    {
                        if(i < wcsPos.size())
                            table.updateAt(wcsColumns[i].getName(), newRow, wcsPos[i]);
                    }
                }


                if(star.PSFparameters().size() > 0)
                    table.updateAt("fits", newRow, star.PSFparameters());

                table.updateAt("maxima", newRow, star.maxima());
                table.updateAt("FWHM_X", newRow, (star.FWHM().size()>=1)? star.FWHM().X() : std::numeric_limits<double>::quiet_NaN());
                table.updateAt("FWHM_Y", newRow, (star.FWHM().size()>=2)? star.FWHM().Y() : std::numeric_limits<double>::quiet_NaN());
                table.updateAt("angle", newRow, star.angle());
                table.updateAt("TotalFlux", newRow, star.totalFlux());
                table.updateAt("SNR", newRow, star.SNR());
                table.updateAt("Background", newRow, star.background());
                table.updateAt("EvalMethod", newRow, star.evalMethodString());
                table.updateAt("n_pixels", newRow, static_cast<uint32_t>(star.FilledPixels().size()));
                table.updateAt("source", newRow, star.srcFile());
            }

            return table;
        }

        /**
         * @brief Write the detected star catalog to a FITS-LDAC file compatible with SCAMP.
         * @param filename Output FITS file path.
         * @param method   PSF evaluation method to use. Pass Stars::evalType::UNDEFINED (default)
         *                 to skip re-evaluation and assume properties are already measured.
         *                 Any other method will re-evaluate (in parallel) every star whose stored
         *                 evalMethod() differs from the requested one.
         */
        void starFinder::toFITS_LDAC(const std::string& filename, Stars::evalType method)
        {
            // ------------------------------------------------------------------
            // 1.  Optionally re-evaluate stellar properties in parallel.
            //     UNDEFINED → caller guarantees properties are ready; skip.
            //     Any other method → re-evaluate stars whose stored method differs.
            // ------------------------------------------------------------------
            if (method != Stars::evalType::UNDEFINED)
            {
                std::vector<std::future<void>> futures;
                for (auto& [id, star] : fdetectedStars)
                {
                    if (star.evalMethod() != method)
                    {
                        futures.push_back(
                            std::async(std::launch::async,
                                       [this, &star, method]()
                                       { EvalStellarProperties(star, method); }));
                    }
                }
                for (auto& f : futures) f.get();
            }

            // ------------------------------------------------------------------
            // 2.  Create the FITS file and write primary HDU keywords.
            // ------------------------------------------------------------------
            DSL::FITSmanager fm = DSL::FITSmanager::Create(filename, /*replace=*/true);

            {
                auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                char buf[32] = {};
                std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", std::gmtime(&t));
                fm.AppendKeyToPrimary("ORIGIN",
                    DSL::FITSkeyword("DeepSkyTools::DSTastrometry", "Library origin",   DSL::key_type::fChar));
                fm.AppendKeyToPrimary("SOFTWARE",
                    DSL::FITSkeyword(DST::gGIT::version(),           "Software version", DSL::key_type::fChar));
                fm.AppendKeyToPrimary("DATE",
                    DSL::FITSkeyword(std::string(buf),               "Creation date",    DSL::key_type::fChar));
            }

            // ------------------------------------------------------------------
            // 3.  LDAC_IMHEAD: 1 row, 1 column = full source image header.
            // ------------------------------------------------------------------
            std::string rawHeader = fsky->HDU().asString();
            while (rawHeader.size() % 80 != 0) rawHeader += ' ';
            const int64_t hdrLen = static_cast<int64_t>(rawHeader.size());

            DSL::FITStable imhead("LDAC_IMHEAD");
            imhead.InsertColumn(std::make_shared< DSL::FITScolumn<std::string> >(
                "Field Header Card", DSL::dtype::tstring, int64_t(1), hdrLen, std::string(""), size_t(0), /*fixedWidth=*/true));
            imhead.push_back("Field Header Card", rawHeader);

            // ------------------------------------------------------------------
            // 4.  LDAC_OBJECTS: one row per detected star.
            // ------------------------------------------------------------------
            DSL::FITStable objects("LDAC_OBJECTS");
            size_t col = 0;
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<int32_t>  >("NUMBER",           DSL::dtype::tint,   "",      col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("X_IMAGE",           DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("Y_IMAGE",           DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("XWIN_IMAGE",        DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("YWIN_IMAGE",        DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("ERRAWIN_IMAGE",     DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("ERRBWIN_IMAGE",     DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("ERRTHETAWIN_IMAGE", DSL::dtype::tfloat, "deg",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("A_IMAGE",           DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("B_IMAGE",           DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("THETA_IMAGE",       DSL::dtype::tfloat, "deg",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("AWIN_IMAGE",        DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("BWIN_IMAGE",        DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("THETAWIN_IMAGE",    DSL::dtype::tfloat, "deg",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("FLUX_AUTO",         DSL::dtype::tfloat, "ADU",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("FLUXERR_AUTO",      DSL::dtype::tfloat, "ADU",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("MAG_AUTO",          DSL::dtype::tfloat, "mag",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("MAGERR_AUTO",       DSL::dtype::tfloat, "mag",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("BACKGROUND",        DSL::dtype::tfloat, "ADU",   col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("FWHM_IMAGE",        DSL::dtype::tfloat, "pixel", col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("SNR_WIN",           DSL::dtype::tfloat, "",      col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("ELONGATION",        DSL::dtype::tfloat, "",      col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<float>    >("CLASS_STAR",        DSL::dtype::tfloat, "",      col++));
            objects.InsertColumn(std::make_shared< DSL::FITScolumn<int16_t>  >("FLAGS",             DSL::dtype::tshort, "",      col++));

            for (const auto& [id, star] : fdetectedStars)
            {
                const float x_img   = static_cast<float>(star.position().X() + 1.0);
                const float y_img   = static_cast<float>(star.position().Y() + 1.0);
                const float fwhm_x  = (star.FWHM().size() >= 1) ? static_cast<float>(star.FWHM().X()) : 0.f;
                const float fwhm_y  = (star.FWHM().size() >= 2) ? static_cast<float>(star.FWHM().Y()) : fwhm_x;
                const float sigma_x = fwhm_x / 2.355f;
                const float sigma_y = fwhm_y / 2.355f;
                const float snr     = static_cast<float>(star.SNR());
                const float flux    = static_cast<float>(star.totalFlux());
                const float fluxErr = (snr > 0.f) ? flux / snr : 0.f;
                const float mag     = (flux > 0.f) ? static_cast<float>(-2.5 * std::log10(static_cast<double>(flux))) : 99.f;
                const float magErr  = (snr > 0.f) ? static_cast<float>(1.0857 / snr) : 99.f;
                const float errA    = (snr > 0.f) ? sigma_x / snr : 0.f;
                const float errB    = (snr > 0.f) ? sigma_y / snr : 0.f;
                const float theta   = static_cast<float>(star.angle() * 180.0 / M_PI);
                const float elong   = (sigma_y > 0.f) ? sigma_x / sigma_y : 1.f;

                size_t row = objects.push_back("NUMBER", static_cast<int32_t>(star.ID()));
                objects.updateAt("X_IMAGE",           row, x_img);
                objects.updateAt("Y_IMAGE",           row, y_img);
                objects.updateAt("XWIN_IMAGE",        row, x_img);
                objects.updateAt("YWIN_IMAGE",        row, y_img);
                objects.updateAt("ERRAWIN_IMAGE",     row, errA);
                objects.updateAt("ERRBWIN_IMAGE",     row, errB);
                objects.updateAt("ERRTHETAWIN_IMAGE", row, 0.f);
                objects.updateAt("A_IMAGE",           row, sigma_x);
                objects.updateAt("B_IMAGE",           row, sigma_y);
                objects.updateAt("THETA_IMAGE",       row, theta);
                objects.updateAt("AWIN_IMAGE",        row, sigma_x);
                objects.updateAt("BWIN_IMAGE",        row, sigma_y);
                objects.updateAt("THETAWIN_IMAGE",    row, theta);
                objects.updateAt("FLUX_AUTO",         row, flux);
                objects.updateAt("FLUXERR_AUTO",      row, fluxErr);
                objects.updateAt("MAG_AUTO",          row, mag);
                objects.updateAt("MAGERR_AUTO",       row, magErr);
                objects.updateAt("BACKGROUND",        row, static_cast<float>(star.background()));
                objects.updateAt("FWHM_IMAGE",        row, (fwhm_x + fwhm_y) * 0.5f);
                objects.updateAt("SNR_WIN",           row, snr);
                objects.updateAt("ELONGATION",        row, elong);
                objects.updateAt("CLASS_STAR",        row, 1.f);
                objects.updateAt("FLAGS",             row, static_cast<int16_t>(0));
            }

            // ------------------------------------------------------------------
            // 5.  Write extensions and flush.
            // ------------------------------------------------------------------
            fm.InsertTable(std::make_shared<DSL::FITStable>(imhead));
            fm.InsertTable(std::make_shared<DSL::FITStable>(objects));
            fm.Write();
        }

#pragma endregion
#pragma region -- processing
        /**
         * @brief Find stars in the input image using the selected detection method
         *
         */
        void starFinder::findStars(const enum detectionMethod& method, double quantizationStep)
        {
            fmethod = method;
            findStars(quantizationStep);
        }

        /**
         * @brief  Find stars in the input image using the selected detection method
         * 
         */
        void starFinder::findStars(double quantizationStep)
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::findStars","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");
            
            if(fsky->GetDimension() < 2 || fsky->GetDimension() > 3)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::findStars","Input FITScube must be 2D or 3D.["+std::to_string(__LINE__)+"]");

            size_t nLayer = (fsky->GetDimension() == 3) ? fsky->Size(3) : 1;
            
            // Process layers in parallel with maximum 4 threads
            const size_t maxThreads = 4;
            const size_t numThreads = std::min(maxThreads, nLayer);
    
            std::vector<std::future<void>> futures;
            futures.reserve(nLayer);

            for(size_t layerIdx=0; layerIdx<nLayer; ++layerIdx)
            {
                futures.push_back(std::async(std::launch::async, [this, layerIdx, &quantizationStep]()
                {
                    std::shared_ptr<FITScube> layer;
                    if(fsky->GetDimension() == 3)
                        layer = fsky->Layer(layerIdx);
                    else
                        layer = fsky;

                    switch (fmethod)
                    {
                        case detectionMethod::MORPHO:
                            _findStarsMorpho(layer,quantizationStep,layerIdx);
                            break;
                        case detectionMethod::THRESHOLD:
                            _findStarsThreshold(layer,quantizationStep,layerIdx);
                            break;
                        case detectionMethod::SNR:
                            _findStarsSNR(layer,quantizationStep,layerIdx);
                            break;
                        case detectionMethod::NOISECHISEL:
                            _findStarsNoiseChisel(layer,quantizationStep,layerIdx);
                            break;
                        case detectionMethod::CUSTOM:
                            _findStarsCustom(layer,quantizationStep,layerIdx);
                            break;
                        default:
                            layer.reset();
                            throw std::invalid_argument("DST::starFinder::findStars: Unknown detection method.["+std::to_string(__LINE__)+"]");
                    }

                if(fsky->GetDimension() == 3)
                    layer.reset();
                }));

                // Limit concurrent threads
                if(futures.size() >= numThreads)
                {
                    futures.front().get();
                    futures.erase(futures.begin());
                }
            }

            // Wait for remaining tasks
            for(auto& future : futures)
                future.get();

            if(foutURL != "")
            {
                flabel->SetName("label_image");
                flabel->Write(foutURL+"/label_image.fits",true);
            }

                        // Perform basicStellar properties assesment
            for(std::map<size_t, Stars>::iterator it=fdetectedStars.begin(); it!=fdetectedStars.end(); ++it)
            {
                if(it->second.FilledPixels().size() == 0)
                {
                    fdetectedStars.erase(it);
                    continue;
                }

                // Compute star properties here (centroid, flux, FWHM, etc.)
                EvalStellarProperties(it->second, Stars::evalType::MOMENT);
            }

        }


#pragma endregion
#pragma region -- private methods

        /**
         * @brief Find stars using morphological operators
         * 
         */
        void starFinder::_findStarsMorpho(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {
            size_t maxThreads = std::thread::hardware_concurrency();
            
            if(maxThreads == 0)
                maxThreads = 1;

            if(img == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::_findStarsMorpho","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(img->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::_findStarsMorpho","Input FITScube must be 2D.["+std::to_string(__LINE__)+"]");

            // Apply morphological opening to enhance star-like features
            const bool gpuOk = fbackend_->available() && isGPUEligible(img);

            std::shared_ptr<FITScube> iK, iN, iS;

            if(gpuOk)
            {
                // GPU path: the Vulkan backend is NOT thread-safe. When findStars() processes
                // a 3D cube with multiple concurrent threads, all GPU calls must be serialized.
                std::lock_guard<std::mutex> gpuLock(fgpuMutex);
                auto kBS_f = convertKernel<float>(fBS);
                iK = Morphology::GPU::folding<float>(img, kBS_f, *fbackend_);
                iK->MaskPixels(img->GetMask());

                auto kBM_f = convertKernel<float>(fBM);
                iN = Morphology::GPU::folding<float>(img, kBM_f, *fbackend_);
                iN->MaskPixels(img->GetMask());

                auto kBE_f = convertKernel<float>(fBE);
                iS = Morphology::GPU::folding<float>(iN, kBE_f, *fbackend_);
                iS->MaskPixels(img->GetMask());
            }
            else
            {
                // CPU path: run iK and iN in parallel, then iS.
                size_t nT = std::thread::hardware_concurrency();
                if(nT > 2) nT = (nT - 2)/2;
                if(nT == 0) nT = 1;

                auto futureIK = std::async(std::launch::async, [&img, this, nT]()
                    {
                        auto r = Morphology::folding<double>(img, fBS, nT);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });
                auto futureIN = std::async(std::launch::async, [&img, this, nT]()
                    {
                        auto r = Morphology::folding<double>(img, fBM, nT);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });

                iK = futureIK.get();
                iN = futureIN.get();
                iS = Morphology::folding<double>(iN, fBE, maxThreads);
                iS->MaskPixels(img->GetMask());
            }

            if(foutURL != "")
            {
                iK->SetName("BS_opening_layer_"+std::to_string(iLayer));
                iN->SetName("BM_dilation_layer_"+std::to_string(iLayer));
                iS->SetName("BE_erosion_layer_"+std::to_string(iLayer));

                iK->Write(foutURL+"/BS_opening_layer_"+std::to_string(iLayer)+".fits",true);
                iN->Write(foutURL+"/BM_dilation_layer_"+std::to_string(iLayer)+".fits",true);
                iS->Write(foutURL+"/BE_erosion_layer_"+std::to_string(iLayer)+".fits",true);
            }

            // Compute significance map
            //*(iK.get()) -= *(iS.get());
            //double minVal = std::min(iK->GetMinimum(),iS->GetMinimum());
            //*(iK.get()) -=  minVal; 

            double lowestValue = std::numeric_limits<double>::max();
            for(size_t iElem = 0; iElem < iK->Nelements(); ++iElem)
            {
                if(iK->Masked(iElem) || iS->Masked(iElem))
                    continue;

                double valK = iK->DoubleValueAtPixel(iElem);
                double valS = iS->DoubleValueAtPixel(iElem);
                // Guard against non-finite values propagating into significance map
                if(!(std::isfinite(valK) && std::isfinite(valS)))
                {
                    iK->SetPixelValue(std::numeric_limits<double>::lowest(), iElem);
                    iK->MaskPixels({iElem});
                    continue;
                }
                double val = std::min(valK,valS);
                iK->SetPixelValue((double) (valK - val),iElem);

                if((valK - val) < lowestValue)
                    lowestValue = valK - val;
            }

            // If no valid pixels updated lowestValue, avoid subtracting max()
            if(lowestValue == std::numeric_limits<double>::max())
                lowestValue = 0.0;
            *(iK.get()) -=  lowestValue;


            //for(size_t i=0; i<iK->Nelements(); ++i)
            //{
            //    if(iS->DoubleValueAtPixel(i) <= std::numeric_limits<double>::min())
            //        (*(iK->GetData<double>()))[i] *= 0;
            //}

            if(foutURL != "")
            {
                iK->SetName("significance_map_layer_"+std::to_string(iLayer));
                iK->Write(foutURL+"/significance_map_layer_"+std::to_string(iLayer)+".fits",true);
            }

            iS.reset();
            iN.reset();

            // _watershedFlooding uses GetData<double>() internally.
            // If iK is a float image (GPU path), promote to double before calling it.
            if(iK->GetEqBitPerPixel() != DOUBLE_IMG)
            {
                auto iK_d = std::make_shared<FITSimg<double>>(std::vector<size_t>{iK->Size(1), iK->Size(2)});
                for(size_t i = 0; i < iK->Nelements(); ++i)
                    iK_d->SetPixelValue(iK->DoubleValueAtPixel(i), i);
                iK_d->MaskPixels(iK->GetMask());
                iK = iK_d;
            }

            double threshold = iK->Getpercentil(fthresholdValue);
            _watershedFlooding(iK, threshold, fquantumStepQueue, iLayer);
            iK.reset();
 
        }

        void starFinder::_findStarsThreshold(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {

            size_t maxThreads = std::thread::hardware_concurrency();
            
            if(maxThreads == 0)
                maxThreads = 1;

            // The first step consists in smoothing the image to reduce the noise level.
            // NOTE: MEAN and MEDIAN are fold operations; the GPU FOLD shader is not reliable.
            // These always use the CPU parallel path.
            fBS.setType(DST::Morphology::ktype::MEAN);
            fBE.setType(DST::Morphology::ktype::MEDIAN);

            std::shared_ptr<FITScube> smoothedImage, bkgMap;

            {
                // CPU path: run both in parallel.
                size_t nT = std::thread::hardware_concurrency();
                if(nT > 3) nT = (nT - 3)/3;
                if(nT == 0) nT = 1;

                auto futureSI = std::async(std::launch::async, [&img, this, nT]()
                    {
                        auto r = Morphology::folding<double>(img, fBS, nT);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });
                auto futureBK = std::async(std::launch::async, [&img, this, nT]()
                    {
                        auto r = Morphology::folding<double>(img, fBE, nT);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });

                smoothedImage = futureSI.get();
                bkgMap        = futureBK.get();
            }

            // Subtract background
            *(smoothedImage.get()) -= *(bkgMap.get());

            if(foutURL != "")
            {
                smoothedImage->Write(foutURL+"/thresholded_image_layer_"+std::to_string(iLayer)+".fits",true);
            }

            // The second step consists in searching for connected components above the threshold value using a watershed algorithm.

            _watershedFlooding(smoothedImage, fthresholdValue, fquantumStepQueue,iLayer);

            smoothedImage.reset();
        }

        void starFinder::_findStarsSNR(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {
            // Define kernel to be used
            fBS.setType(DST::Morphology::ktype::RMSE);
            fBE.setType(DST::Morphology::ktype::MEDIAN);


            // The first step consists in estimating the background and noise level of the image.
            // NOTE: RMSE, MEDIAN and SUM are fold operations; the GPU FOLD shader is not reliable.
            // These always use the CPU parallel path.
            std::shared_ptr<FITScube> errMap, bkgMap, snrMap;

            {
                // CPU path: run all three in parallel.
                size_t nT = std::thread::hardware_concurrency();
                if(nT > 3) nT = (nT - 3)/3;
                if(nT == 0) nT = 1;

                auto futureErr = std::async(std::launch::async, [&img, this, nT]()
                    {
                        auto r = Morphology::folding<double>(img, fBS, nT);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });
                auto futureBkg = std::async(std::launch::async, [&img, this, nT]()
                    {
                        auto r = Morphology::folding<double>(img, fBE, nT);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });
                auto futureSNR = std::async(std::launch::async, [&img]()
                    {
                        Morphology::kernel<double> k = DST::Morphology::kernel<double>::one({1,1});
                        k.setType(DST::Morphology::ktype::SUM);
                        auto r = Morphology::folding<double>(img, k);
                        r->MaskPixels(img->GetMask());
                        return r;
                    });

                errMap = futureErr.get();
                bkgMap = futureBkg.get();
                snrMap = futureSNR.get();
            }

            // Compute Li & Ma SNR map
            //*(snrMap.get()) -= *(bkgMap.get());
            //*(errMap.get()) += 1e-10;               // Avoid division by zero
            //*(snrMap.get()) /= *(errMap.get());

            for(size_t ik=0; ik<snrMap->Nelements(); ++ik)
            {
                double N_on  = snrMap->DoubleValueAtPixel(ik);
                double N_off = bkgMap->DoubleValueAtPixel(ik);

                if(N_off <= 0 || N_on <= 0)
                {
                    (*(snrMap->GetData<double>()))[ik] = 0.0;
                    continue;
                }

                double term1 = N_on * std::log10( N_on/ N_off);
                double term2 = N_on - N_off;

                double S = std::sqrt(2.0 * term1 - term2);

                (*(snrMap->GetData<double>()))[ik] = S;
            }

            if(foutURL != "")
            {
                errMap->Write(foutURL+"/error_map_layer_"+std::to_string(iLayer)+".fits",true);
                bkgMap->Write(foutURL+"/background_map_layer_"+std::to_string(iLayer)+".fits",true);
                snrMap->Write(foutURL+"/snr_map_layer_"+std::to_string(iLayer)+".fits",true);
            }

            errMap.reset();
            bkgMap.reset();

            // The third step consists in grouping connected pixels above the SNR threshold into star candidates using a watershed algorithm.


            _watershedFlooding(snrMap, fthresholdValue, fquantumStepQueue, iLayer);
            snrMap.reset();
        }

        /**
         * @brief Find extended sources using the NoiseChisel algorithm
         *        (Akhlaghi & Ichikawa 2015, ApJS 220, 1).
         *
         * Pipeline:
         *  1. Gaussian convolution of the input image.
         *  2. Per-tile quantile thresholding → binary detection map.
         *  3. Morphological erosion to remove noise spikes.
         *  4. Connected component labeling of surviving blobs.
         *  5. Per-blob S/N filter: reject blobs below fNC.snrThreshold.
         *  6. Dilation to recover the true extent of surviving detections.
         *  7. Sky background estimation on the source-masked image.
         *  8. Sky subtraction.
         *  9. Optional re-detection on the sky-subtracted image.
         * 10. Watershed flooding to label individual sources.
         *
         * All tunable parameters are read from the fNC struct; none are hard-coded.
         * GPU/CPU dispatch follows the same pattern as _findStarsMorpho().
         */
        void starFinder::_findStarsNoiseChisel(const std::shared_ptr<FITScube>& img,
                                                const double& fquantumStepQueue,
                                                const size_t& iLayer)
        {
            if(img == nullptr)
                throw FITSexception(BAD_DIMEN, "DST",
                    "starFinder::_findStarsNoiseChisel",
                    "Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");
            if(img->GetDimension() != 2)
                throw FITSexception(BAD_DIMEN, "DST",
                    "starFinder::_findStarsNoiseChisel",
                    "Input FITScube must be 2D.["+std::to_string(__LINE__)+"]");

            size_t maxThreads = std::thread::hardware_concurrency();
            if(maxThreads == 0) maxThreads = 1;

            const size_t W    = img->Size(1);
            const size_t H    = img->Size(2);
            const bool   gpuOk = fbackend_->available() && isGPUEligible(img);

            // ------------------------------------------------------------------
            // Inner lambda: runs detection steps 1–6 (convolve → binarise →
            // erode → CCL → S/N filter → dilate) on an arbitrary source image
            // and returns the dilated detection mask (FITScube with values 0/1).
            // ------------------------------------------------------------------
            auto runDetection = [&](const std::shared_ptr<FITScube>& src)
                -> std::shared_ptr<FITScube>
            {
                // ---- Step 1: Gaussian convolution ----------------------------
                const size_t kDim = 2 * fNC.convRadius + 1;
                std::valarray<double> gw(kDim * kDim);
                double gsum = 0.0;
                for(size_t j = 0; j < kDim; ++j)
                {
                    const double dy = static_cast<double>(j) - static_cast<double>(fNC.convRadius);
                    for(size_t i = 0; i < kDim; ++i)
                    {
                        const double dx = static_cast<double>(i) - static_cast<double>(fNC.convRadius);
                        const double v  = std::exp(-(dx*dx + dy*dy)
                                                    / (2.0 * fNC.convSigma * fNC.convSigma));
                        gw[j * kDim + i] = v;
                        gsum += v;
                    }
                }
                gw /= gsum; // normalise manually → use CORRELATE (not nCORRELATE)

                DST::Morphology::kernel<double> gaussKernel(
                    gw, {kDim, kDim}, DST::Morphology::ktype::CORRELATE);

                std::shared_ptr<FITScube> convImg;
                if(gpuOk)
                {
                    std::lock_guard<std::mutex> gpuLock(fgpuMutex);
                    auto gk_f = convertKernel<float>(gaussKernel);
                    convImg = Morphology::GPU::correlate<float>(
                        src, gk_f, *fbackend_, maxThreads, false);
                }
                else
                {
                    convImg = Morphology::correlate<double>(
                        src, gaussKernel, maxThreads, false);
                }
                convImg->MaskPixels(src->GetMask());

                if(foutURL != "")
                {
                    convImg->SetName("nc_convolved_layer_"+std::to_string(iLayer));
                    convImg->Write(foutURL+"/nc_convolved_layer_"
                                   +std::to_string(iLayer)+".fits", true);
                }

                // ---- Step 2: Per-tile quantile threshold map -----------------
                const size_t nTX = (W + fNC.tileSize - 1) / fNC.tileSize;
                const size_t nTY = (H + fNC.tileSize - 1) / fNC.tileSize;
                std::vector<double> tileThresh(nTX * nTY,
                                               std::numeric_limits<double>::max());

                for(size_t ty = 0; ty < nTY; ++ty)
                {
                    for(size_t tx = 0; tx < nTX; ++tx)
                    {
                        const size_t x0 = tx * fNC.tileSize;
                        const size_t y0 = ty * fNC.tileSize;
                        const size_t x1 = std::min(x0 + fNC.tileSize, W);
                        const size_t y1 = std::min(y0 + fNC.tileSize, H);

                        std::vector<double> tileVals;
                        tileVals.reserve((x1-x0)*(y1-y0));
                        for(size_t y = y0; y < y1; ++y)
                            for(size_t x = x0; x < x1; ++x)
                            {
                                const size_t idx = y * W + x;
                                if(!convImg->Masked(idx))
                                    tileVals.push_back(convImg->DoubleValueAtPixel(idx));
                            }

                        if(tileVals.empty()) continue;

                        const size_t qIdx = static_cast<size_t>(
                            fNC.detectionQuantile
                            * static_cast<double>(tileVals.size() - 1));
                        std::nth_element(tileVals.begin(),
                                         tileVals.begin() + qIdx,
                                         tileVals.end());
                        tileThresh[ty * nTX + tx] = tileVals[qIdx];
                    }
                }

                // ---- Step 3: Binarise ----------------------------------------
                auto binary = std::make_shared<FITSimg<uint8_t>>(
                    std::vector<size_t>{W, H});
                std::vector<size_t> binaryMasked;
                
                for(size_t y = 0; y < H; ++y)
                {
                    for(size_t x = 0; x < W; ++x)
                    {
                        const size_t idx = y * W + x;
                        if(src->Masked(idx))
                        {
                            binary->SetPixelValue(0.0, idx);
                            binaryMasked.push_back(idx);
                        }
                        else
                        {
                            const size_t tx     = std::min(x / fNC.tileSize, nTX - 1);
                            const size_t ty     = std::min(y / fNC.tileSize, nTY - 1);
                            const double thresh = tileThresh[ty * nTX + tx];
                            binary->SetPixelValue(
                                (convImg->DoubleValueAtPixel(idx) > thresh) ? 1.0 : 0.0,
                                idx);
                        }
                    }
                }

                
                if(!binaryMasked.empty())
                    binary->MaskPixels(binaryMasked);
                convImg.reset();

                if(foutURL != "")
                {
                    binary->SetName("nc_binary_layer_"+std::to_string(iLayer));
                    binary->Write(foutURL+"/nc_binary_layer_"
                                   +std::to_string(iLayer)+".fits", true);
                }

                // ---- Step 4: Erosion -----------------------------------------
                const size_t eDim = 2 * fNC.erosionRadius + 1;
                std::valarray<double> ew(1.0, eDim * eDim);
                DST::Morphology::kernel<double> erodeKernel(
                    ew, {eDim, eDim}, DST::Morphology::ktype::ERODE);

                std::shared_ptr<FITScube> erodedBinary;
                if(gpuOk)
                {
                    std::lock_guard<std::mutex> gpuLock(fgpuMutex);
                    auto ek_f = convertKernel<float>(erodeKernel);
                    erodedBinary = Morphology::GPU::erosion<float>(
                        binary, ek_f, *fbackend_, maxThreads);
                }
                else
                {
                    erodedBinary = Morphology::erosion<double>(
                        binary, erodeKernel, maxThreads);
                }
                erodedBinary->MaskPixels(src->GetMask());
                binary.reset();

                if(foutURL != "")
                {
                    erodedBinary->SetName("nc_eroded_layer_"+std::to_string(iLayer));
                    erodedBinary->Write(foutURL+"/nc_eroded_layer_"
                                         +std::to_string(iLayer)+".fits", true);
                }

                // ---- Step 5: Connected component labeling -------------------
                // connectedComponents() requires double-typed input for DoubleValueAtPixel.
                std::shared_ptr<FITScube> erodedDbl;
                if(erodedBinary->GetEqBitPerPixel() != DOUBLE_IMG)
                {
                    auto tmp = std::make_shared<FITSimg<double>>(
                        std::vector<size_t>{W, H});
                    for(size_t i = 0; i < erodedBinary->Nelements(); ++i)
                        tmp->SetPixelValue(erodedBinary->DoubleValueAtPixel(i), i);
                    tmp->MaskPixels(erodedBinary->GetMask());
                    erodedDbl = tmp;
                }
                else
                {
                    erodedDbl = erodedBinary;
                }
                erodedBinary.reset();

                std::shared_ptr<FITScube> labelImg =
                    DST::Morphology::connectedComponents(erodedDbl);
                erodedDbl.reset();

                if(foutURL != "")
                {
                    labelImg->SetName("nc_labels_layer_"+std::to_string(iLayer));
                    labelImg->Write(foutURL+"/nc_labels_layer_"
                                     +std::to_string(iLayer)+".fits", true);
                }

                // ---- Step 6: Per-blob S/N filter ----------------------------
                // Collect label → pixel index mapping
                std::map<int32_t, std::vector<size_t>> blobPixels;
                for(size_t i = 0; i < W * H; ++i)
                {
                    if(labelImg->Masked(i)) continue;
                    const int32_t lbl =
                        static_cast<int32_t>(labelImg->DoubleValueAtPixel(i));
                    if(lbl == 0) continue;
                    blobPixels[lbl].push_back(i);
                }

                // Estimate per-tile sky and noise from the source image
                const size_t nTX2 = (W + fNC.tileSize - 1) / fNC.tileSize;
                const size_t nTY2 = (H + fNC.tileSize - 1) / fNC.tileSize;
                DST::Background::skyEstimator snrSkyEst;
                snrSkyEst.setTileSize(fNC.tileSize);
                snrSkyEst.setSigmaClipIter(fNC.skyIterations);
                snrSkyEst.setSigmaClipThreshold(fNC.skyClipSigma);
                snrSkyEst.estimate(src);
                const auto& tileStats = snrSkyEst.tileStats();

                std::set<int32_t> rejectedLabels;
                for(const auto& kv : blobPixels)
                {
                    const int32_t lbl       = kv.first;
                    const auto&   pixels    = kv.second;
                    const double  nPix      = static_cast<double>(pixels.size());

                    // Blob centroid (pixel coords)
                    size_t cx_sum = 0, cy_sum = 0;
                    double sumSig = 0.0;
                    for(size_t idx : pixels)
                    {
                        sumSig  += src->DoubleValueAtPixel(idx);
                        cx_sum  += idx % W;
                        cy_sum  += idx / W;
                    }
                    const double meanSig = sumSig / nPix;
                    const size_t cx = cx_sum / pixels.size();
                    const size_t cy = cy_sum / pixels.size();

                    // Tile containing the centroid
                    const size_t tx  = std::min(cx / fNC.tileSize, nTX2 - 1);
                    const size_t ty  = std::min(cy / fNC.tileSize, nTY2 - 1);
                    const auto&  ts  = tileStats[ty * nTX2 + tx];
                    const double sky   = ts.valid ? ts.sky : 0.0;
                    const double noise = ts.valid ? ts.rms : 1.0;

                    if((meanSig - sky) / noise < fNC.snrThreshold)
                        rejectedLabels.insert(lbl);
                }

                // Build surviving-blob mask
                auto survivingMask = std::make_shared<FITSimg<uint8_t>>(
                    std::vector<size_t>{W, H});
                std::vector<size_t> survMasked;
                for(size_t i = 0; i < W * H; ++i)
                {
                    if(src->Masked(i)) { survMasked.push_back(i); continue; }
                    const int32_t lbl =
                        static_cast<int32_t>(labelImg->DoubleValueAtPixel(i));
                    survivingMask->SetPixelValue(
                        (lbl != 0 && rejectedLabels.find(lbl) == rejectedLabels.end())
                            ? 1.0 : 0.0,
                        i);
                }
                if(!survMasked.empty())
                    survivingMask->MaskPixels(survMasked);
                labelImg.reset();

                if(foutURL != "")
                {
                    survivingMask->SetName("nc_surviving_layer_"+std::to_string(iLayer));
                    survivingMask->Write(foutURL+"/nc_surviving_layer_"
                                          +std::to_string(iLayer)+".fits", true);
                }

                // ---- Step 7: Dilation to recover true source extent ----------
                const size_t dDim = 2 * fNC.dilationRadius + 1;
                std::valarray<double> dw(1.0, dDim * dDim);
                DST::Morphology::kernel<double> dilateKernel(
                    dw, {dDim, dDim}, DST::Morphology::ktype::DILATE);

                std::shared_ptr<FITScube> detectionMask;
                if(gpuOk)
                {
                    std::lock_guard<std::mutex> gpuLock(fgpuMutex);
                    auto dk_f = convertKernel<float>(dilateKernel);
                    detectionMask = Morphology::GPU::dilation<float>(
                        survivingMask, dk_f, *fbackend_, maxThreads);
                }
                else
                {
                    detectionMask = Morphology::dilation<double>(
                        survivingMask, dilateKernel, maxThreads);
                }
                detectionMask->MaskPixels(src->GetMask());
                survivingMask.reset();

                if(foutURL != "")
                {
                    detectionMask->SetName("nc_detection_mask_layer_"+std::to_string(iLayer));
                    detectionMask->Write(foutURL+"/nc_detection_mask_layer_"
                                          +std::to_string(iLayer)+".fits", true);
                }

                return detectionMask;
            }; // end runDetection lambda

            // ------------------------------------------------------------------
            // First detection pass (on the raw input image)
            // ------------------------------------------------------------------
            std::shared_ptr<FITScube> detectionMask = runDetection(img);

            // ------------------------------------------------------------------
            // Step 8: Sky estimation on source-masked image
            // ------------------------------------------------------------------
            auto maskedForSky = std::make_shared<FITSimg<double>>(
                std::vector<size_t>{W, H});
            std::vector<size_t> skyMaskIdx;
            for(size_t i = 0; i < W * H; ++i)
            {
                if(img->Masked(i) || detectionMask->DoubleValueAtPixel(i) > 0.0)
                {
                    maskedForSky->SetPixelValue(0.0, i);
                    skyMaskIdx.push_back(i);
                }
                else
                {
                    maskedForSky->SetPixelValue(img->DoubleValueAtPixel(i), i);
                }
            }
            if(!skyMaskIdx.empty())
                maskedForSky->MaskPixels(skyMaskIdx);
            skyMaskIdx.clear();

            DST::Background::skyEstimator skyEst;
            skyEst.setTileSize(fNC.tileSize);
            skyEst.setSigmaClipIter(fNC.skyIterations);
            skyEst.setSigmaClipThreshold(fNC.skyClipSigma);
            skyEst.estimate(maskedForSky);
            maskedForSky.reset();

            auto skyMap = skyEst.getSkyMap();
            auto rmsMap = skyEst.getRMSMap();

            if(foutURL != "")
            {
                skyMap->SetName("nc_skymap_layer_"+std::to_string(iLayer));
                skyMap->Write(foutURL+"/nc_skymap_layer_"
                               +std::to_string(iLayer)+".fits", true);
                rmsMap->SetName("nc_rmsmap_layer_"+std::to_string(iLayer));
                rmsMap->Write(foutURL+"/nc_rmsmap_layer_"
                               +std::to_string(iLayer)+".fits", true);
            }

            // ------------------------------------------------------------------
            // Step 9: Sky subtraction
            // ------------------------------------------------------------------
            auto skySubtracted = std::make_shared<FITSimg<double>>(
                std::vector<size_t>{W, H});
            for(size_t i = 0; i < W * H; ++i)
            {
                skySubtracted->SetPixelValue(
                    img->Masked(i)
                        ? 0.0
                        : img->DoubleValueAtPixel(i) - skyMap->DoubleValueAtPixel(i),
                    i);
            }
            skySubtracted->MaskPixels(img->GetMask());
            skyMap.reset();
            rmsMap.reset();

            if(foutURL != "")
            {
                skySubtracted->SetName("nc_skysubtracted_layer_"+std::to_string(iLayer));
                skySubtracted->Write(foutURL+"/nc_skysubtracted_layer_"
                                      +std::to_string(iLayer)+".fits", true);
            }

            // ------------------------------------------------------------------
            // Step 10: Optional re-detection on sky-subtracted image.
            // The result replaces the first-pass detection mask so that the
            // watershed only floods within confirmed source footprints.
            // ------------------------------------------------------------------
            if(fNC.reDetect)
                detectionMask = runDetection(skySubtracted);

            // ------------------------------------------------------------------
            // Step 11: Restrict watershed to detected source footprints.
            // Pixels outside the detection mask are masked out so that the
            // watershed cannot create labels for background noise peaks.
            // The threshold is 0.0: after sky subtraction the background is
            // near zero; only positive (source) pixels are flooded.
            // ------------------------------------------------------------------
            std::vector<size_t> outsideFootprint;
            for(size_t i = 0; i < W * H; ++i)
            {
                if(skySubtracted->Masked(i)) continue;
                if(detectionMask->DoubleValueAtPixel(i) == 0.0)
                    outsideFootprint.push_back(i);
            }
            if(!outsideFootprint.empty())
                skySubtracted->MaskPixels(outsideFootprint);
            

            double threshold = skySubtracted->Getpercentil(fthresholdValue);
            _watershedFlooding(skySubtracted, threshold, fquantumStepQueue, iLayer);
            skySubtracted.reset();
            detectionMask.reset();
        }

        void starFinder::_findStarsCustom(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {
            std::shared_ptr<FITScube> customMap = _customPreprocessing(img, iLayer);
            if(customMap == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFinder","_findStarsCustom","Custom preprocessing returned a null FITScube pointer.["+std::to_string(__LINE__)+"]");

            if(foutURL != "")
            {
                customMap->Write(foutURL+"/userSignificance_map_layer_"+std::to_string(iLayer)+".fits",true);
            }

            _watershedFlooding(customMap, fthresholdValue, fquantumStepQueue, iLayer);
            customMap.reset();
        }

        std::shared_ptr<FITScube> starFinder::_customPreprocessing(const std::shared_ptr<FITScube>& img, const size_t& iLayer)
        {
            throw std::invalid_argument("DST::starFinder::_customPreprocessing: Custom preprocessing require to be implemented in derived class.["+std::to_string(__LINE__)+"]");
        }

        /**
         * @brief Watershed flooding algorithm to detect stars
         * 
         */
        void starFinder::_watershedFlooding(const std::shared_ptr<FITScube>& img, const double& _thresh,  const double& fquantumStepQueue, const size_t& iLayer)
        {
            if(img == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::_watershedFlooding","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            size_t width  = img->Size(1);
            size_t height = img->Size(2);

            std::valarray<size_t> labels=std::valarray<size_t>(img->Nelements());
            labels *= 0;

            double maxVal = img->GetMaximum()-fquantumStepQueue;
            const std::valarray<double>* data = img->GetData<double>();

            if((verbose & verboseLevel::VERBOSE_BASIC)== verboseLevel::VERBOSE_BASIC)
            {
                std::cout << "DST::starFinder::_watershedFlooding: Starting watershed flooding on layer "<<iLayer<<" with max value "<<maxVal<<" and threshold value "<<_thresh<<"."<<std::endl;
            }

            for(double threshold=maxVal; threshold >= _thresh; threshold-=fquantumStepQueue)
            {
                // Find pixels above the current threshold
                std::vector<size_t> pixelIndices;

                const size_t N = img->Nelements();
                const size_t maxThreads = std::max<size_t>(1, std::thread::hardware_concurrency());
                const size_t numThreads = std::min<size_t>(maxThreads, 8); // cap if you want
                std::vector<std::thread> threads;
                threads.reserve(numThreads);
                
                // thread-local buffers
                std::vector< std::vector<size_t> > locals(numThreads);
                for (size_t t = 0; t < numThreads; ++t)
                {
                    size_t begin = (N * t) / numThreads;
                    size_t end   = (N * (t + 1)) / numThreads;
                
                    threads.emplace_back([&, t, begin, end] {
                        auto& out = locals[t];
                        out.reserve((end - begin) / 16 + 64); // heuristic
                        for (size_t idx = begin; idx < end; ++idx)
                        {
                            if (!img->Masked(idx) && (*data)[idx] >= threshold && labels[idx] == 0)
                                out.push_back(idx);
                        }
                    });
                }
                for (auto& th : threads) th.join();
            
                // merge in block order to keep ascending order
                size_t total = 0;
                for (auto& v : locals) total += v.size();
                pixelIndices.clear();
                pixelIndices.reserve(total);
                for (auto& v : locals) pixelIndices.insert(pixelIndices.end(), v.begin(), v.end());

                // Process each pixel above the threshold
                for(const size_t& pixIdx : pixelIndices)
                {
                    // Check neighboring pixels for existing labels
                    std::set<size_t> neighborLabels;
                    std::vector<size_t> neighbors = _getNeighborIndices(pixIdx, width, height);
                    for(const size_t& nIdx : neighbors)
                    {
                        if(labels[nIdx] != 0)
                            neighborLabels.insert(labels[nIdx]);
                    }

                    if(neighborLabels.empty())
                    {
                        size_t newLabel = 0;
                        {
                            std::lock_guard<std::mutex> lock(flabelMutex); 

                            // New star detected
                            newLabel += fdetectedStars.size() + 1;
                            labels[pixIdx] = newLabel;
                        
                            fdetectedStars.insert(std::pair<size_t,Stars>(newLabel,Stars(newLabel)));
                        }

                        std::vector<point> belong2Star = _floodingPlateau(pixIdx, labels, img, threshold, newLabel);
                        if(is3D)
                        {
                            for(std::vector<point>::iterator p=belong2Star.begin(); p!=belong2Star.end(); ++p)
                                p->SetZ(static_cast<double>(iLayer));
                        }

                        {
                            std::lock_guard<std::mutex> lock(flabelMutex); 
                            fdetectedStars[newLabel].addPixel(belong2Star);
                        }

                    }
                    else if(neighborLabels.size() == 1)
                    {
                        // Assign existing star label
                        size_t existingLabel = *neighborLabels.begin();
                        labels[pixIdx] = existingLabel;
                        std::vector<point> belong2Star = _floodingPlateau(pixIdx, labels, img, threshold, existingLabel);
                        if(is3D)
                        {
                            for(std::vector<point>::iterator p=belong2Star.begin(); p!=belong2Star.end(); ++p)
                                p->SetZ(static_cast<double>(iLayer));
                        }

                        {
                            std::lock_guard<std::mutex> lock(flabelMutex); 
                            fdetectedStars[existingLabel].addPixel(belong2Star);
                        }
                    }
                }
            }

            std::unique_lock<std::mutex> lock(flabelMutex);
            if(flabel == nullptr)
            {
                size_t iS = 0;
                if(!is3D)
                {
                    flabel = std::shared_ptr< FITScube >(new FITSimg<int64_t>(2,{img->Size(1),img->Size(2)}));
                    iS=0;
                }
                else
                {
                    flabel = std::shared_ptr< FITScube >(new FITSimg<int64_t>(3,{img->Size(1),img->Size(2),fsky->Size(3)}));
                    iS=iLayer*img->Size(1)*img->Size(2);
                }
                
                flabel->SetName("label_image");
                *(flabel->GetData<int64_t>()) *= 0;
                for(size_t i=0; i<labels.size(); ++i)
                {
                    if(labels[i] < 1)
                        continue;

                    (*(flabel->GetData<int64_t>()))[i+iS]=static_cast<int64_t>(labels[i]);
                }
            }
            else
            {
                size_t iS = (!is3D) ? 0 : iLayer*img->Size(1)*img->Size(2);
                for(size_t i=0; i<labels.size(); ++i)
                {
                    if(labels[i] < 1)
                        continue;
                    (*(flabel->GetData<int64_t>()))[i+iS]=static_cast<int64_t>(labels[i]);
                }
            }
            lock.unlock(); // explicit unlock
        }


        /**
         * @brief Helper function to get the 8 neighbor pixel indices for a given pixel index
         * @param pixIdx starting pixel index
         * @return \c std::vector<size_t> vector of neighbor pixel indices
         */
        std::vector<size_t> starFinder::_getNeighborIndices(const size_t& pixIdx, const size_t& w, const size_t& h) const
        {
            std::vector<size_t> neighbors;
            neighbors.reserve(8);

            size_t width = w;   // size along x-axis
            size_t height = h;  // size along y-axis

            // Calculate x,y from linear index: pixIdx = x + y*width
            std::vector<size_t> coo = fsky->PixelCoordinates(pixIdx);

            // Check all 8 neighbors using linear index arithmetic
            // Top row (y-1)
            if(coo[1] > 0)
            {
                if(coo[0] > 0)          neighbors.push_back(pixIdx - width - 1);  // bottom-left
                                        neighbors.push_back(pixIdx - width);      // bottom
                if((coo[0]+1) < width)  neighbors.push_back(pixIdx - width + 1);  // bottom-right
            }

            // Same row (y)
            if(coo[0] > 0)          neighbors.push_back(pixIdx - 1);          // left
            if((coo[0]+1) < width)  neighbors.push_back(pixIdx + 1);          // right

            // Bottom row (y+1)
            if((coo[1]+1) < height)
            {
                if(coo[0] > 0)          neighbors.push_back(pixIdx + width - 1);  // top-left
                                        neighbors.push_back(pixIdx + width);      // top
                if((coo[0]+1) < width)  neighbors.push_back(pixIdx + width + 1);  // top-right
            }

            return neighbors;
        }

        /**
         * @brief Helper function to flood fill a plateau region with the same label
         * @param pixIdx Starting pixel index
         * @param labels Label array
         * @param data Image data array
         * @param threshold Current threshold value
         * @param label Label to assign to plateau region
         */
        std::vector<point> starFinder::_floodingPlateau(const size_t& pixIdx, std::valarray<size_t>& labels, 
                                       const std::shared_ptr<FITScube> img, 
                                       const double& threshold,
                                       const size_t& label)
        {
            size_t width = img->Size(1);
            size_t height = img->Size(2);

            std::vector<point> plateauPixels;
            std::queue<size_t> toProcess;
            toProcess.push(pixIdx);
            labels[pixIdx] = label;
            if(!is3D)
                plateauPixels.push_back(point::pointFromVector<size_t>(img->PixelCoordinates(pixIdx)));
            else
            {
                std::vector<size_t> coords = img->PixelCoordinates(pixIdx);
                coords.push_back(0); // Z=0 for 2D image in 3D context
                plateauPixels.push_back(point::pointFromVector<size_t>(coords));
            }

            while(!toProcess.empty())
            {
                size_t currentIdx = toProcess.front();
                toProcess.pop();

                std::vector<size_t> neighbors = _getNeighborIndices(currentIdx, width, height);
                for(const size_t& nIdx : neighbors)
                {
                    // Check if neighbor is at same level (plateau) and unlabeled
                    if(!img->Masked(nIdx) && 
                       labels[nIdx] == 0 && 
                       std::abs(img->DoubleValueAtPixel(nIdx) - threshold) < 1e-10)
                    {
                        labels[nIdx] = label;
                        toProcess.push(nIdx);
                        if(!is3D)
                            plateauPixels.push_back(point::pointFromVector<size_t>(img->PixelCoordinates(nIdx)));
                        else
                        {
                            std::vector<size_t> coords = img->PixelCoordinates(nIdx);
                            coords.push_back(0); // Z=0 for 2D image in 3D context
                            plateauPixels.push_back(point::pointFromVector<size_t>(coords));
                        }
                    }
                }
            }

            return plateauPixels;
        }
#pragma endregion
#pragma region -- stellar properties evaluation
        /**
         * @brief Evaluate stellar properties using the specified method
         * @param star Star object to evaluate
         * @param method Evaluation method
         */
        void starFinder::EvalStellarProperties(Stars& star, const enum Stars::evalType& method)
        {
            switch (method)
            {
                case Stars::evalType::MOMENT:
                    _evalStellarPropertiesMoment(star);
                    break;
                default:
                    _evalStellarPropertiesModel(star,method);
            }

            _evalApperturePhotometry(star);
        }

        void starFinder::EvalStellarProperties(const size_t& starID, const enum Stars::evalType& method)
        {
            std::map<size_t, Stars>::iterator it = fdetectedStars.find(starID);
            if(it == fdetectedStars.end())
                throw std::invalid_argument("DST::starFinder::EvalStellarProperties: Star ID not found.["+std::to_string(__LINE__)+"]");

            EvalStellarProperties(it->second, method);
        }

        /**
         * @brief Evaluate stellar properties using moment method
         * @param star Star object to evaluate
         */
        void starFinder::_evalStellarPropertiesMoment(Stars& star)
        {
            if(fsky == nullptr || star.FilledPixels().empty())
                return;

            double fluxSum   = 0.0;
            double xWeighted = 0.0;
            double yWeighted = 0.0;
            double maxima    = -std::numeric_limits<double>::infinity();
            size_t maxIndex  = 0;
            double iL        = 0.0;

            // Loop only over star pixels
            for(std::vector<point>::const_iterator pix = star.FilledPixels().cbegin(); pix!=star.FilledPixels().cend(); ++pix)
            {
                size_t idx = fsky->PixelIndex(pix->GetPoint());
                
                if(fsky->Masked(idx))
                    continue;

                if(is3D && iL == 0.0)
                    iL = pix->Z();

                double val = fsky->DoubleValueAtPixel(idx); // type-agnostic access
                fluxSum   += val;
                xWeighted += pix->X() * val;
                yWeighted += pix->Y() * val;

                if(val > maxima)
                {
                    maxima   = val;
                    maxIndex = idx;
                }
            }

            if(fluxSum > 0.0)
            {
                double cx = xWeighted / fluxSum;
                double cy = yWeighted / fluxSum;
                if(is3D)
                    star.setPosition(point({cx, cy,iL}));
                else
                    star.setPosition(point({cx, cy}));

                // Second-order intensity-weighted moments → FWHM estimate
                double M20 = 0.0, M02 = 0.0;
                for(std::vector<point>::const_iterator pix = star.FilledPixels().cbegin(); pix!=star.FilledPixels().cend(); ++pix)
                {
                    size_t idx = fsky->PixelIndex(pix->GetPoint());
                    if(fsky->Masked(idx))
                        continue;
                    double val = fsky->DoubleValueAtPixel(idx);
                    double dx  = pix->X() - cx;
                    double dy  = pix->Y() - cy;
                    M20 += dx * dx * val;
                    M02 += dy * dy * val;
                }
                M20 /= fluxSum;
                M02 /= fluxSum;
                constexpr double k = 2.354820045; // 2*sqrt(2*ln2)
                star.setFWHM(point({k * std::sqrt(M20), k * std::sqrt(M02)}));
            }

            if(std::isfinite(maxima))
            {
                star.setMaxima(maxima);
                star.setMaximaPosition(point::pointFromVector<size_t>(fsky->PixelCoordinates(maxIndex)));
            }
        }

        /**
         * @brief Evaluate stellar properties using Gaussian Fit of the star profile
         * @param star Star object to evaluate
         * @param integrated Whether to use integrated Gaussian fit within pixels surface or to use the value at pixel centers (default)
         */
        void starFinder::_evalStellarPropertiesModel(Stars& star, const DST::Stars::evalType& how)
        {
            if(fsky == nullptr || star.FilledPixels().empty())
                return;

            //double fluxSum   = 0.0;
            //double xWeighted = 0.0;
            //double yWeighted = 0.0;
            double iL        = 0.0;

            //1- Crop a small region arround the star to speed up the fitting
            //1.a- Find the centroid position relying on the moment methods.
            _evalStellarPropertiesMoment(star);

            //1.b- Define a small bounding box the accapsulata all pixel which belong to that star + e small margine etending by all the width on each side of the bounding box.
            
            double xMin = std::numeric_limits<double>::infinity();
            double xMax = -std::numeric_limits<double>::infinity();
            double yMin = std::numeric_limits<double>::infinity();
            double yMax = -std::numeric_limits<double>::infinity();

            for(std::vector<point>::const_iterator pix = star.FilledPixels().cbegin(); pix!=star.FilledPixels().cend(); ++pix)
            {
                if(pix->X() < xMin) xMin = pix->X();
                if(pix->X() > xMax) xMax = pix->X();
                if(pix->Y() < yMin) yMin = pix->Y();
                if(pix->Y() > yMax) yMax = pix->Y();

                if(is3D && iL == 0.0)
                    iL = pix->Z();

                size_t idx = fsky->PixelIndex(pix->GetPoint());
                if(fsky->Masked(idx))
                    continue;

                //double val = fsky->DoubleValueAtPixel(idx);
            }

            xMin = std::max(0.0, xMin -(xMax - xMin)/2.);
            xMax = std::min(static_cast<double>(fsky->Size(1)), xMax +(xMax - xMin)/2.);
            yMin = std::max(0.0, yMin -(yMax - yMin)/2.);
            yMax = std::min(static_cast<double>(fsky->Size(2)), yMax +(yMax - yMin)/2.);

            //1.c- crop the 2D image to thise bounding box.
            std::shared_ptr< FITSimg<double> > cpy = nullptr;
            std::shared_ptr<FITScube> croppedImage = nullptr;
            std::shared_ptr<FITScube> croppedLabel = nullptr;
            if(is3D)
            {
                croppedImage = fsky->Window(xMin, yMin, (xMax-xMin), (yMax-yMin))->Layer(static_cast<size_t>(iL+0.5));
                croppedLabel = flabel->Window(xMin, yMin, (xMax-xMin), (yMax-yMin))->Layer(static_cast<size_t>(iL+0.5));
            }
            else
            {
                croppedImage = fsky->Window(xMin, yMin, (xMax-xMin), (yMax-yMin));
                croppedLabel = flabel->Window(xMin, yMin, (xMax-xMin), (yMax-yMin));
            }

            for(size_t i=0; i<croppedImage->Nelements(); ++i)
            {
                if(croppedLabel->Int64ValueAtPixel(i) != static_cast<int64_t>(star.ID()))
                    croppedImage->MaskPixel(i);
            }
            croppedLabel.reset();

            try
            {
                cpy = std::dynamic_pointer_cast<FITSimg<double>>(croppedImage->cast<double>());
            }
            catch(FITSexception& e)
            {
                if(e.errorCode() != SHARED_BADARG)  // or whatever identifies "already same type"
                    throw;
                cpy = std::make_shared< FITSimg<double> >(FITSimg<double>(*std::dynamic_pointer_cast<FITSimg<double>>(croppedImage)));
            }
            croppedImage.reset();


            //2- Fits the image with the corresponding model (Gaussian, Moffat, etc..) and extract the stellar properties from the fitted model parameters.
            star.setEvalMethod(how);
            starFitter fitter(cpy, how);


            ROOT::Minuit2::MnUserParameters upar;
            upar.Add("flux", star.maxima(), 0.1);
            upar.Add("x0", star.position().X()-xMin, 0.01);
            upar.Add("y0", star.position().Y()-yMin, 0.01);

            switch (how)
            {
                case Stars::evalType::GAUSS:
                case Stars::evalType::iGAUSS:
                    upar.Add("sigma", 1.0, 0.01);
                    break;

                case Stars::evalType::GAUSS2D:
                case Stars::evalType::iGAUSS2D:
                    upar.Add("sigma_x", 1.0, 0.01);
                    upar.Add("sigma_y", 1.0, 0.01);
                    upar.Add("theta"  , 0.0, 0.001);
                    break;
                
                case Stars::evalType::MOFFAT:
                    upar.Add("alpha", 1.5, 0.01);
                    upar.Add("beta",  1.5, 0.01);
                    break;
                
                case Stars::evalType::MOFFAT2D:
                    upar.Add("alpha", 1.5, 0.001);
                    upar.Add("beta" , 1.5, 0.001);
                    upar.Add("gamma", 0.8, 0.001);
                    upar.Add("theta", 0.0, 0.01);
                    break;

                case Stars::evalType::LORENTZIAN:
                    upar.Add("gamma", 1.0, 0.01, 0.01, 50.0);
                    break;

                case Stars::evalType::LORENTZIAN2D:
                    upar.Add("gamma", 1.0, 0.01, 0.01, 50.0);
                    upar.Add("q",     1.0, 0.01, 0.01,  2.0);
                    upar.Add("theta", 0.0, 0.01);
                    break;

                case Stars::evalType::VOIGT:
                    upar.Add("fG", 1.0, 0.01, 0.01, 50.0);
                    upar.Add("fL", 1.0, 0.01, 0.01, 50.0);
                    break;

                case Stars::evalType::VOIGT2D:
                    upar.Add("fG",    1.0, 0.01, 0.01, 50.0);
                    upar.Add("fL",    1.0, 0.01, 0.01, 50.0);
                    upar.Add("q",     1.0, 0.01, 0.01,  2.0);
                    upar.Add("theta", 0.0, 0.01);
                    break;

                case Stars::evalType::SERSIC:
                    upar.Add("re", 1.0, 0.01, 0.01, 50.0);
                    upar.Add("n",  1.0, 0.01, 0.01, 10.0);
                    break;

                case Stars::evalType::SERSIC2D:
                    upar.Add("re",    1.0, 0.01, 0.01, 50.0);
                    upar.Add("n",     1.0, 0.01, 0.01, 10.0);
                    upar.Add("q",     1.0, 0.01, 0.01,  2.0);
                    upar.Add("theta", 0.0, 0.01);
                    break;

                case Stars::evalType::KING:
                    upar.Add("rc", 1.0, 0.01, 0.01, 50.0);
                    upar.Add("rt", 3.0, 0.01, 0.01, 50.0);
                    break;

                case Stars::evalType::KING2D:
                    upar.Add("rc",    1.0, 0.01, 0.01, 50.0);
                    upar.Add("rt",    3.0, 0.01, 0.01, 50.0);
                    upar.Add("q",     1.0, 0.01, 0.01,  2.0);
                    upar.Add("theta", 0.0, 0.01);
                    break;

                case Stars::evalType::AIRY:
                    upar.Add("a", 1.0, 0.01, 0.01, 50.0);
                    break;

                default:
                    break;
            }

            upar.Add("background", 0.0, 0.01);

            ROOT::Minuit2::MnMigrad migrad(fitter, upar);
            //migrad.SetPrintLevel(-1); // Suppress output
            ROOT::Minuit2::FunctionMinimum min = migrad();

            if(!min.IsValid())
                std::cerr << "DST::starFinder::_evalStellarPropertiesModel: Fit did not converge for star ID "<<star.ID()<<"."<<std::endl;

            std::vector<double> best_fit_params;
            best_fit_params.push_back(min.UserState().Value("flux"));
            best_fit_params.push_back(min.UserState().Value("x0") + xMin);
            best_fit_params.push_back(min.UserState().Value("y0") + yMin);

            switch (how)
            {
                case Stars::evalType::GAUSS:
                case Stars::evalType::iGAUSS:
                    best_fit_params.push_back(min.UserState().Value("sigma"));
                    star.setFWHM(point({2.354820045 * min.UserState().Value("sigma"),
                                        2.354820045 * min.UserState().Value("sigma")})); // FWHM = 2.355 * sigma for Gaussian
                    star.setAngle(0.0);
                    break;

                case Stars::evalType::GAUSS2D:
                case Stars::evalType::iGAUSS2D:
                    best_fit_params.push_back(min.UserState().Value("sigma_x"));
                    best_fit_params.push_back(min.UserState().Value("sigma_y"));
                    best_fit_params.push_back(min.UserState().Value("theta"  ));
                    star.setFWHM(point({2.354820045 * min.UserState().Value("sigma_x"),
                                        2.354820045 * min.UserState().Value("sigma_y")}));
                    star.setAngle(min.UserState().Value("theta"));
                    break;
                
                case Stars::evalType::MOFFAT:
                    best_fit_params.push_back(min.UserState().Value("alpha"));
                    best_fit_params.push_back(min.UserState().Value("beta" ));
                    star.setFWHM(point({2.0 * min.UserState().Value("alpha") * std::sqrt(std::pow(2.0, 1.0/min.UserState().Value("beta")) - 1.0),
                                        2.0 * min.UserState().Value("alpha") * std::sqrt(std::pow(2.0, 1.0/min.UserState().Value("beta")) - 1.0)})); // FWHM = 2*alpha*sqrt(2^(1/beta)-1) for Moffat
                    star.setAngle(0.0);
                    break;
                
                case Stars::evalType::MOFFAT2D:
                    best_fit_params.push_back(min.UserState().Value("alpha"));
                    best_fit_params.push_back(min.UserState().Value("beta" ));
                    best_fit_params.push_back(min.UserState().Value("gamma"));
                    best_fit_params.push_back(min.UserState().Value("theta"));
                    star.setFWHM(point({2.0 * min.UserState().Value("alpha") * std::sqrt(std::pow(2.0, 1.0/min.UserState().Value("beta")) - 1.0),
                                        2.0 * min.UserState().Value("alpha") * std::sqrt(std::pow(2.0, 1.0/min.UserState().Value("beta")) - 1.0)})); // FWHM = 2*alpha*sqrt(2^(1/beta)-1) for Moffat
                    star.setAngle(min.UserState().Value("theta"));
                    break;

                case Stars::evalType::LORENTZIAN:
                    best_fit_params.push_back(min.UserState().Value("gamma"));
                    star.setFWHM(point({2.0*min.UserState().Value("gamma"),
                                        2.0*min.UserState().Value("gamma")}));
                    star.setAngle(0.0);
                    break;

                case Stars::evalType::LORENTZIAN2D:
                    best_fit_params.push_back(min.UserState().Value("gamma"));
                    best_fit_params.push_back(min.UserState().Value("q"));
                    best_fit_params.push_back(min.UserState().Value("theta"));
                    star.setFWHM(point({2.0*min.UserState().Value("gamma"),
                                        2.0*min.UserState().Value("gamma")*min.UserState().Value("q")}));
                    star.setAngle(min.UserState().Value("theta"));
                    break;

                case Stars::evalType::VOIGT:
                {
                    double fG = min.UserState().Value("fG");
                    double fL = min.UserState().Value("fL");
                    double fwhm = 0.5346*fL + std::sqrt(0.2166*fL*fL + fG*fG);
                    best_fit_params.push_back(fG);
                    best_fit_params.push_back(fL);
                    star.setFWHM(point({fwhm, fwhm}));
                    star.setAngle(0.0);
                    break;
                }

                case Stars::evalType::VOIGT2D:
                {
                    double fG = min.UserState().Value("fG");
                    double fL = min.UserState().Value("fL");
                    double q  = min.UserState().Value("q");
                    double fwhm = 0.5346*fL + std::sqrt(0.2166*fL*fL + fG*fG);
                    best_fit_params.push_back(fG);
                    best_fit_params.push_back(fL);
                    best_fit_params.push_back(q);
                    best_fit_params.push_back(min.UserState().Value("theta"));
                    star.setFWHM(point({fwhm, fwhm*q}));
                    star.setAngle(min.UserState().Value("theta"));
                    break;
                }

                case Stars::evalType::SERSIC:
                    best_fit_params.push_back(min.UserState().Value("re"));
                    best_fit_params.push_back(min.UserState().Value("n"));
                    star.setFWHM(point({2.0*min.UserState().Value("re"),
                                        2.0*min.UserState().Value("re")}));
                    star.setAngle(0.0);
                    break;

                case Stars::evalType::SERSIC2D:
                    best_fit_params.push_back(min.UserState().Value("re"));
                    best_fit_params.push_back(min.UserState().Value("n"));
                    best_fit_params.push_back(min.UserState().Value("q"));
                    best_fit_params.push_back(min.UserState().Value("theta"));
                    star.setFWHM(point({2.0*min.UserState().Value("re"),
                                        2.0*min.UserState().Value("re")*min.UserState().Value("q")}));
                    star.setAngle(min.UserState().Value("theta"));
                    break;

                case Stars::evalType::KING:
                    best_fit_params.push_back(min.UserState().Value("rc"));
                    best_fit_params.push_back(min.UserState().Value("rt"));
                    star.setFWHM(point({2.0*min.UserState().Value("rc"),
                                        2.0*min.UserState().Value("rc")}));
                    star.setAngle(0.0);
                    break;

                case Stars::evalType::KING2D:
                    best_fit_params.push_back(min.UserState().Value("rc"));
                    best_fit_params.push_back(min.UserState().Value("rt"));
                    best_fit_params.push_back(min.UserState().Value("q"));
                    best_fit_params.push_back(min.UserState().Value("theta"));
                    star.setFWHM(point({2.0*min.UserState().Value("rc"),
                                        2.0*min.UserState().Value("rc")*min.UserState().Value("q")}));
                    star.setAngle(min.UserState().Value("theta"));
                    break;

                case Stars::evalType::AIRY:
                    best_fit_params.push_back(min.UserState().Value("a"));
                    star.setFWHM(point({2.0*min.UserState().Value("a"),
                                        2.0*min.UserState().Value("a")}));
                    star.setAngle(0.0);
                    break;

                default:
                    break;
            }
            best_fit_params.push_back(min.UserState().Value("background"));
            star.setPSFparameters(best_fit_params);


            cpy.reset();
        }

#pragma region -- aperture photometry

        /**
         * @brief Measure total flux, background and SNR via aperture + sky-annulus photometry.
         * @details Aperture radius = 1.5 × FWHM_max; annulus from 2.0 to 3.5 × FWHM_max.
         *          Background per pixel = median of annulus pixels (robust against neighbours).
         *          σ_sky = 1.4826 × MAD (normalised MAD, equivalent Gaussian σ).
         *          SNR computed with the CCD equation assuming gain = 1 ADU/e⁻.
         *          If flabel is available, pixels belonging to other stars are excluded from
         *          both the aperture and the annulus.
         * @param star Star object whose flux, background and SNR are set.
         */
        void starFinder::_evalApperturePhotometry(Stars& star)
        {
            if(fsky == nullptr || star.FilledPixels().empty())
                return;

            const double cx = star.position().X();
            const double cy = star.position().Y();

            // -- Determine FWHM_max ---------------------------------------------------
            double fwhm_max = 0.0;
            if(star.FWHM().size() >= 2)
                fwhm_max = std::max(star.FWHM().X(), star.FWHM().Y());
            else if(star.FWHM().size() == 1)
                fwhm_max = star.FWHM().X();

            if(fwhm_max <= 0.0)
            {
                // Safety fallback: use bounding-box diagonal of FilledPixels / 2
                double xMin =  std::numeric_limits<double>::infinity();
                double xMax = -std::numeric_limits<double>::infinity();
                double yMin =  std::numeric_limits<double>::infinity();
                double yMax = -std::numeric_limits<double>::infinity();
                for(const auto& p : star.FilledPixels())
                {
                    if(p.X() < xMin) xMin = p.X();
                    if(p.X() > xMax) xMax = p.X();
                    if(p.Y() < yMin) yMin = p.Y();
                    if(p.Y() > yMax) yMax = p.Y();
                }
                fwhm_max = 0.5 * std::sqrt((xMax-xMin)*(xMax-xMin) + (yMax-yMin)*(yMax-yMin));
                if(fwhm_max <= 0.0) return;
            }

            const double r_ap  = 1.5 * fwhm_max;
            const double r_in  = 2.0 * fwhm_max;
            const double r_out = 3.5 * fwhm_max;

            // -- Layer for 3D images --------------------------------------------------
            double iL = 0.0;
            if(is3D)
                iL = star.FilledPixels().front().Z();

            // -- Bounding box clamped to image bounds ---------------------------------
            const size_t x0 = static_cast<size_t>(std::max(0.0, cx - r_out - 1.0));
            const size_t x1 = static_cast<size_t>(std::min(static_cast<double>(fsky->Size(1) - 1), cx + r_out + 1.0));
            const size_t y0 = static_cast<size_t>(std::max(0.0, cy - r_out - 1.0));
            const size_t y1 = static_cast<size_t>(std::min(static_cast<double>(fsky->Size(2) - 1), cy + r_out + 1.0));

            std::vector<double> aperture_vals;
            std::vector<double> annulus_vals;

            for(size_t ix = x0; ix <= x1; ++ix)
            {
                for(size_t iy = y0; iy <= y1; ++iy)
                {
                    const size_t idx = is3D
                        ? fsky->PixelIndex(std::vector<double>({static_cast<double>(ix), static_cast<double>(iy), iL}))
                        : fsky->PixelIndex(std::vector<double>({static_cast<double>(ix), static_cast<double>(iy)}));

                    if(fsky->Masked(idx))
                        continue;

                    // Exclude pixels belonging to other detected stars
                    if(flabel != nullptr)
                    {
                        const size_t lbl = static_cast<size_t>(flabel->DoubleValueAtPixel(idx));
                        if(lbl != 0 && lbl != star.ID())
                            continue;
                    }

                    const double dx = static_cast<double>(ix) - cx;
                    const double dy = static_cast<double>(iy) - cy;
                    const double r  = std::sqrt(dx*dx + dy*dy);
                    const double val = fsky->DoubleValueAtPixel(idx);

                    if(r <= r_ap)
                        aperture_vals.push_back(val);
                    else if(r >= r_in && r < r_out)
                        annulus_vals.push_back(val);
                }
            }

            if(annulus_vals.empty()) return;

            // -- Background: median of annulus ----------------------------------------
            std::sort(annulus_vals.begin(), annulus_vals.end());
            const size_t n_ann = annulus_vals.size();
            const double background = (n_ann % 2 == 0)
                ? 0.5 * (annulus_vals[n_ann/2 - 1] + annulus_vals[n_ann/2])
                : annulus_vals[n_ann/2];

            // -- σ_sky: normalised MAD -------------------------------------------------
            std::vector<double> deviations;
            deviations.reserve(n_ann);
            for(double v : annulus_vals)
                deviations.push_back(std::abs(v - background));
            std::sort(deviations.begin(), deviations.end());
            const double sigma_sky = 1.4826 * ((n_ann % 2 == 0)
                ? 0.5 * (deviations[n_ann/2 - 1] + deviations[n_ann/2])
                : deviations[n_ann/2]);

            // -- Aperture flux (background-subtracted sum) ----------------------------
            double flux = 0.0;
            for(double v : aperture_vals)
                flux += v - background;
            const size_t n_ap = aperture_vals.size();

            // -- SNR: CCD equation (gain = 1 ADU/e⁻ assumed) -------------------------
            const double snr = (n_ap > 0 && sigma_sky > 0.0)
                ? flux / std::sqrt(std::abs(flux) + static_cast<double>(n_ap) * sigma_sky * sigma_sky)
                : 0.0;

            star.setBackground(background);
            star.setTotalFlux(flux);
            star.setSNR(snr);
        }

#pragma endregion

#pragma region -- ploting

#ifdef __WITH_QT__

        /**
         * @brief Plot detected stars on the input image
         * @param starColor Color of the star markers
         * @param starSize Size of the star markers
         */
        void starFinder::Plot(QCPAxisRect* axisRect, size_t iL, const QPen& starPen) const
        {
            if(fsky == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::show","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            // Ligne MANQUANTE :
            if(flabel == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","starFinder::plot","Label image is null.["+std::to_string(__LINE__)+"]");

            if(fsky->GetDimension() == 3 && iL > fsky->Size(3))
                throw std::out_of_range("DST::starFinder::show: Layer index out of range.["+std::to_string(__LINE__)+"]");

            QCustomPlot* qPlot = axisRect->parentPlot();
            QCPAxis* xAxis = axisRect->axis(QCPAxis::atBottom);
            QCPAxis* yAxis = axisRect->axis(QCPAxis::atLeft);

            QMetaObject::invokeMethod(qPlot, [this,iL,starPen,qPlot,xAxis,yAxis,axisRect]() mutable
            {
                {
                    // Start plotting the sky image for the corresponding layer

                    QCPGraph* point = new QCPGraph(xAxis, yAxis);
                    point->setLineStyle(QCPGraph::lsNone);
                    QVector<double> x;
                    QVector<double> y;

                    QCPScatterStyle StarPos(QCPScatterStyle::ssCross, starPen.color(), starPen.color(), 10);

                    struct StarInfo {
                        size_t id;
                        double px, py;
                        bool hasWCS;
                        std::vector<double> wcs;
                        std::vector<std::string> wcsLabels;
                        bool hasFWHM;
                        double fwhmX, fwhmY;
                        bool hasAngle;
                        double angle;
                        bool hasSNR;
                        double snr;
                    };
                    std::vector<StarInfo> tooltipList;
                    bool skyHasWCS = fsky->getNumberOfWCS() > 0;
                    std::vector<std::string> wcsAxisLabels;
                    if (skyHasWCS)
                        wcsAxisLabels = fsky->getWCS().CTYPE(0);

                    for(auto iS: fdetectedStars)
                    {
                        unsigned int n_curve_point = 500;
                        QCPCurve* cell = nullptr;

                        DST::Math::point center = iS.second.position();
                        bool skipe = false;

                        if(center.size() >=3)
                        {
                            if(static_cast<size_t>(center.Z()+1.5) != iL)
                                skipe |= true;
                        }

                        if(skipe)
                            continue;

                        x.push_back(center.X());
                        y.push_back(center.Y());

                        StarInfo info;
                        info.id        = iS.second.ID();
                        info.px        = center.X();
                        info.py        = center.Y();
                        info.hasWCS    = skyHasWCS;
                        if (skyHasWCS)
                        {
                            info.wcs       = fsky->WorldCoordinates(iS.second.position().GetPoint());
                            info.wcsLabels = wcsAxisLabels;
                        }
                        info.hasFWHM   = iS.second.FWHM().size() >= 2;
                        info.fwhmX     = info.hasFWHM ? iS.second.FWHM().X() : 0.0;
                        info.fwhmY     = info.hasFWHM ? iS.second.FWHM().Y() : 0.0;
                        info.hasAngle  = !std::isnan(iS.second.angle());
                        info.angle     = info.hasAngle ? iS.second.angle() : 0.0;
                        info.hasSNR    = !std::isnan(iS.second.SNR());
                        info.snr       = info.hasSNR ? iS.second.SNR() : 0.0;
                        tooltipList.push_back(info);

                        if(iS.second.FWHM().size() < 1)
                            continue;

                        cell = new QCPCurve(xAxis, yAxis);

                        QVector<double> x_ell(n_curve_point);
                        QVector<double> y_ell(n_curve_point);

                        double a = iS.second.FWHM().X()/2.;
                        double b = iS.second.FWHM().Y()/2.;

                        double phi = std::isnan(iS.second.angle()) ? 0.0 : iS.second.angle()*acos(-1)/180.;

                        for(unsigned int n = 0; n < n_curve_point; n++)
                        {
                            double t = 2.*acos(-1)*static_cast<double>(n)/static_cast<double>(n_curve_point-1);
                            DST::Math::point pt({a*cos(t)*cos(phi) - b*sin(t)*sin(phi), b*sin(t)*cos(phi)+a*cos(t)*sin(phi)});
                            pt += center;

                            x_ell[n] = pt.X();
                            y_ell[n] = pt.Y();
                        }
                        cell->setData(x_ell, y_ell);
                        cell->setPen(starPen);

                        // Aperture-photometry radii: 1.5, 2.0, 3.5 × FWHM_max (dashed)
                        double fwhm_max = std::max(iS.second.FWHM().X(), iS.second.FWHM().Y());
                        if(fwhm_max > 0.0)
                        {
                            QPen dashedPen(starPen);
                            dashedPen.setStyle(Qt::DashLine);
                            for(double factor : {1.5, 2.0, 3.5})
                            {
                                QCPCurve* ring = new QCPCurve(xAxis, yAxis);
                                double r = factor * fwhm_max;
                                QVector<double> xr(n_curve_point), yr(n_curve_point);
                                for(unsigned int n = 0; n < n_curve_point; n++)
                                {
                                    double t = 2.*acos(-1)*static_cast<double>(n)/static_cast<double>(n_curve_point-1);
                                    xr[n] = center.X() + r * cos(t);
                                    yr[n] = center.Y() + r * sin(t);
                                }
                                ring->setData(xr, yr);
                                ring->setPen(dashedPen);
                            }
                        }
                    }

                    point->setData(x, y);
                    point->setScatterStyle(StarPos);

                    // Hover tooltip — replace any previous handler for this axis rect
                    QObject* oldHandler = axisRect->findChild<QObject*>("_starTooltipHandler_");
                    if (oldHandler) delete oldHandler;
                    QObject* handler = new QObject(axisRect);
                    handler->setObjectName("_starTooltipHandler_");
                    QObject::connect(qPlot, &QCustomPlot::mouseMove, handler,
                        [tooltipList, qPlot, xAxis, yAxis, axisRect](QMouseEvent* event)
                        {
                            // Ignore events outside this axis rect
                            if (!axisRect->rect().contains(event->pos()))
                            {
                                QToolTip::hideText();
                                return;
                            }

                            double mx = xAxis->pixelToCoord(event->pos().x());
                            double my = yAxis->pixelToCoord(event->pos().y());

                            double minDist = std::numeric_limits<double>::max();
                            const StarInfo* nearest = nullptr;
                            for (const auto& info : tooltipList)
                            {
                                double dx = info.px - mx;
                                double dy = info.py - my;
                                double dist = dx*dx + dy*dy;
                                if (dist < minDist) { minDist = dist; nearest = &info; }
                            }

                            if (!nearest) { QToolTip::hideText(); return; }

                            // Only show within 20 px of the marker
                            double xPix = xAxis->coordToPixel(nearest->px);
                            double yPix = yAxis->coordToPixel(nearest->py);
                            double dxPix = xPix - event->pos().x();
                            double dyPix = yPix - event->pos().y();
                            if (dxPix*dxPix + dyPix*dyPix > 400.0) { QToolTip::hideText(); return; }

                            QString text = QString("Star #%1").arg(nearest->id);
                            text += QString("\nPixel: (%1, %2)").arg(nearest->px, 0, 'f', 1).arg(nearest->py, 0, 'f', 1);
                            if (nearest->hasWCS && nearest->wcs.size() >= 2)
                                for (size_t i = 0; i < nearest->wcs.size() && i < nearest->wcsLabels.size(); ++i)
                                    text += QString("\n%1: %2").arg(QString::fromStdString(nearest->wcsLabels[i])).arg(nearest->wcs[i], 0, 'f', 6);
                            if (nearest->hasFWHM)
                                text += QString("\nFWHM: %1 x %2 px").arg(nearest->fwhmX, 0, 'f', 2).arg(nearest->fwhmY, 0, 'f', 2);
                            if (nearest->hasAngle)
                                text += QString("\nAngle: %1 deg").arg(nearest->angle, 0, 'f', 1);
                            if (nearest->hasSNR)
                                text += QString("\nSNR: %1").arg(nearest->snr, 0, 'f', 2);

                            QToolTip::showText(QCursor::pos(), text, qPlot);
                        }
                    );

                    qPlot->replot();
                }

                {
                        std::lock_guard<std::mutex> dk(QtFITS::msg_m);
                        if (!QtFITS::finished)
                             QtFITS::finished |= true;
                        QtFITS::msg_cv.notify_one(); // notify while holding msg_m (or move notify outside if preferred)
                        QtFITS::finished &= false;
                } // msg_m released here

            });
            
        }

#endif

#pragma endregion
    }
#pragma endregion

} // namespace name
#pragma endregion
#pragma endregion

// ---------------------------------------------------------------------------
// starFitter plotting methods (Qt-only)
// ---------------------------------------------------------------------------

#ifdef __WITH_QT__

namespace DST { namespace algorithms {

void starFitter::Plot(QCPAxisRect* axisRect,
                      const std::vector<double>& params,
                      size_t oversample,
                      const std::string& saveFile) const
{
    if(!fsky)
        throw FITSexception(BAD_DIMEN, "DST", "starFitter::Plot",
                            "fsky is null. [" + std::to_string(__LINE__) + "]");

    oversample = std::clamp(oversample, size_t(1), size_t(100));

    const size_t nx  = fsky->Size(1);
    const size_t ny  = fsky->Size(2);
    const size_t onx = nx * oversample;
    const size_t ony = ny * oversample;

    // Generate model on oversampled grid outside Qt thread (can be expensive)
    // Zmodel[row][col], fractional fsky pixel coords: x = col/oversample, y = row/oversample
    const double inv_s = 1.0 / static_cast<double>(oversample);
    std::vector<std::vector<double>> Zmodel(ony, std::vector<double>(onx));
    for(size_t iy = 0; iy < ony; iy++)
        for(size_t ix = 0; ix < onx; ix++)
            Zmodel[iy][ix] = Eval(static_cast<double>(ix) * inv_s,
                                  static_cast<double>(iy) * inv_s, params);

    // Compute contour geometry outside Qt thread
    // QCustomContours uses row/col indices as coordinates → must scale by inv_s when rendering
    QCustomContours contours(Zmodel, 8);
    contours.setPen(QPen(Qt::white, 1.0));
    contours.setLabelFontSize(0);

    QCustomPlot* qPlot = axisRect->parentPlot();
    QCPAxis*     xAxis = axisRect->axis(QCPAxis::atBottom);
    QCPAxis*     yAxis = axisRect->axis(QCPAxis::atLeft);

    QMetaObject::invokeMethod(qPlot,
        [=, this]() mutable
        {
            // ---- 1. Heatmap of fsky with Viridis colormap ----
            QCPColorMap* colorMap = new QCPColorMap(xAxis, yAxis);
            colorMap->data()->setSize(static_cast<int>(nx), static_cast<int>(ny));
            colorMap->data()->setRange(
                QCPRange(0.0, static_cast<double>(nx) - 1.0),
                QCPRange(0.0, static_cast<double>(ny) - 1.0));
            colorMap->data()->fill(0.0);

            for(size_t i = 0; i < fsky->Size(); i++) {
                if(fsky->isMasked(i)) continue;
                const auto pix = fsky->PixelCoordinates(i);
                colorMap->data()->setCell(static_cast<int>(pix[0]),
                                          static_cast<int>(pix[1]),
                                          (*fsky)[i]);
            }

            QCPColorScale* colorScale = new QCPColorScale(qPlot);
            qPlot->plotLayout()->addElement(0, 1, colorScale);
            colorScale->setType(QCPAxis::atRight);
            colorMap->setColorScale(colorScale);
            colorMap->setGradient(QtFITS::makeViridisGradient());
            colorMap->rescaleDataRange();
            colorMap->setInterpolate(false);

            // ---- 2. Centroid cross at (x0, y0) = (params[1], params[2]) ----
            if(params.size() >= 3) {
                QCPGraph* centroid = new QCPGraph(xAxis, yAxis);
                centroid->setLineStyle(QCPGraph::lsNone);
                centroid->setScatterStyle(
                    QCPScatterStyle(QCPScatterStyle::ssCross,
                                   Qt::yellow, Qt::yellow, 14));
                centroid->addData(params[1], params[2]);
            }

            // ---- 3. Model contours scaled from model grid to fsky pixel coords ----
            // QCustomContours coords are in [0..ony-1][0..onx-1]; divide by oversample
            // to map to fsky pixel space [0..ny-1][0..nx-1].
            // NaN-separated polylines are split into individual QCPCurve objects.
            const auto& lines = contours.lines();
            for(int k = 0; k < static_cast<int>(lines.size()); k++) {
                const auto& [xs, ys] = lines[k];
                QVector<double> segX, segY, segT;
                double t = 0.0;
                auto flush = [&]() {
                    if(segX.size() > 1) {
                        QCPCurve* curve = new QCPCurve(xAxis, yAxis);
                        curve->setData(segT, segX, segY);
                        curve->setPen(QPen(Qt::white, 1.0));
                    }
                    segX.clear(); segY.clear(); segT.clear(); t = 0.0;
                };
                for(int j = 0; j < xs.size(); j++) {
                    if(std::isnan(xs[j])) {
                        flush();
                    } else {
                        segX.push_back(xs[j] * inv_s);
                        segY.push_back(ys[j] * inv_s);
                        segT.push_back(t++);
                    }
                }
                flush();
            }

            // ---- 4. Rescale axes, equal aspect ratio, render ----
            qPlot->rescaleAxes();
            qPlot->yAxis->setScaleRatio(qPlot->xAxis, 1.0);
            qPlot->replot();

            // ---- 5. Optional file save ----
            if(!saveFile.empty()) {
                namespace fs = std::filesystem;
                auto ext = fs::path(saveFile).extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                const auto qs = QString::fromStdString(saveFile);
                if(ext == ".jpg" || ext == ".jpeg")
                    qPlot->saveJpg(qs);
                else if(ext == ".pdf")
                    qPlot->savePdf(qs);
                else
                    qPlot->savePng(qs.endsWith(".png") ? qs : qs + ".png");
            }
        });
}

void starFitter::Plot(QCPAxisRect* axisRect,
                      const ROOT::Minuit2::MnUserParameters& mnParams,
                      size_t oversample,
                      const std::string& saveFile) const
{
    std::vector<double> vals;
    vals.reserve(mnParams.Parameters().size());
    for(const auto& p : mnParams.Parameters())
        vals.push_back(p.Value());
    Plot(axisRect, vals, oversample, saveFile);
}

}} // namespace DST::algorithms

#endif // __WITH_QT__
