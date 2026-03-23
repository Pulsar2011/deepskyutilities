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
#pragma region - DST::starFinder class implementation
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

#pragma endregion

#pragma region -- Accessor
    /** 
     * @details Get the evaluation method as string
     * @return \c std::string Evaluation method string
     */
    const std::string Stars::evalMethodString() const
    {
        switch (fevalMethod)
        {
            case evalType::MOMENT:
                return "MOMENT";
            case evalType::GAUSS:
                return "GAUSS-FIT";
            case evalType::GAUSS2D:
                return "2D-GAUSS-FIT";
            case evalType::MOFFAT:
                return "MOFFAT-FIT";
            case evalType::MOFFAT2D:
                return "2D-MOFFAT-FIT";
            default:
                return "UNKNOWN";
        }
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
        FITScolumn<uint64_t>  ID      ("ID", tulonglong, "", 0);
        FITScolumn<std::vector<double>> coo    ("coord_px", tdouble, "pixels", 1);
        FITScolumn<std::vector<double>> coo_max("maxcoord_px", tdouble, "pixels", 2);

        std::vector< FITScolumn<double> > wcsColumns;
        for(size_t i=0; i<wcsCoord.size(); i++)
            wcsColumns.push_back( FITScolumn<double>( wcsCoord[i], tdouble, wcsUnits[i], static_cast<int>(3+i) ) );

        FITScolumn<double>       maxima        ("maxima", tdouble, "", 4);
        FITScolumn<double>       FWHM_X        ("FWHM_X", tdouble, "pixels", 5);
        FITScolumn<double>       FWHM_Y        ("FWHM_Y", tdouble, "pixels", 6);
        FITScolumn<double>       angle         ("angle", tdouble, "radians", 7);
        FITScolumn<double>       TotalFlux     ("TotalFlux", tdouble, pixel_units, 8);
        FITScolumn<double>       SNR           ("SNR", tdouble, "", 9);
        FITScolumn<double>       Background    ("Background", tdouble, "", 10);
        FITScolumn<std::string>  EvalMethod    ("EvalMethod", tstring, "", 11);
        FITScolumn<uint32_t>     n_pixels      ("n_pixels", tuint, "", 12);
        FITScolumn<std::string>  source        ("source", tstring, "", 13);

        table.InsertColumn(std::make_shared< FITScolumn<uint64_t> >(ID));
        table.InsertColumn(std::make_shared< FITScolumn<std::vector<double>> >(coo));
        table.InsertColumn(std::make_shared< FITScolumn<std::vector<double>> >(coo_max));

        for(size_t i=0; i<wcsColumns.size(); i++)
            table.InsertColumn(std::make_shared< FITScolumn<double> >(wcsColumns[i]));
        
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

                iK->Write(foutURL+"/BS_opening.fits",true);
                iN->Write(foutURL+"/BM_dilation.fits",true);
                iS->Write(foutURL+"/BE_erosion.fits",true);
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
                iK->Write(foutURL+"/significance_map.fits",true);
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
                smoothedImage->Write(foutURL+"/thresholded_image.fits",true);
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
                errMap->Write(foutURL+"/error_map.fits",true);
                bkgMap->Write(foutURL+"/background_map.fits",true);
                snrMap->Write(foutURL+"/snr_map.fits",true);
            }

            errMap.reset();
            bkgMap.reset();

            // The third step consists in grouping connected pixels above the SNR threshold into star candidates using a watershed algorithm.


            _watershedFlooding(snrMap, fthresholdValue, fquantumStepQueue, iLayer);
            snrMap.reset();
        }

        void starFinder::_findStarsCustom(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {
            std::shared_ptr<FITScube> customMap = _customPreprocessing(img, iLayer);
            if(customMap == nullptr)
                throw FITSexception(BAD_DIMEN,"DST::starFinder","_findStarsCustom","Custom preprocessing returned a null FITScube pointer.["+std::to_string(__LINE__)+"]");

            if(foutURL != "")
            {
                customMap->Write(foutURL+"/userSignificance_map.fits",true);
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
                    throw std::invalid_argument("DST::starFinder::EvalStellarProperties: Unknown evaluation method.["+std::to_string(__LINE__)+"]");
            }
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
            }

            if(std::isfinite(maxima))
            {
                star.setMaxima(maxima);
                star.setMaximaPosition(point::pointFromVector<size_t>(fsky->PixelCoordinates(maxIndex)));
            }
        }

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

                        double phi = iS.second.angle()*acos(-1)/180.;

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
