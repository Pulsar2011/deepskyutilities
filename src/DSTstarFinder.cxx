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

#include <queue>
#include <set>

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
        fBackground(std::numeric_limits<double>::quiet_NaN())
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
                                    fmethod(detectionMethod::MORPHO),fdetectedStars(),
                                    fthresholdValue(std::numeric_limits<double>::min()),
                                    fBS(Morphology::kernel<double>::one({3,3})),
                                    fBE(Morphology::kernel<double>::one({9,9})),
                                    fBM(Morphology::kernel<double>::one({13,13})),
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
        }

        starFinder::~starFinder()
        {
            fdetectedStars.clear();
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
                            _findStarsMorpho(layer,quantizationStep);
                            break;
                        case detectionMethod::THRESHOLD:
                            _findStarsThreshold(layer,quantizationStep);
                            break;
                        case detectionMethod::SNR:
                            _findStarsSNR(layer,quantizationStep);
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
            // Compute iK and iN in parallel
            std::future<std::shared_ptr<FITScube>> futureIK = std::async(std::launch::async, [&img, this]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 2)
                        nThreads = (nThreads - 2)/2;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    std::shared_ptr<FITScube> iK = Morphology::opening<double>(img, fBS,nThreads);
                    iK->MaskPixels(img->GetMask());
                    return iK;
                });

            std::future<std::shared_ptr<FITScube>> futureIN = std::async(std::launch::async, [&img, this]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 2)
                        nThreads = (nThreads - 2)/2;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    std::shared_ptr<FITScube> iN = Morphology::dilation<double>(img, fBM,nThreads);
                    iN->MaskPixels(img->GetMask());
                    return iN;
                });

            // Wait for both to complete
            std::shared_ptr<FITScube> iK = futureIK.get();
            std::shared_ptr<FITScube> iN = futureIN.get();

            std::shared_ptr<FITScube> iS = Morphology::erosion<double>(iN, fBE, maxThreads);
            iS->MaskPixels(img->GetMask());

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
            *(iK.get()) -= *(iS.get());
            double minVal = std::min(iK->GetMinimum(),iS->GetMinimum());
            *(iK.get()) -=  minVal; 

            for(size_t i=0; i<iK->Nelements(); ++i)
            {
                if(iS->DoubleValueAtPixel(i) <= std::numeric_limits<double>::min())
                    (*(iK->GetData<double>()))[i] *= 0;
            }

            if(foutURL != "")
            {
                iK->SetName("significance_map_layer_"+std::to_string(iLayer));
                iK->Write(foutURL+"/significance_map.fits",true);
            }

            iS.reset();
            iN.reset();

            fthresholdValue = iK->Getpercentil(fthresholdValue);

            // Search stars using watershedalgorithm on the significance map
            _watershedFlooding(iK, fquantumStepQueue, iLayer);
 
        }

        void starFinder::_findStarsThreshold(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {

            size_t maxThreads = std::thread::hardware_concurrency();
            
            if(maxThreads == 0)
                maxThreads = 1;

            // The first stepconsists in smoothing the image to reduses the noise level applying a mean filter on the image using the BS kernel.
            std::future<std::shared_ptr<FITScube>> futuresmoothedImage = std::async(std::launch::async, [&img, this]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 3)
                        nThreads = (nThreads - 3)/3;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    fBS.setType(DST::Morphology::ktype::MEAN);
                    std::shared_ptr<FITScube> smoothedImage = Morphology::folding<double>(img, fBS, nThreads);
                    smoothedImage->MaskPixels(img->GetMask());
                    return smoothedImage;
                });
                        
            std::future<std::shared_ptr<FITScube>> futurebkgMap = std::async(std::launch::async, [&img, this]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 3)
                        nThreads = (nThreads - 3)/3;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    fBE.setType(DST::Morphology::ktype::MEDIAN);

                    std::shared_ptr<FITScube> bkgMap = Morphology::folding<double>(img, fBE, nThreads);
                    bkgMap->MaskPixels(img->GetMask());
                    return bkgMap;
                });

            // Wait for both to complete
            std::shared_ptr<FITScube> smoothedImage = futuresmoothedImage.get();
            std::shared_ptr<FITScube> bkgMap        = futurebkgMap.get();

            // Subtract background
            *(smoothedImage.get()) -= *(bkgMap.get());

            if(foutURL != "")
            {
                smoothedImage->Write(foutURL+"/thresholded_image.fits",true);
            }

            // The second step consists in searching for connected components above the threshold value using a watershed algorithm.

            _watershedFlooding(smoothedImage, fquantumStepQueue,iLayer);

            smoothedImage.reset();
        }

        void starFinder::_findStarsSNR(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
        {
            // Define kernel to be used
            fBS.setType(DST::Morphology::ktype::RMSE);
            fBE.setType(DST::Morphology::ktype::MEDIAN);


            // The first step consists in estimating the background and noise level of the image.
            std::future<std::shared_ptr<FITScube>> futureErrMap = std::async(std::launch::async, [&img, this]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 3)
                        nThreads = (nThreads - 3)/3;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    std::shared_ptr<FITScube> errMap = Morphology::folding<double>(img, fBS, nThreads);
                    errMap->MaskPixels(img->GetMask());

                    return errMap;
                });

            
            std::future<std::shared_ptr<FITScube>> futurebkgMap = std::async(std::launch::async, [&img, this]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 3)
                        nThreads = (nThreads - 3)/3;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    std::shared_ptr<FITScube> bkgMap = Morphology::folding<double>(img, fBE, nThreads);
                    bkgMap->MaskPixels(img->GetMask());
                    return bkgMap;
                });

            // The second step consists in searching for pixels with SNR above a given threshold.

            std::future<std::shared_ptr<FITScube>> futuresnrMap = std::async(std::launch::async, [&img]()
                {
                    size_t nThreads = std::thread::hardware_concurrency();
                    if(nThreads > 3)
                        nThreads = (nThreads - 3)/3;
                    
                    if(nThreads == 0)
                        nThreads = 1;

                    Morphology::kernel<double> snrKernel = DST::Morphology::kernel<double>::one({1,1});
                    snrKernel.setType(DST::Morphology::ktype::SUM);

                    std::shared_ptr<FITScube> snrMap = Morphology::folding<double>(img, snrKernel);
                    snrMap->MaskPixels(img->GetMask());
                    return snrMap;
                });

            // Wait for both to complete
            std::shared_ptr<FITScube> errMap = futureErrMap.get();
            std::shared_ptr<FITScube> bkgMap = futurebkgMap.get();
            std::shared_ptr<FITScube> snrMap = futuresnrMap.get();

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


            _watershedFlooding(snrMap, fquantumStepQueue, iLayer);
            snrMap.reset();
        }

        /**
         * @brief Watershed flooding algorithm to detect stars
         * 
         */
        void starFinder::_watershedFlooding(const std::shared_ptr<FITScube>& img, const double& fquantumStepQueue, const size_t& iLayer)
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
                std::cout << "DST::starFinder::_watershedFlooding: Starting watershed flooding on layer "<<iLayer<<" with max value "<<maxVal<<" and threshold value "<<fthresholdValue<<"."<<std::endl;
            }

            for(double threshold=maxVal; threshold >= fthresholdValue; threshold-=fquantumStepQueue)
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
                        // New star detected
                        size_t newLabel = fdetectedStars.size() + 1;
                        labels[pixIdx] = newLabel;
                        fdetectedStars.insert(std::pair<size_t,Stars>(newLabel,Stars(newLabel)));
                        std::vector<point> belong2Star = _floodingPlateau(pixIdx, labels, img, threshold, newLabel);
                        if(is3D)
                        {
                            for(std::vector<point>::iterator p=belong2Star.begin(); p!=belong2Star.end(); ++p)
                                p->SetZ(static_cast<double>(iLayer));
                        }
                        fdetectedStars[newLabel].addPixel(belong2Star);

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
                        fdetectedStars[existingLabel].addPixel(belong2Star);
                    }
                }
            }

            if(foutURL != "")
            {
                FITSimg<int64_t> labelImg(img->GetDimension(),{img->Size(1),img->Size(2)});
                labelImg.SetName("label_image_layer_"+std::to_string(iLayer));
                *(labelImg.template GetData<int64_t>()) *= 0;
                for(size_t i=0; i<labels.size(); ++i)
                {
                    if(labels[i] <= 0)
                        continue;
                    (*(labelImg.GetData<int64_t>()))[i]=static_cast<int64_t>(labels[i]);
                }
                labelImg.Write(foutURL+"/label_image.fits",true);
            }
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
    }
#pragma endregion
    
} // namespace name
