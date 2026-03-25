#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <memory>
#include <random>
#include <cmath>
#include <vector>
#include <functional>
#include <utility>
#include <DSTfits/FITSimg.h>
#include <DSTools/DSTstarFinder.h>

namespace TestUtils
{
    /**
     * @brief Generate a 2D Gaussian at a given position
     * 
     * @param img Image to add the Gaussian to
     * @param star Star object containing parameters (position, intensity, width)
     */
    template<typename T>
    void addGaussianStar(DSL::FITSimg<T>& img, const DST::Stars& star)
    {
        const size_t nx = img.Size(1);
        const size_t ny = img.Size(2);
        
        // Extract star parameters
        const double x = star.position().X();
        const double y = star.position().Y();
        const double intensity = star.maxima();
        const double sigma = star.FWHM().X() / 2.355; // Convert FWHM to sigma
        
        // Define region of influence (5-sigma around center)
        const double extent = 5.0 * sigma;
        const int xmin = std::max(0, static_cast<int>(x - extent));
        const int xmax = std::min(static_cast<int>(nx), static_cast<int>(x + extent + 1));
        const int ymin = std::max(0, static_cast<int>(y - extent));
        const int ymax = std::min(static_cast<int>(ny), static_cast<int>(y + extent + 1));
        
        const double sigma2 = sigma * sigma;
        
        for(int iy = ymin; iy < ymax; ++iy)
        {
            for(int ix = xmin; ix < xmax; ++ix)
            {
                const double dx = static_cast<double>(ix) - x;
                const double dy = static_cast<double>(iy) - y;
                const double r2 = dx*dx + dy*dy;
                const double gauss = intensity * std::exp(-0.5 * r2 / sigma2);
                
                // Get current pixel value and add Gaussian
                std::vector<size_t> coords = {static_cast<size_t>(ix), static_cast<size_t>(iy)};
                const size_t idx = img.PixelIndex(coords);
                const T current = img.DoubleValueAtPixel(idx);
                img.SetPixelValue(current + static_cast<T>(gauss), idx);
            }
        }
    }

    /**
     * @brief Generate a synthetic FITS image with random Gaussian stars
     * 
     * @param width Image width in pixels
     * @param height Image height in pixels
     * @param nStars Number of stars to generate
     * @param sigma_stars Gaussian width of stars (in pixels)
     * @param min_intensity Minimum star intensity
     * @param max_intensity Maximum star intensity
     * @param background Background level
     * @param noise_sigma Standard deviation of Gaussian noise
     * @param seed Random seed (0 = use time-based seed)
     * @param catalog Output parameter: vector to store generated star catalog
     * @return std::shared_ptr<DSL::FITSimg<double>> Generated image
     */
    inline std::shared_ptr< DSL::FITSimg<double> > generateStarField(
        std::vector<DST::Stars>& catalog,
        size_t width = 512,
        size_t height = 512,
        size_t nStars = 50,
        double sigma_stars = 0.7,
        double min_intensity = 10.0,
        double max_intensity = 10000.0,
        double background = 10.0,
        double noise_sigma = 5.0,
        unsigned int seed = 0)
    {
        // Create empty image
        auto img = std::make_shared<DSL::FITSimg<double>>(std::vector<size_t>{width, height});
        
        // Initialize random number generator
        std::mt19937_64 rng(seed == 0 ? 
            static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) : seed);
        
        // Distributions
        std::uniform_real_distribution<double> x_dist(5.0, static_cast<double>(width) - 5.0);
        std::uniform_real_distribution<double> y_dist(5.0, static_cast<double>(height) - 5.0);
        std::uniform_real_distribution<double> intensity_dist(min_intensity, max_intensity);
        std::normal_distribution<double> noise_dist(background, noise_sigma);
        
        // Add background with noise
        for(size_t i = 0; i < img->Nelements(); ++i)
        {
            const double pixel_value = noise_dist(rng);
            std::vector<size_t> coords = img->PixelCoordinates(i);

            double x0=-25, y0=-105;
            const double diffuse_background = background * 100 * exp(-0.5*((coords[0]-x0)*(coords[0]-x0)+(coords[1]-y0)*(coords[1]-y0))/(430.0*430.0));

            img->SetPixelValue(pixel_value+diffuse_background, i);
        }
        
        // Generate random stars
        catalog.clear();
        catalog.reserve(nStars);
        
        const double fwhm = sigma_stars * 2.355; // Convert sigma to FWHM
        
        for(size_t i = 0; i < nStars; ++i)
        {
            double x = x_dist(rng);
            double y = y_dist(rng);
            double intensity = intensity_dist(rng);
            
            // Create Star object with generated parameters
            DST::Stars star(i + 1, DST::Math::point({x, y}));
            star.setMaxima(intensity);
            star.setFWHM(fwhm, fwhm);
            star.setBackground(background);
            
            catalog.push_back(star);
            addGaussianStar(*img, star);
        }
        
        // Update FITS header
        img->HDU().ValueForKey("OBJECT", "Synthetic Star Field", DSL::fChar, "Test image");
        img->HDU().ValueForKey("NSTARS", static_cast<int>(nStars), "Number of synthetic stars");
        img->HDU().ValueForKey("STARSIG", sigma_stars, "Star Gaussian sigma [pixels]");
        img->HDU().ValueForKey("BKGND", background, "Background level");
        img->HDU().ValueForKey("NOISE", noise_sigma, "Noise sigma");
        
        return img;
    }

    /**
     * @brief Overload without catalog output for backward compatibility
     */
    inline std::shared_ptr<DSL::FITSimg<double>> generateStarField(
        size_t width = 512,
        size_t height = 512,
        size_t nStars = 50,
        double sigma_stars = 0.7,
        double min_intensity = 10.0,
        double max_intensity = 10000.0,
        double background = 10.0,
        double noise_sigma = 5.0,
        unsigned int seed = 0)
    {
        std::vector<DST::Stars> dummy_catalog;
        return generateStarField(dummy_catalog, width, height, nStars, sigma_stars, 
                                min_intensity, max_intensity, background, noise_sigma, seed);
    }

    /**
     * @brief Create a flat-background image with stars rendered by a model function.
     *
     * @param width   Image width in pixels
     * @param height  Image height in pixels
     * @param background Uniform background level added to every pixel
     * @param positions  List of (x0, y0) star centres in pixel coordinates
     * @param starModel  Callable(x, y, x0, y0) returning the star value ABOVE background
     * @return Shared pointer to the generated FITSimg<double>
     */
    inline std::shared_ptr<DSL::FITSimg<double>>
    makeModelImage(size_t width, size_t height,
                   double background,
                   const std::vector<std::pair<double,double>>& positions,
                   std::function<double(double /*x*/, double /*y*/,
                                        double /*x0*/, double /*y0*/)> starModel)
    {
        auto img = std::make_shared<DSL::FITSimg<double>>(std::vector<size_t>{width, height});
        for (size_t k = 0; k < img->Nelements(); ++k)
        {
            auto pix = img->PixelCoordinates(k);
            double x = static_cast<double>(pix[0]);
            double y = static_cast<double>(pix[1]);
            double val = background;
            for (const auto& pos : positions)
                val += starModel(x, y, pos.first, pos.second);
            img->SetPixelValue(val, k);
        }
        return img;
    }

} // namespace TestUtils

#endif // TEST_UTILS_H
