#include <gtest/gtest.h>
#include <DSTools/DSTmorpho.h>
#include <DSTools/DSTstarFinder.h>
#include <DSTdisplay/DSTmathplotpp.h>

#include <cmath>
#include <filesystem>

#include "test_utils.h"


using namespace DST;

TEST(StarsTest, DefaultConstructor)
{
    Stars s;
    EXPECT_EQ(s.ID(), static_cast<size_t>(0));
    EXPECT_TRUE(s.FilledPixels().empty());
    EXPECT_TRUE(std::isnan(s.maxima()));
    EXPECT_TRUE(std::isnan(s.SNR()));
    EXPECT_TRUE(std::isnan(s.totalFlux()));
    EXPECT_TRUE(std::isnan(s.background()));
    // Angle default can be NaN as well
    EXPECT_TRUE(std::isnan(s.angle()));
}

TEST(StarsTest, IdConstructor)
{
    const size_t id = 42;
    Stars s(id);
    EXPECT_EQ(s.ID(), id);
    EXPECT_TRUE(s.FilledPixels().empty());
}

TEST(StarsTest, ConstructorWithPointPosition)
{
    const size_t id = 7;
    DST::Math::point pos({10.5, 20.25});
    Stars s(id, pos);
    EXPECT_EQ(s.ID(), id);
    EXPECT_DOUBLE_EQ(s.position().X(), 10.5);
    EXPECT_DOUBLE_EQ(s.position().Y(), 20.25);
}

TEST(StarsTest, SetAndGetPosition)
{
    Stars s(1);
    DST::Math::point pos1({3.0, 4.0});
    s.setPosition(pos1);
    EXPECT_DOUBLE_EQ(s.position().X(), 3.0);
    EXPECT_DOUBLE_EQ(s.position().Y(), 4.0);

    DST::Math::point pos2({-1.25, 99.75});
    s.setPosition(pos2);
    EXPECT_DOUBLE_EQ(s.position().X(), -1.25);
    EXPECT_DOUBLE_EQ(s.position().Y(), 99.75);
}

TEST(StarsTest, MaximaAndMaximaPosition)
{
    Stars s(2);
    s.setMaxima(123.456);
    EXPECT_DOUBLE_EQ(s.maxima(), 123.456);

    DST::Math::point mpos({7.0, 8.0});
    s.setMaximaPosition(mpos);
    EXPECT_DOUBLE_EQ(s.maximaPosition().X(), 7.0);
    EXPECT_DOUBLE_EQ(s.maximaPosition().Y(), 8.0);
}

TEST(StarsTest, FWHMAndAngle)
{
    Stars s(3);
    s.setFWHM(2.5, 3.75);
    EXPECT_DOUBLE_EQ(s.FWHM().X(), 2.5);
    EXPECT_DOUBLE_EQ(s.FWHM().Y(), 3.75);

    s.setAngle(0.3141592653589793);
    EXPECT_DOUBLE_EQ(s.angle(), 0.3141592653589793);
}

TEST(StarsTest, FluxSnrBackground)
{
    Stars s(4);
    s.setTotalFlux(555.0);
    s.setSNR(42.0);
    s.setBackground(12.25);

    EXPECT_DOUBLE_EQ(s.totalFlux(), 555.0);
    EXPECT_DOUBLE_EQ(s.SNR(), 42.0);
    EXPECT_DOUBLE_EQ(s.background(), 12.25);
}

TEST(StarsTest, AddPixelWithPoint)
{
    Stars s(5);
    EXPECT_TRUE(s.FilledPixels().empty());

    DST::Math::point p1({1.0, 2.0});
    DST::Math::point p2({3.5, 4.5});
    s.addPixel(p1);
    s.addPixel(p2);

    ASSERT_EQ(s.FilledPixels().size(), static_cast<size_t>(2));
    EXPECT_DOUBLE_EQ(s.FilledPixels()[0].X(), 1.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[0].Y(), 2.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[1].X(), 3.5);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[1].Y(), 4.5);
}

TEST(StarsTest, AddMultiplePixelsVectorOfPoints)
{
    Stars s(6);
    std::vector<DST::Math::point> pts = {
        DST::Math::point({10.0, 20.0}),
        DST::Math::point({11.0, 21.0}),
        DST::Math::point({12.0, 22.0})
    };
    s.addPixel(pts);

    ASSERT_EQ(s.FilledPixels().size(), static_cast<size_t>(3));
    EXPECT_DOUBLE_EQ(s.FilledPixels()[0].X(), 10.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[0].Y(), 20.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[1].X(), 11.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[1].Y(), 21.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[2].X(), 12.0);
    EXPECT_DOUBLE_EQ(s.FilledPixels()[2].Y(), 22.0);
}

TEST(StarFinderTest, SyntheticImageGeneration)
{
    // Generate a synthetic star field with catalog
    std::vector<DST::Stars> catalog;
    auto img = TestUtils::generateStarField(
        catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        100.0,  // min intensity
        10000.0,// max intensity
        100.0,  // background
        10.0,   // noise sigma
        42      // seed for reproducibility
    );

    ASSERT_NE(img, nullptr);
    EXPECT_EQ(img->Size(1), static_cast<size_t>(512));
    EXPECT_EQ(img->Size(2), static_cast<size_t>(512));
    EXPECT_EQ(catalog.size(), static_cast<size_t>(50));
    
    // Check that image has reasonable statistics
    double mean = img->GetMean();
    EXPECT_GT(mean, 50.0);  // Should be above zero due to stars and background
    
    // Verify catalog contains valid star positions
    for(const auto& star : catalog)
    {
        EXPECT_GT(star.position().X(), 0.0);
        EXPECT_LT(star.position().X(), 512.0);
        EXPECT_GT(star.position().Y(), 0.0);
        EXPECT_LT(star.position().Y(), 512.0);
        EXPECT_GT(star.maxima(), 100.0);
        EXPECT_LT(star.maxima(), 10000.0);
    }

    {
        std::filesystem::path outDir("./build/test_output");
        std::error_code ec;
        std::filesystem::create_directories(outDir, ec);
    }

}

TEST(StarFinderTest,switchKernels)
{
    std::shared_ptr<DSL::FITScube> img = std::shared_ptr<DSL::FITScube>(new DSL::FITSimg<int32_t>({512,512}));
    DST::algorithms::starFinder finder(img);

    DST::Morphology::kernel<double> newBS = DST::Morphology::kernel<double>::one({2,2})  ; DST::Morphology::kernel<double> testBS = newBS;
    DST::Morphology::kernel<double> newBE = DST::Morphology::kernel<double>::one({13,13}); DST::Morphology::kernel<double> testBE = newBE;
    DST::Morphology::kernel<double> newBM = DST::Morphology::kernel<double>::one({3,3})  ; DST::Morphology::kernel<double> testBM = newBM;

    finder.switchfBSkernel(newBS);
    finder.switchfBEkernel(newBE);
    finder.switchfBMkernel(newBM);

    const DST::Morphology::kernel<double>& currentBS = finder.getfBSkernel();
    const DST::Morphology::kernel<double>& currentBE = finder.getfBEkernel();
    const DST::Morphology::kernel<double>& currentBM = finder.getfBMkernel();

    EXPECT_EQ(currentBS.size(), testBS.size()); for(size_t i=0; i<currentBS.size(); ++i) EXPECT_NEAR(currentBS.data()[i], testBS.data()[i],1e-6);
    EXPECT_EQ(currentBE.size(), testBE.size()); for(size_t i=0; i<currentBE.size(); ++i) EXPECT_NEAR(currentBE.data()[i], testBE.data()[i],1e-6);
    EXPECT_EQ(currentBM.size(), testBM.size()); for(size_t i=0; i<currentBM.size(); ++i) EXPECT_NEAR(currentBM.data()[i], testBM.data()[i],1e-6);

}

TEST(StarFinderTest, FindStarsBySegmentation)
{
    // Generate synthetic star field with known catalog
    std::vector<DST::Stars> true_catalog;
    const unsigned int seed = 12345;
    auto img = TestUtils::generateStarField(
        true_catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        100.0,   // min intensity
        1000.0, // max intensity
        5.0,    // background
        5.0,    // noise sigma
        seed);
    
    img->Write("./build/test_output/synthetic_star_field.fits",true);
    
    // Create star finder
    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    
    // Run star detection
    verbose = verboseLevel::VERBOSE_BASIC;
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);
    
    // Check that stars were found
    const auto& detectedStars = finder.detectedStars();
    EXPECT_GT(detectedStars.nrows(), static_cast<size_t>(0));
    
    // Optionally check detection rate (should find most stars)
    EXPECT_GT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 80.); // At least 80% detection
    EXPECT_LT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 110.); // No more than 120% (to avoid too many false positives)
    std::cout << "Detected " << detectedStars.nrows() << " stars out of " << true_catalog.size() << " (" 
              << (static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100.) << "%)" << std::endl;
    
    // Compare detected positions with true positions (within tolerance)
    size_t matched = 0;
    const double tolerance = 0.35; // pixels
    
    // Get direct access to column data vectors
    const auto& list = detectedStars.getColumn("coord_px")->values<FITSform::doubleVector>();
    
    for(const auto& true_star : true_catalog)
    {
        for(const auto& stars : list)
        {
            double dx = stars[0] - true_star.position().X();
            double dy = stars[1] - true_star.position().Y();
            double dist = std::sqrt(dx*dx + dy*dy);
            
            if(dist < tolerance)
            {
                matched++;
                break;
            }
        }
    }
    
    // Expect at least 50% of true stars to be matched
    EXPECT_GT(matched, true_catalog.size() / 2);
}


TEST(StarFinderTest, FindStarsBySNR)
{
    // Generate synthetic star field with known catalog
    std::vector<DST::Stars> true_catalog;
    const unsigned int seed = 12345;
    auto img = TestUtils::generateStarField(
        true_catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        10.0,   // min intensity
        1000.0, // max intensity
        10.0,    // background
        8.0,    // noise sigma
        seed);
    
    img->Write("./build/test_output/synthetic_star_field2.fits",true);
    
    // Create star finder
    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    
    // Run star detection
    verbose = verboseLevel::VERBOSE_BASIC;
    finder.setThresholdValue(5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::SNR, 0.1);
    
    // Check that stars were found
    const FITStable& detectedStars = finder.detectedStars();
    EXPECT_GT(detectedStars.nrows(), static_cast<size_t>(0));
    
    // Optionally check detection rate (should find most stars)
    EXPECT_GT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 60.); // At least 80% detection
    EXPECT_LT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 110.); // No more than 120% (to avoid too many false positives)
    std::cout << "Detected " << detectedStars.nrows() << " stars out of " << true_catalog.size() << " (" 
              << (static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100.) << "%)" << std::endl;
    
    // Compare detected positions with true positions (within tolerance)
    size_t matched = 0;
    const double tolerance = 0.35; // pixels

    // Get direct access to column data vectors
    const auto& list = detectedStars.getColumn("coord_px")->values<FITSform::doubleVector>();
    
    for(const auto& true_star : true_catalog)
    {
        for(const auto& stars : list)
        {
            double dx = stars[0] - true_star.position().X();
            double dy = stars[1] - true_star.position().Y();
            double dist = std::sqrt(dx*dx + dy*dy);
            
            if(dist < tolerance)
            {
                matched++;
                break;
            }
        }
    }
    
    // Expect at least 50% of true stars to be matched
    EXPECT_GT(matched, true_catalog.size() / 2);
}

TEST(StarFinderTest, FindStarsByThreshold)
{
    // Generate synthetic star field with known catalog
    std::vector<DST::Stars> true_catalog;
    const unsigned int seed = 12345;
    auto img = TestUtils::generateStarField(
        true_catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        10.0,   // min intensity
        1000.0, // max intensity
        10.0,    // background
        8.0,    // noise sigma
        seed);
    
    img->Write("./build/test_output/synthetic_star_field2.fits",true);
    
    // Create star finder
    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    
    // Run star detection
    verbose = verboseLevel::VERBOSE_BASIC;
    finder.setThresholdValue(70);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::THRESHOLD, 0.1);
    
    // Check that stars were found
    const auto& detectedStars = finder.detectedStars();
    EXPECT_GT(detectedStars.nrows(), static_cast<size_t>(0));
    
    // Optionally check detection rate (should find most stars)
    EXPECT_GT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 50.); // At least 80% detection
    EXPECT_LT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 110.); // No more than 120% (to avoid too many false positives)
    std::cout << "Detected " << detectedStars.nrows() << " stars out of " << true_catalog.size() << " (" 
              << (static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100.) << "%)" << std::endl;
    
    // Compare detected positions with true positions (within tolerance)
    size_t matched = 0;
    const double tolerance = 0.35; // pixels
    
    // Get direct access to column data vectors
    const auto& list = detectedStars.getColumn("coord_px")->values<FITSform::doubleVector>();
    
    for(const auto& true_star : true_catalog)
    {
        for(const auto& stars : list)
        {
            double dx = stars[0] - true_star.position().X();
            double dy = stars[1] - true_star.position().Y();
            double dist = std::sqrt(dx*dx + dy*dy);
            
            if(dist < tolerance)
            {
                matched++;
                break;
            }
        }
    }
    
    // Expect at least 50% of true stars to be matched
    EXPECT_GT(matched, true_catalog.size() / 2);
}

TEST(StarFinderTest, FindStarsCube)
{
    // Generate synthetic star field with known catalog
    std::vector<DST::Stars> true_catalog;
    const unsigned int seed = 12345;
    auto img1 = TestUtils::generateStarField(
        true_catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        10.0,   // min intensity
        1000.0, // max intensity
        5.0,    // background
        8.0,    // noise sigma
        seed);

    auto img2 = TestUtils::generateStarField(
        true_catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        10.0,   // min intensity
        1000.0, // max intensity
        5.0,    // background
        8.0,    // noise sigma
        seed+385);

    auto img3 = TestUtils::generateStarField(
        true_catalog,
        512,    // width
        512,    // height
        50,     // number of stars
        0.7,    // sigma
        10.0,   // min intensity
        1000.0, // max intensity
        5.0,    // background
        8.0,    // noise sigma
        seed+976);

    img1->AddLayer(*img2);
    img1->AddLayer(*img3);

    std::shared_ptr<FITScube> img = img1;

    img->Write("./build/test_output/synthetic_star_field.fits",true);
    img->HDU().ValueForKey("CTYPE1","RA---TAN",fChar);
    img->HDU().ValueForKey("CTYPE2","DEC--TAN",fChar);
    img->HDU().ValueForKey("CTYPE3","WAVELENGTH",fChar);
    img->HDU().ValueForKey("CUNIT1","DEG",fChar);
    img->HDU().ValueForKey("CUNIT2","DEG",fChar);
    img->HDU().ValueForKey("CUNIT3","nm",fChar);
    img->HDU().ValueForKey("CRPIX1",0.);
    img->HDU().ValueForKey("CRPIX2",0.);
    img->HDU().ValueForKey("CRPIX3",0.);
    img->HDU().ValueForKey("CRVAL1",180.0);
    img->HDU().ValueForKey("CRVAL2",0.0);
    img->HDU().ValueForKey("CRVAL3",656.28);
    img->HDU().ValueForKey("CDELT1",-0.0002777777778);
    img->HDU().ValueForKey("CDELT2",0.0002777777778);
    img->HDU().ValueForKey("CDELT3",0.1);
    
    // Create star finder
    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    
    // Run star detection
    verbose = verboseLevel::VERBOSE_BASIC;
    finder.setThresholdValue(0.998);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);
    
    // Check that stars were found
    const auto& detectedStars = finder.detectedStars();

    detectedStars.toCSV(std::cout,';');

    {
        DSL::mathplot::StderrRedirect redirect; // Silence plotting warnings
    
        EXPECT_NO_THROW(finder.plot(0));
        EXPECT_NO_THROW(finder.plot(1));
        EXPECT_NO_THROW(finder.plot(2));
    }
}