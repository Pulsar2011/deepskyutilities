#include <gtest/gtest.h>
#include <DSTools/DSTmorpho.h>
#include <DSTools/DSTstarFinder.h>
#include <DSTools/DSTbackground.h>
#include <DSTfits/FITSmanager.h>

#ifdef __WITH_QT__
#include <DSTdisplay/QtFITSviewer.h>
#endif

#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>

#include "test_utils.h"


using namespace DST;

// ---------------------------------------------------------------------------
// Global test environment: wipe and recreate the output directory once before
// any test runs so stale files from previous runs never cause failures.
// ---------------------------------------------------------------------------
class TestOutputEnv : public ::testing::Environment
{
public:
    void SetUp() override
    {
        const std::filesystem::path outDir = "./build/test_output";
        std::filesystem::remove_all(outDir);
        std::filesystem::create_directories(outDir);
    }
};

const ::testing::Environment* const kTestOutputEnv =
    ::testing::AddGlobalTestEnvironment(new TestOutputEnv);

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

// ---------------------------------------------------------------------------
// starFitter::Eval() / operator() consistency
// ---------------------------------------------------------------------------

TEST(StarFitterFunction, EvalConsistency_gauss)
{
    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS params: A=500, x0=5, y0=5, sigma=1.5, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3])
                         + params[4];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::GAUSS);

    // Perfect data → chi2 should be 0
    EXPECT_NEAR(fitter(params), 0.0, 1e-8);

    // At the centroid gauss2D = A → Eval = A + background
    for(size_t i = 0; i < W; i++)
    {
        for(size_t j = 0; j < H; j++)
        {
            EXPECT_NEAR(fitter.Eval(static_cast<double>(i), static_cast<double>(j), params), (*img)[img->PixelIndex({i,j})] , 1e-8);
        }
    }
}

TEST(StarFitterFunction, EvalConsistency_igauss)
{
    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS params: A=500, x0=5, y0=5, sigma=1.5, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::igauss2D(x, y, params[0], params[1], params[2], params[3], (double)1., (double) 1.)
                         + params[4];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::iGAUSS);

    // Perfect data → chi2 should be 0
    EXPECT_NEAR(fitter(params), 0.0, 1e-8);

    // At the centroid gauss2D = A → Eval = A + background
    for(size_t i = 0; i < W; i++) {
        for(size_t j = 0; j < H; j++) {
            EXPECT_NEAR(fitter.Eval(static_cast<double>(i), static_cast<double>(j), params), (*img)[img->PixelIndex({i,j})] , 1e-8);
        }
    }
}

TEST(StarFitterFunction, EvalConsistency_gauss2D)
{
    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS2D params: A=500, x0=5, y0=5, sigmaX=1.5, sigmaY=2.5, theta=pi/3, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, DST::Math::MathCore::Pi()/3., 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5])
                         + params[6];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::GAUSS2D);

    // Perfect data → chi2 should be 0
    EXPECT_NEAR(fitter(params), 0.0, 1e-8);

    // At the centroid gauss2D = A → Eval = A + background
    for(size_t i = 0; i < W; i++) {
        for(size_t j = 0; j < H; j++) {
            EXPECT_NEAR(fitter.Eval(static_cast<double>(i), static_cast<double>(j), params), (*img)[img->PixelIndex({i,j})] , 1e-8);
        }
    }
}

TEST(StarFitterFunction, EvalConsistency_igauss2D)
{
    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS2D params: A=500, x0=5, y0=5, sigmaX=1.5, sigmaY=2.5, theta=pi/3, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, DST::Math::MathCore::Pi()/3., 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::igauss2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5],(double) 1., (double) 1.)
                         + params[6];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::iGAUSS2D);

    // Perfect data → chi2 should be 0
    EXPECT_NEAR(fitter(params), 0.0, 1e-8);

    // At the centroid gauss2D = A → Eval = A + background
    for(size_t i = 0; i < W; i++) {
        for(size_t j = 0; j < H; j++) {
            EXPECT_NEAR(fitter.Eval(static_cast<double>(i), static_cast<double>(j), params), (*img)[img->PixelIndex({i,j})] , 1e-8);
        }
    }
}

TEST(StarFitterFunction, EvalConsistency_moffat)
{
    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // MOFFAT params: A=500, x0=5, y0=5, alfa=1.5, beta=2.5, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::moffat2D(x, y, params[0], params[1], params[2], params[3], params[4])
                         + params[5];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::MOFFAT);

    // Perfect data → chi2 should be 0
    EXPECT_NEAR(fitter(params), 0.0, 1e-8);

    // At the centroid gauss2D = A → Eval = A + background
    for(size_t i = 0; i < W; i++) {
        for(size_t j = 0; j < H; j++) {
            EXPECT_NEAR(fitter.Eval(static_cast<double>(i), static_cast<double>(j), params), (*img)[img->PixelIndex({i,j})] , 1e-8);
        }
    }
}

TEST(StarFitterFunction, EvalConsistency_moffat2D)
{
    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // MOFFAT params: A=500, x0=5, y0=5, alfa=1.5, beta=2.5, gamma=0.75, theta =pi/3, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, 0.75, DST::Math::MathCore::Pi()/3., 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::moffat2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5],params[6])
                         + params[7];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::MOFFAT2D);

    // Perfect data → chi2 should be 0
    EXPECT_NEAR(fitter(params), 0.0, 1e-8);

    // At the centroid gauss2D = A → Eval = A + background
    for(size_t i = 0; i < W; i++) {
        for(size_t j = 0; j < H; j++) {
            EXPECT_NEAR(fitter.Eval(static_cast<double>(i), static_cast<double>(j), params), (*img)[img->PixelIndex({i,j})] , 1e-8);
        }
    }
}

// ---------------------------------------------------------------------------
// starFitter::Fitting() consistency
// ---------------------------------------------------------------------------

TEST(StarFitterFit, EvalConsistency_gauss)
{
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    std::default_random_engine re;

    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS params: A=500, x0=5, y0=5, sigma=1.5, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3])
                         + params[4];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::GAUSS);

    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("A",          params[0]*(1+unif(re)), 1e-6);
    upar.Add("x0",         params[1]*(1+unif(re)), 1e-6);
    upar.Add("y0",         params[2]*(1+unif(re)), 1e-6);
    upar.Add("sigma",      params[3]*(1+unif(re)), 1e-6);
    upar.Add("background", params[4]*(1+unif(re)), 1e-6);

    ROOT::Minuit2::MnMigrad migrad(fitter, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad();
    EXPECT_TRUE(result.IsValid());

    for(size_t i = 0; i < params.size(); i++)
    {
        EXPECT_LT(std::abs((result.UserState().Value(i)-params[i])/params[i]), 0.1);
    }
}

TEST(StarFitterFit, EvalConsistency_igauss)
{
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    std::default_random_engine re;

    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS params: A=500, x0=5, y0=5, sigma=1.5, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::igauss2D(x, y, params[0], params[1], params[2], params[3], (double)1., (double) 1.)
                         + params[4];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::iGAUSS);

    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("A",          params[0]*(1+unif(re)), 1e-6);
    upar.Add("x0",         params[1]*(1+unif(re)), 1e-6);
    upar.Add("y0",         params[2]*(1+unif(re)), 1e-6);
    upar.Add("sigma",      params[3]*(1+unif(re)), 1e-6);
    upar.Add("background", params[4]*(1+unif(re)), 1e-6);

    ROOT::Minuit2::MnMigrad migrad(fitter, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad();
    EXPECT_TRUE(result.IsValid());

    for(size_t i = 0; i < params.size(); i++)
    {
        EXPECT_LT(std::abs((result.UserState().Value(i)-params[i])/params[i]), 0.1);
    }
}

TEST(StarFitterFit, EvalConsistency_gauss2D)
{
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    std::default_random_engine re;

    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS2D params: A=500, x0=5, y0=5, sigmaX=1.5, sigmaY=2.5, theta=pi/3, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, DST::Math::MathCore::Pi()/3., 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5])
                         + params[6];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::GAUSS2D);

    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("A",          params[0]*(1+unif(re)), 1e-6);
    upar.Add("x0",         params[1]*(1+unif(re)), 1e-6);
    upar.Add("y0",         params[2]*(1+unif(re)), 1e-6);
    upar.Add("sigmaX",     params[3]*(1+unif(re)), 1e-6);
    upar.Add("sigmaY",     params[4]*(1+unif(re)), 1e-6);
    upar.Add("Theta" ,     params[5]*(1+unif(re)), 1e-6);
    upar.Add("background", params[6]*(1+unif(re)), 1e-6);

    ROOT::Minuit2::MnMigrad migrad(fitter, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad();
    EXPECT_TRUE(result.IsValid());

    for(size_t i = 0; i < params.size(); i++)
    {
        EXPECT_LT(std::abs((result.UserState().Value(i)-params[i])/params[i]), 0.1);
    }
}

TEST(StarFitterFit, EvalConsistency_igauss2D)
{
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    std::default_random_engine re;

    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // GAUSS2D params: A=500, x0=5, y0=5, sigmaX=1.5, sigmaY=2.5, theta=pi/3, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, DST::Math::MathCore::Pi()/3., 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::igauss2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5],(double) 1., (double) 1.)
                         + params[6];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::iGAUSS2D);

    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("A",          params[0]*(1+unif(re)), 1e-6);
    upar.Add("x0",         params[1]*(1+unif(re)), 1e-6);
    upar.Add("y0",         params[2]*(1+unif(re)), 1e-6);
    upar.Add("sigmaX",     params[3]*(1+unif(re)), 1e-6);
    upar.Add("sigmaY",     params[4]*(1+unif(re)), 1e-6);
    upar.Add("Theta" ,     params[5]*(1+unif(re)), 1e-6);
    upar.Add("background", params[6]*(1+unif(re)), 1e-6);

    ROOT::Minuit2::MnMigrad migrad(fitter, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad();
    EXPECT_TRUE(result.IsValid());

    for(size_t i = 0; i < params.size(); i++)
    {
        EXPECT_LT(std::abs((result.UserState().Value(i)-params[i])/params[i]), 0.1);
    }
}

TEST(StarFitterFit, EvalConsistency_moffat)
{
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    std::default_random_engine re;

    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // MOFFAT params: A=500, x0=5, y0=5, alfa=1.5, beta=2.5, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::moffat2D(x, y, params[0], params[1], params[2], params[3], params[4])
                         + params[5];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::MOFFAT);

    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("A",          params[0]*(1+unif(re)), 1e-6);
    upar.Add("x0",         params[1]*(1+unif(re)), 1e-6);
    upar.Add("y0",         params[2]*(1+unif(re)), 1e-6);
    upar.Add("alpha",      params[3]*(1+unif(re)), 1e-6);
    upar.Add("beta",       params[4]*(1+unif(re)), 1e-6);
    upar.Add("background", params[5]*(1+unif(re)), 1e-6);

    ROOT::Minuit2::MnMigrad migrad(fitter, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad();
    EXPECT_TRUE(result.IsValid());

    for(size_t i = 0; i < params.size(); i++)
    {
        EXPECT_LT(std::abs((result.UserState().Value(i)-params[i])/params[i]), 0.1);
    }
}

TEST(StarFitterFit, EvalConsistency_moffat2D)
{
    std::uniform_real_distribution<double> unif(-0.1,0.1);
    std::default_random_engine re;

    // Build an 11×11 image filled with exact GAUSS values (no noise → chi2 must be ~0)
    const size_t W = 11, H = 11;
    auto img = std::make_shared< DSL::FITSimg<double> >(std::vector<size_t>{W, H});

    // MOFFAT params: A=500, x0=5, y0=5, alfa=1.5, beta=2.5, gamma=0.75, theta =pi/3, background=10
    const std::vector<double> params = {500., 5., 5., 1.5, 2.5, 0.75, DST::Math::MathCore::Pi()/3., 10.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::moffat2D(x, y, params[0], params[1], params[2], params[3], params[4], params[5],params[6])
                         + params[7];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::MOFFAT2D);

    ROOT::Minuit2::MnUserParameters upar;
    upar.Add("A",          params[0]*(1+unif(re)), 1e-6);
    upar.Add("x0",         params[1]*(1+unif(re)), 1e-6);
    upar.Add("y0",         params[2]*(1+unif(re)), 1e-6);
    upar.Add("alpha",      params[3]*(1+unif(re)), 1e-6);
    upar.Add("beta",       params[4]*(1+unif(re)), 1e-6);
    upar.Add("gamma",      params[5]*(1+unif(re)), 1e-6);
    upar.Add("theta",      params[6]*(1+unif(re)), 1e-6);
    upar.Add("background", params[7]*(1+unif(re)), 1e-6);

    ROOT::Minuit2::MnMigrad migrad(fitter, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad();
    EXPECT_TRUE(result.IsValid());

    for(size_t i = 0; i < params.size(); i++)
    {
        EXPECT_LT(std::abs((result.UserState().Value(i)-params[i])/params[i]), 0.1);
    }
}

// ---------------------------------------------------------------------------
// starFitter::Plot() — visual / file output
// ---------------------------------------------------------------------------

TEST(StarFitter, Plot)
{
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    const size_t W = 11, H = 11;
    auto img = std::make_shared<DSL::FITSimg<double>>(std::vector<size_t>{W, H});

    const std::vector<double> params = {500., 5., 5., 1.5, 2.5,DST::Math::MathCore::Pi()/3.,0.};

    for(size_t i = 0; i < img->Nelements(); i++) {
        auto pix = img->PixelCoordinates(i);
        const double x = static_cast<double>(pix[0]);
        const double y = static_cast<double>(pix[1]);
        const double v = DST::Math::function::gauss2D(x, y, params[0], params[1], params[2], params[3],params[4],params[5])
                         + params[6];
        img->SetPixelValue(v, i);
    }

    DST::algorithms::starFitter fitter(img, DST::Stars::evalType::GAUSS2D);

    int arg = 0;
    QApplication app(arg, nullptr);

    QtFITSviewer* mw = new QtFITSviewer(250, 250, 600, 600);
    mw->show();
    QCustomPlot* customPlot = mw->GetPlot();
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    const std::string outFile = "./build/test_output/starfitter_plot.png";
    EXPECT_NO_THROW(fitter.Plot(customPlot->axisRect(), params, 10, outFile));

    std::ifstream file(outFile);
    EXPECT_TRUE(file.good());

    //app.exec();
#endif
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
    finder.setOutputUrl("./build/test_output/FindStarsBySegmentation");
    
    // Run star detection
    verbose = verboseLevel::VERBOSE_BASIC;
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);
    
    // Check that stars were found
    const auto& detectedStars = finder.detectedStars();
    EXPECT_GT(detectedStars.nrows(), static_cast<size_t>(0));
    
    // Optionally check detection rate (should find most stars)
    EXPECT_GT(static_cast<double>(detectedStars.nrows())/static_cast<double>(true_catalog.size())*100., 75.); // At least 80% detection
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
    finder.setOutputUrl("./build/test_output/FindStarsBySNR");
    
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
    finder.setOutputUrl("./build/test_output/FindStarsByThreshold");
    
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

TEST(StarFinderTest, FindStarsCustom)
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
        5.0,    // background
        8.0,    // noise sigma
        seed);

        // Create star finder
    DST::algorithms::starFinder finder(img);
    EXPECT_ANY_THROW(finder.findStars(DST::algorithms::starFinder::detectionMethod::CUSTOM, 0.1));

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
    finder.setOutputUrl("./build/test_output/FindStarsCube");
    
    // Run star detection
    verbose = verboseLevel::VERBOSE_BASIC;
    finder.setThresholdValue(0.998);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);
    
    // Check that stars were found
    const auto& detectedStars = finder.detectedStars();

    detectedStars.toCSV(std::cout,';');
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        for(size_t i=1; i<4; ++i)
        {
            QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
            mw->setWindowTitle(img->GetName().c_str());
            mw->show();

            QCustomPlot *customPlot = mw->GetPlot();
            customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

            QtFITSviewer::QtDisplayLayer(customPlot,img,i);
            EXPECT_NO_THROW(finder.Plot(customPlot->axisRect(), i));
            customPlot->savePng(("./build/test_output/detected_stars_layer"+std::to_string(i)+".png").c_str());

            std::ifstream file("./build/test_output/detected_stars_layer"+std::to_string(i)+".png");
            ASSERT_TRUE(file.good());
        }

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
// Helpers shared by all StarPropertiesEstimator tests
// ---------------------------------------------------------------------------
namespace
{
    // 10 well-separated positions: ≥80 px from any edge of 512×512, ≥120 px apart
    const std::vector<std::pair<double,double>> kStarPositions = {
        { 80.0,  80.0}, { 80.0, 200.0}, { 80.0, 320.0}, { 80.0, 440.0},
        {200.0, 140.0}, {200.0, 260.0}, {200.0, 380.0},
        {320.0,  80.0}, {320.0, 200.0}, {320.0, 320.0}
    };
    const double kBkg = 100.0;
    const double kTol = 0.05; // 5% relative tolerance

    // Run fitting for one model on all already-detected stars and return the updated table.
    DSL::FITStable runFit(DST::algorithms::starFinder& finder,
                          DST::Stars::evalType evalType)
    {
        const auto& detected = finder.detectedStars();
        const auto& ids = detected.getColumn("ID")->values<uint64_t>();
        for (size_t i = 0; i < detected.nrows(); ++i)
            finder.EvalStellarProperties(ids[i], evalType);
        return finder.detectedStars();
    }
} // anonymous namespace

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_GAUSS)
{
    // fits = [A, x0, y0, sigma, background]
    const double kA     = 5000.0;
    const double kSigma = 2.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::gauss2D(x, y, kA, x0, y0, kSigma);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_GAUSS");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::GAUSS);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)     / kA,     kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kSigma) / kSigma, kTol) << "star " << i << " sigma";
        EXPECT_LT(std::abs(fits[i][4] - kBkg)   / kBkg,   kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif

}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_iGAUSS)
{
    // fits = [A, x0, y0, sigma, background]
    const double kA     = 5000.0;
    const double kSigma = 2.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::igauss2D(x, y, kA, x0, y0, kSigma, 1.0, 1.0);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_iGAUSS");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::iGAUSS);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)     / kA,     kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kSigma) / kSigma, kTol) << "star " << i << " sigma";
        EXPECT_LT(std::abs(fits[i][4] - kBkg)   / kBkg,   kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif

}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_GAUSS2D)
{
    // fits = [A, x0, y0, sigma_x, sigma_y, theta, background]
    // Symmetric (sigma_x == sigma_y) to avoid theta degeneracy — theta not checked.
    const double kA      = 5000.0;
    const double kSigmaX = 2.0;
    const double kSigmaY = 2.0;
    const double kTheta  = DST::Math::MathCore::Pi() / 3.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::gauss2D(x, y, kA, x0, y0, kSigmaX, kSigmaY, kTheta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_GAUSS2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::GAUSS2D);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)      / kA,      kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kSigmaX) / kSigmaX, kTol) << "star " << i << " sigma_x";
        EXPECT_LT(std::abs(fits[i][4] - kSigmaY) / kSigmaY, kTol) << "star " << i << " sigma_y";
        EXPECT_LT(std::abs(fits[i][6] - kBkg)    / kBkg,    kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif

}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_iGAUSS2D)
{
    // fits = [A, x0, y0, sigma_x, sigma_y, theta, background]
    const double kA      = 5000.0;
    const double kSigmaX = 2.0;
    const double kSigmaY = 2.0;
    const double kTheta  = 0.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::igauss2D(x, y, kA, x0, y0, kSigmaX, kSigmaY, kTheta, 1.0, 1.0);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_iGAUSS2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::iGAUSS2D);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)      / kA,      kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kSigmaX) / kSigmaX, kTol) << "star " << i << " sigma_x";
        EXPECT_LT(std::abs(fits[i][4] - kSigmaY) / kSigmaY, kTol) << "star " << i << " sigma_y";
        EXPECT_LT(std::abs(fits[i][6] - kBkg)    / kBkg,    kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_MOFFAT)
{
    // fits = [A, x0, y0, alpha, beta, background]
    const double kA     = 5000.0;
    const double kAlpha = 2.0;
    const double kBeta  = 3.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::moffat2D(x, y, kA, x0, y0, kAlpha, kBeta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_MOFFAT");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::MOFFAT);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)     / kA,     kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kAlpha) / kAlpha, kTol) << "star " << i << " alpha";
        EXPECT_LT(std::abs(fits[i][4] - kBeta)  / kBeta,  kTol) << "star " << i << " beta";
        EXPECT_LT(std::abs(fits[i][5] - kBkg)   / kBkg,   kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_MOFFAT2D)
{
    // fits = [A, x0, y0, alpha, beta, gamma, theta, background]
    // gamma=1 (circular) to avoid theta degeneracy — theta not checked.
    const double kA     = 5000.0;
    const double kAlpha = 2.0;
    const double kBeta  = 3.0;
    const double kGamma = 1.0;
    const double kTheta = 0.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::moffat2D(x, y, kA, x0, y0, kAlpha, kBeta, kGamma, kTheta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_MOFFAT2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::MOFFAT2D);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)     / kA,     kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kAlpha) / kAlpha, kTol) << "star " << i << " alpha";
        EXPECT_LT(std::abs(fits[i][4] - kBeta)  / kBeta,  kTol) << "star " << i << " beta";
        EXPECT_LT(std::abs(fits[i][5] - kGamma) / kGamma, kTol) << "star " << i << " gamma";
        EXPECT_LT(std::abs(fits[i][7] - kBkg)   / kBkg,   kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_LORENTZIAN)
{
    // fits = [A, x0, y0, gamma, background]
    const double kA     = 5000.0;
    const double kGamma = 2.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::lorentzian2D(x, y, kA, x0, y0, kGamma);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_LORENTZIAN");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::LORENTZIAN);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)     / kA,     kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kGamma) / kGamma, kTol) << "star " << i << " gamma";
        EXPECT_LT(std::abs(fits[i][4] - kBkg)   / kBkg,   kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_LORENTZIAN2D)
{
    // fits = [A, x0, y0, gamma, q, theta, background]
    // q=1 (circular) to avoid theta degeneracy — q and theta not checked.
    const double kA      = 5000.0;
    const double kGamma2 = 2.0;
    const double kQ      = 1.0;
    const double kTheta  = 0.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::lorentzian2D(x, y, kA, x0, y0, kGamma2, kQ, kTheta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_LORENTZIAN2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::LORENTZIAN);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)      / kA,      kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kGamma2) / kGamma2, kTol) << "star " << i << " gamma";
        EXPECT_LT(std::abs(fits[i][4] - kBkg)    / kBkg,    kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_VOIGT)
{
    // fits = [A, x0, y0, fG, fL, background]
    const double kA  = 5000.0;
    const double kFG = 1.5;
    const double kFL = 1.5;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::pseudo_voigt2D(x, y, kA, x0, y0, kFG, kFL);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_VOIGT");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::VOIGT);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)   / kA,   kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kFG)  / kFG,  kTol) << "star " << i << " fG";
        EXPECT_LT(std::abs(fits[i][4] - kFL)  / kFL,  kTol) << "star " << i << " fL";
        EXPECT_LT(std::abs(fits[i][5] - kBkg) / kBkg, kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_VOIGT2D)
{
    // fits = [A, x0, y0, fG, fL, q, theta, background]
    // q=1 (circular) to avoid theta degeneracy — q and theta not checked.
    const double kA     = 5000.0;
    const double kFG2   = 1.5;
    const double kFL2   = 1.5;
    const double kQ     = 1.0;
    const double kTheta = 0.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::pseudo_voigt2D(x, y, kA, x0, y0, kFG2, kFL2, kQ, kTheta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_VOIGT2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::VOIGT2D);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)   / kA,   kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kFG2) / kFG2, kTol) << "star " << i << " fG";
        EXPECT_LT(std::abs(fits[i][4] - kFL2) / kFL2, kTol) << "star " << i << " fL";
        EXPECT_LT(std::abs(fits[i][7] - kBkg) / kBkg, kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_SERSIC)
{
    // fits = [A, x0, y0, re, n, background]
    const double kA  = 5000.0;
    const double kRe = 2.0;
    const double kN  = 1.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::sersic2D(x, y, kA, x0, y0, kRe, kN);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_SERSIC");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::SERSIC);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)   / kA,   kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kRe)  / kRe,  kTol) << "star " << i << " re";
        EXPECT_LT(std::abs(fits[i][4] - kN)   / kN,   kTol) << "star " << i << " n";
        EXPECT_LT(std::abs(fits[i][5] - kBkg) / kBkg, kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_SERSIC2D)
{
    // fits = [A, x0, y0, re, n, q, theta, background]
    // q=1 (circular) to avoid theta degeneracy — q and theta not checked.
    const double kA     = 5000.0;
    const double kRe2   = 2.0;
    const double kN2    = 1.0;
    const double kQ     = 1.0;
    const double kTheta = 0.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::sersic2D(x, y, kA, x0, y0, kRe2, kN2, kQ, kTheta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_SERSIC2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::SERSIC2D);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)   / kA,   kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kRe2) / kRe2, kTol) << "star " << i << " re";
        EXPECT_LT(std::abs(fits[i][4] - kN2)  / kN2,  kTol) << "star " << i << " n";
        EXPECT_LT(std::abs(fits[i][7] - kBkg) / kBkg, kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_KING)
{
    // fits = [A, x0, y0, rc, rt, background]
    const double kA  = 5000.0;
    const double kRc = 1.5;
    const double kRt = 6.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::king2D(x, y, kA, x0, y0, kRc, kRt);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_KING");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::KING);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)   / kA,   kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kRc)  / kRc,  kTol) << "star " << i << " rc";
        EXPECT_LT(std::abs(fits[i][4] - kRt)  / kRt,  kTol) << "star " << i << " rt";
        EXPECT_LT(std::abs(fits[i][5] - kBkg) / kBkg, kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_KING2D)
{
    // fits = [A, x0, y0, rc, rt, q, theta, background]
    // q=1 (circular) to avoid theta degeneracy — q and theta not checked.
    const double kA     = 5000.0;
    const double kRc2   = 1.5;
    const double kRt2   = 6.0;
    const double kQ     = 1.0;
    const double kTheta = 0.0;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::king2D(x, y, kA, x0, y0, kRc2, kRt2, kQ, kTheta);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_KING2D");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::KING2D);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)   / kA,   kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kRc2) / kRc2, kTol) << "star " << i << " rc";
        EXPECT_LT(std::abs(fits[i][4] - kRt2) / kRt2, kTol) << "star " << i << " rt";
        EXPECT_LT(std::abs(fits[i][7] - kBkg) / kBkg, kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }
    
#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
TEST(StarFinderTest, StarPropertiesEstimator_AIRY)
{
    // fits = [A, x0, y0, a, background]
    const double kA      = 5000.0;
    const double kAscale = 1.5;

    auto img = TestUtils::makeModelImage(512, 512, kBkg, kStarPositions,
        [&](double x, double y, double x0, double y0) {
            return DST::Math::function::airy2D(x, y, kA, x0, y0, kAscale);
        });

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/StarPropertiesEstimator_AIRY");
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    ASSERT_EQ(finder.detectedStars().nrows(), static_cast<size_t>(10));

    const auto& table = runFit(finder, Stars::evalType::AIRY);
    const auto& fits  = table.getColumn("fits")->values<FITSform::doubleVector>();
    for (size_t i = 0; i < table.nrows(); ++i)
    {
        EXPECT_LT(std::abs(fits[i][0] - kA)      / kA,      kTol) << "star " << i << " amplitude";
        EXPECT_LT(std::abs(fits[i][3] - kAscale) / kAscale, kTol) << "star " << i << " a";
        EXPECT_LT(std::abs(fits[i][4] - kBkg)    / kBkg,    kTol) << "star " << i << " background";
        EXPECT_GT(fits[i][1], 0.0) << "star " << i << " x0";
        EXPECT_GT(fits[i][2], 0.0) << "star " << i << " y0";
    }

#if !defined(Linuxx86_64) && !defined(Linuxaarch64)
    {
        int arg=0;
        QApplication app(arg, nullptr);

        QtFITSviewer *mw = new QtFITSviewer(250, 250, 600,600);
        mw->setWindowTitle(img->GetName().c_str());
        mw->show();

        QCustomPlot *customPlot = mw->GetPlot();
        customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );

        QtFITSviewer::QtDisplay(customPlot,img);
        EXPECT_NO_THROW(finder.Plot(customPlot->axisRect()));

        //app.exec();
    }
#endif
}

// ---------------------------------------------------------------------------
// connectedComponents — unit tests
// ---------------------------------------------------------------------------

TEST(ConnectedComponentsTest, EmptyImage_AllZero)
{
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{8, 8});
    for(size_t i = 0; i < img->Nelements(); ++i)
        img->SetPixelValue(0.0, i);

    auto labels = DST::Morphology::connectedComponents(img);
    ASSERT_NE(labels, nullptr);
    EXPECT_EQ(labels->Nelements(), img->Nelements());
    for(size_t i = 0; i < labels->Nelements(); ++i)
        EXPECT_EQ(labels->DoubleValueAtPixel(i), 0.0) << "pixel " << i;
}

TEST(ConnectedComponentsTest, SingleBlob_AllForeground)
{
    // 4×4 image fully foreground → all pixels share the same label
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{4, 4});
    for(size_t i = 0; i < img->Nelements(); ++i)
        img->SetPixelValue(1.0, i);

    auto labels = DST::Morphology::connectedComponents(img);
    double ref = labels->DoubleValueAtPixel(0);
    EXPECT_GT(ref, 0.0);
    for(size_t i = 0; i < labels->Nelements(); ++i)
        EXPECT_EQ(labels->DoubleValueAtPixel(i), ref) << "pixel " << i;
}

TEST(ConnectedComponentsTest, TwoSeparateBlobs)
{
    // 1×9 row: [1 1 1 0 0 0 1 1 1] → two blobs
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{9, 1});
    for(size_t i = 0; i < 9; ++i)
        img->SetPixelValue((i < 3 || i > 5) ? 1.0 : 0.0, i);

    auto labels = DST::Morphology::connectedComponents(img);
    double l0 = labels->DoubleValueAtPixel(0);
    double l6 = labels->DoubleValueAtPixel(6);
    EXPECT_GT(l0, 0.0);
    EXPECT_GT(l6, 0.0);
    EXPECT_NE(l0, l6) << "two disconnected blobs must receive different labels";
    for(size_t i = 0; i < 3; ++i)  EXPECT_EQ(labels->DoubleValueAtPixel(i), l0);
    for(size_t i = 3; i < 6; ++i)  EXPECT_EQ(labels->DoubleValueAtPixel(i), 0.0);
    for(size_t i = 6; i < 9; ++i)  EXPECT_EQ(labels->DoubleValueAtPixel(i), l6);
}

TEST(ConnectedComponentsTest, MaskedPixelsAreBackground)
{
    // 3×3 fully foreground but centre masked → centre label must be 0
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{3, 3});
    for(size_t i = 0; i < 9; ++i)
        img->SetPixelValue(1.0, i);
    img->MaskPixels({4}); // centre pixel (1,1)

    auto labels = DST::Morphology::connectedComponents(img);
    EXPECT_EQ(labels->DoubleValueAtPixel(4), 0.0) << "masked centre must have label 0";
}

TEST(ConnectedComponentsTest, NullInput_Throws)
{
    EXPECT_THROW(DST::Morphology::connectedComponents(nullptr),
                 DSL::FITSexception);
}

TEST(ConnectedComponentsTest, ThreeDInput_Throws)
{
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{4, 4, 4});
    EXPECT_THROW(DST::Morphology::connectedComponents(img),
                 DSL::FITSexception);
}

// ---------------------------------------------------------------------------
// skyEstimator — unit tests
// ---------------------------------------------------------------------------

TEST(SkyEstimatorTest, FlatBackground_Recovered)
{
    // Uniform sky = 50.0, no noise → every pixel of the sky map should be ~50
    const size_t W = 128, H = 128;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W, H});
    for(size_t i = 0; i < img->Nelements(); ++i)
        img->SetPixelValue(50.0, i);

    DST::Background::skyEstimator est;
    est.setTileSize(32);
    est.setSigmaClipIter(3);
    est.setSigmaClipThreshold(3.0);
    EXPECT_NO_THROW(est.estimate(img));

    auto skyMap = est.getSkyMap();
    ASSERT_NE(skyMap, nullptr);
    EXPECT_EQ(skyMap->Nelements(), img->Nelements());
    for(size_t i = 0; i < skyMap->Nelements(); ++i)
        EXPECT_NEAR(skyMap->DoubleValueAtPixel(i), 50.0, 1.0) << "pixel " << i;
}

TEST(SkyEstimatorTest, NoisyBackground_SkyReasonable)
{
    // Gaussian noise around sky=100 → per-tile sky estimate within 3 counts of truth
    const size_t W = 256, H = 256;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W, H});
    std::mt19937_64 rng(999);
    std::normal_distribution<double> noise(100.0, 5.0);
    for(size_t i = 0; i < img->Nelements(); ++i)
        img->SetPixelValue(noise(rng), i);

    DST::Background::skyEstimator est;
    est.setTileSize(64);
    est.setSigmaClipIter(3);
    est.setSigmaClipThreshold(3.0);
    est.estimate(img);

    for(const auto& ts : est.tileStats())
    {
        if(!ts.valid) continue;
        EXPECT_NEAR(ts.sky, 100.0, 3.0) << "tile sky too far from truth";
        EXPECT_GT(ts.rms, 0.0)          << "tile rms must be positive";
    }
}

TEST(SkyEstimatorTest, SourceMasked_SkyUnaffected)
{
    // Sky = 50.  Inject a bright source, mask it, check sky still ~50
    const size_t W = 128, H = 128;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W, H});
    for(size_t i = 0; i < img->Nelements(); ++i)
        img->SetPixelValue(50.0, i);

    std::vector<size_t> srcPixels;
    for(size_t y = 60; y < 68; ++y)
        for(size_t x = 60; x < 68; ++x)
        {
            const size_t idx = y * W + x;
            img->SetPixelValue(5000.0, idx);
            srcPixels.push_back(idx);
        }
    img->MaskPixels(srcPixels);

    DST::Background::skyEstimator est;
    est.setTileSize(32);
    est.estimate(img);

    auto skyMap = est.getSkyMap();
    for(size_t i = 0; i < skyMap->Nelements(); ++i)
    {
        if(img->Masked(i)) continue;
        EXPECT_NEAR(skyMap->DoubleValueAtPixel(i), 50.0, 5.0) << "pixel " << i;
    }
}

TEST(SkyEstimatorTest, RMSMapPositive)
{
    // RMS map values must always be strictly positive
    const size_t W = 64, H = 64;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W, H});
    std::mt19937_64 rng(42);
    std::normal_distribution<double> noise(200.0, 10.0);
    for(size_t i = 0; i < img->Nelements(); ++i)
        img->SetPixelValue(noise(rng), i);

    DST::Background::skyEstimator est;
    est.setTileSize(16);
    est.estimate(img);

    auto rmsMap = est.getRMSMap();
    ASSERT_NE(rmsMap, nullptr);
    for(size_t i = 0; i < rmsMap->Nelements(); ++i)
        EXPECT_GT(rmsMap->DoubleValueAtPixel(i), 0.0) << "rms must be > 0 at pixel " << i;
}

TEST(SkyEstimatorTest, NullInput_Throws)
{
    DST::Background::skyEstimator est;
    EXPECT_THROW(est.estimate(nullptr), DSL::FITSexception);
}

TEST(SkyEstimatorTest, ThreeDInput_Throws)
{
    DST::Background::skyEstimator est;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{8, 8, 8});
    EXPECT_THROW(est.estimate(img), DSL::FITSexception);
}

TEST(SkyEstimatorTest, GetSkyMapBeforeEstimate_Throws)
{
    DST::Background::skyEstimator est;
    EXPECT_THROW(est.getSkyMap(), std::runtime_error);
    EXPECT_THROW(est.getRMSMap(), std::runtime_error);
}

// ---------------------------------------------------------------------------
// NoiseChisel detection — integration tests (CPU / double path)
// ---------------------------------------------------------------------------

TEST(StarFinderTest, FindStarsByNoiseChisel_DefaultParams)
{
    std::vector<DST::Stars> true_catalog;
    const unsigned int seed = 12345;
    auto img = TestUtils::generateStarField(
        true_catalog,
        256, 256,  // smaller than MORPHO tests for speed
        30,
        1.5,       // sigma — slightly broader (extended-source regime)
        50.0, 800.0,
        20.0, 5.0,
        seed);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/FindStarsByNoiseChisel_DefaultParams");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setThresholdValue(0.85); // NoiseChisel is more aggressive by default, so increase threshold to avoid huge number of detections

    EXPECT_NO_THROW(
        finder.findStars(DST::algorithms::starFinder::detectionMethod::NOISECHISEL, 0.1)
    );

    const auto& detected = finder.detectedStars();
    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));

    const double rate = static_cast<double>(detected.nrows())
                      / static_cast<double>(true_catalog.size()) * 100.0;
    std::cout << "[NOISECHISEL default] " << detected.nrows() << "/"
              << true_catalog.size() << " (" << rate << "%)" << std::endl;

    EXPECT_GT(rate,  30.0); // conservative floor: NoiseChisel prioritises purity
    EXPECT_LT(rate, 200.0);
}

TEST(StarFinderTest, FindStarsByNoiseChisel_TunedParams)
{
    std::vector<DST::Stars> true_catalog;
    auto img = TestUtils::generateStarField(
        true_catalog, 256, 256, 30, 1.5, 50.0, 800.0, 20.0, 5.0, 42);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/FindStarsByNoiseChisel_TunedParams");
    verbose = verboseLevel::VERBOSE_NONE;

    finder.setNC_convSigma(1.5);
    finder.setNC_convRadius(3);
    finder.setNC_detectionQuantile(0.65);
    finder.setNC_tileSize(32);
    finder.setNC_erosionRadius(1);
    finder.setNC_snrThreshold(2.0);
    finder.setNC_dilationRadius(3);
    finder.setNC_skyIterations(3);
    finder.setNC_skyClipSigma(3.0);
    finder.setNC_reDetect(true);
    finder.setThresholdValue(0.945);

    EXPECT_NO_THROW(
        finder.findStars(DST::algorithms::starFinder::detectionMethod::NOISECHISEL, 0.1)
    );

    const auto& detected = finder.detectedStars();
    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));

    const double rate = static_cast<double>(detected.nrows())
                      / static_cast<double>(true_catalog.size()) * 100.0;
    std::cout << "[NOISECHISEL tuned]   " << detected.nrows() << "/"
              << true_catalog.size() << " (" << rate << "%)" << std::endl;

    EXPECT_GT(rate,  30.0);
    EXPECT_LT(rate, 200.0);
}

TEST(StarFinderTest, FindStarsByNoiseChisel_NoReDetect)
{
    // reDetect=false still produces detections and does not throw
    std::vector<DST::Stars> true_catalog;
    auto img = TestUtils::generateStarField(
        true_catalog, 256, 256, 20, 1.5, 50.0, 800.0, 20.0, 5.0, 7777);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/FindStarsByNoiseChisel_NoReDetect");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setNC_reDetect(false);
    finder.setThresholdValue(0.85);

    EXPECT_NO_THROW(
        finder.findStars(DST::algorithms::starFinder::detectionMethod::NOISECHISEL, 0.1)
    );
    EXPECT_GT(finder.detectedStars().nrows(), static_cast<size_t>(0));
}

TEST(StarFinderTest, FindStarsByNoiseChisel_SkySubtractedHasLowerMean)
{
    // The sky-subtracted intermediate FITS must have a lower mean than the raw input
    std::vector<DST::Stars> true_catalog;
    auto img = TestUtils::generateStarField(
        true_catalog, 128, 128, 10, 1.5, 50.0, 500.0, 30.0, 5.0, 2024);

    const double inputMean = [&]{
        double s = 0.0;
        for(size_t i = 0; i < img->Nelements(); ++i)
            s += img->DoubleValueAtPixel(i);
        return s / static_cast<double>(img->Nelements());
    }();

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output/SkySubtractedHasLowerMean");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setNC_tileSize(32);
    finder.setNC_skyIterations(2);
    finder.setThresholdValue(0.85);

    ASSERT_NO_THROW(
        finder.findStars(DST::algorithms::starFinder::detectionMethod::NOISECHISEL, 0.1)
    );

    const std::string skyPath = "./build/test_output/SkySubtractedHasLowerMean/nc_skysubtracted_layer_0.fits";
    std::ifstream f(skyPath);
    if(f.good())
    {
        f.close();
        FITSmanager fm(skyPath);
        std::shared_ptr<FITScube> skySubt = fm.GetPrimary();
        //skySubt->Read(skyPath);

        double subtMean = 0.0;
        size_t nValid   = 0;
        for(size_t i = 0; i < skySubt->Nelements(); ++i)
        {
            if(skySubt->Masked(i)) continue;
            subtMean += skySubt->DoubleValueAtPixel(i);
            ++nValid;
        }
        if(nValid > 0) subtMean /= static_cast<double>(nValid);

        std::cout << "[NOISECHISEL sky] input mean=" << inputMean
                  << "  sky-subtracted mean=" << subtMean << std::endl;
        EXPECT_LT(subtMean, inputMean)
            << "sky-subtracted image must have a lower mean than the raw input";
    }
    else
    {
        std::cout << "[NOISECHISEL sky] output file not found; skipping sky check" << std::endl;
    }
}

// ---------------------------------------------------------------------------
// FITS-LDAC export: artificial stellar field
// ---------------------------------------------------------------------------
TEST(starFinderTest, toFITS_LDAC_ArtificialField)
{
    std::vector<DST::Stars> true_catalog;
    auto sky = TestUtils::generateStarField(
        true_catalog, 128, 128, 10, 1.5, 50.0, 500.0, 30.0, 5.0, 2024);

    // 2. Detect stars (findStars runs MOMENT eval internally)
    DST::algorithms::starFinder finder(sky);
    finder.setThresholdValue(0.97);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);
    
    ASSERT_GT(finder.detectedStars().nrows(), size_t(0))
        << "No stars detected in artificial field";

    // 3. Write FITS-LDAC — re-evaluate with MOMENT (ensures all properties are set)
    const std::string outFile = "./build/test_output/artificial_field.ldac";
    ASSERT_NO_THROW(finder.toFITS_LDAC(outFile, DST::Stars::evalType::MOMENT));

    // 4. Validate file structure
    DSL::FITSmanager fm(outFile);
    EXPECT_EQ(fm.NumberOfHeader(), 3)   // primary + LDAC_IMHEAD + LDAC_OBJECTS
        << "Expected 3 HDUs in FITS-LDAC file";

    // HDU 2: LDAC_IMHEAD
    auto imhead = fm.GetTableAtIndex(2);
    ASSERT_NE(imhead, nullptr);
    EXPECT_EQ(imhead->GetName(), "LDAC_IMHEAD");
    EXPECT_EQ(imhead->nrows(),   size_t(1));
    EXPECT_EQ(imhead->ncols(),   size_t(1));

    // HDU 3: LDAC_OBJECTS
    auto objTbl = fm.GetTableAtIndex(3);
    ASSERT_NE(objTbl, nullptr);
    EXPECT_EQ(objTbl->GetName(), "LDAC_OBJECTS");
    EXPECT_GT(objTbl->nrows(),   size_t(0));

    // All SCAMP-required columns must be present
    for (const char* name : {"XWIN_IMAGE", "YWIN_IMAGE",
                              "ERRAWIN_IMAGE", "ERRBWIN_IMAGE", "ERRTHETAWIN_IMAGE",
                              "FLUX_AUTO", "FLUXERR_AUTO",
                              "MAG_AUTO", "MAGERR_AUTO", "FLAGS"})
        EXPECT_NO_THROW(objTbl->getColumn(name)) << "Missing SCAMP column: " << name;

    // Pixel coordinates must be positive (1-based FITS convention)
    const std::vector<float>& xVals = objTbl->column<float>("XWIN_IMAGE").data();
    for (size_t r = 0; r < xVals.size(); ++r)
        EXPECT_GT(xVals[r], 0.f) << "XWIN_IMAGE must be > 0 (1-based) at row " << r;
}
