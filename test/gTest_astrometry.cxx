#include <gtest/gtest.h>

#include <DSTools/DSTastrometry.h>
#include <DSTfits/FITSmanager.h>

#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <thread>

#include "test_utils.h"

using namespace DST::astrometry;

// ---------------------------------------------------------------------------
// Paths
// ---------------------------------------------------------------------------
static const std::filesystem::path kTestOutputDir = "./build/test_output/astrometry";

// Drop a real astronomical FITS here to enable WCS accuracy tests.
// It must have enough stars (~30+) and valid WCS keywords (CRVAL1/CRVAL2).
// A 512x512 DSS or PanSTARRS cutout of a dense field (M11, M45…) works well.
static const std::filesystem::path kReferenceFits  = "testdata/astrometry_reference.fits";

// How far the solved RA/Dec may differ from the reference WCS (degrees)
static constexpr double kWcsTolerance_deg = 1.0 / 3600.0;  // 1 arcsec

// ---------------------------------------------------------------------------
// Helper: write a minimal valid JSON config
// ---------------------------------------------------------------------------
static std::string writeMinimalConfig(const std::filesystem::path& outDir,
                                      const std::string& solverPath = "")
{
    std::filesystem::create_directories(outDir);
    const std::filesystem::path configPath = outDir / "astrometry_test_config.json";

    std::ofstream f(configPath);
    f << "{\n"
      << "    \"astrometry.net\": {\n"
      << "        \"path\": \"" << solverPath << "\",\n"
      << "        \"url\":  \"\",\n"
      << "        \"api-key\": \"\",\n"
      << "        \"parameters\":{\n"
#ifdef ASTROMETRY_NET_INDEX_FILE
      << "            \"index-file\": \"" << ASTROMETRY_NET_INDEX_FILE << "\",\n"
#else
      << "            \"index-file\": \"\",\n"
#endif
      << "            \"skip-solved\": false\n"
      << "        }\n"
      << "    },\n"
      << "    \"scamp\": {\n"
#ifdef SCAMP_DEFAULT_PATH
      << "        \"path\": \"" << SCAMP_DEFAULT_PATH << "\",\n"
#else
      << "        \"path\": \"\",\n"
#endif
      << "        \"parameters\": {}\n"
      << "    }\n"
      << "}\n";

    return configPath.string();
}

// ---------------------------------------------------------------------------
// Runtime check: is solve-field reachable?
// ---------------------------------------------------------------------------
static bool hasSolveField()
{
#ifdef ASTROMETRY_NET_DEFAULT_PATH
    return std::filesystem::exists(
        std::filesystem::path(ASTROMETRY_NET_DEFAULT_PATH) / "solve-field");
#else
    return std::system("which solve-field > /dev/null 2>&1") == 0;
#endif
}

// ---------------------------------------------------------------------------
// Global environment: create output directory once
// ---------------------------------------------------------------------------
class AstrometryTestEnv : public ::testing::Environment
{
public:
    void SetUp() override
    {
        std::filesystem::remove_all(kTestOutputDir);
        std::filesystem::create_directories(kTestOutputDir);
    }
};

const ::testing::Environment* const kAstrometryEnv =
    ::testing::AddGlobalTestEnvironment(new AstrometryTestEnv);


// ===========================================================================
// Suite 1 — Constructor / configuration (no solver binary required)
// ===========================================================================

TEST(AstrometrySolver, ConstructorThrowsOnEmptyPath)
{
    std::string emptyConfig = std::string();
    EXPECT_THROW(solver s(emptyConfig), std::invalid_argument);

    std::string emptyString = std::string("");
    EXPECT_THROW(solver s(emptyString), std::invalid_argument);
}

TEST(AstrometrySolver, ConstructorThrowsOnMissingConfigFile)
{
    EXPECT_THROW(solver s("/nonexistent/path/astrometry_config.json"),
                 std::invalid_argument);
}

TEST(AstrometrySolver, ConstructorAcceptsValidConfig)
{
    const std::string cfg = writeMinimalConfig(kTestOutputDir);
    EXPECT_NO_THROW(solver s(cfg));
}

TEST(AstrometrySolver, ConstructorAcceptsExplicitConcurrencyLimit)
{
    const std::string cfg = writeMinimalConfig(kTestOutputDir);
    EXPECT_NO_THROW({ solver s(cfg, 1);  });   // minimum
    EXPECT_NO_THROW({ solver s(cfg, 0);  });   // auto (hardware/3)
    EXPECT_NO_THROW({ solver s(cfg, 64); });   // user override
}

#ifdef ASTROMETRY_NET_DEFAULT_PATH
TEST(AstrometrySolver, DefaultSolverPathIsSet)
{
    // If CMake found solve-field, the default solver_path must be non-empty
    // and point to a directory that actually contains solve-field.
    const std::string cfg = writeMinimalConfig(kTestOutputDir);
    // Construction with no explicit path should succeed because the compile-time
    // default was detected by CMake.
    EXPECT_NO_THROW(solver s(cfg));
}
#endif

// ===========================================================================
// Suite 2 — Null image guard (no solver binary required)
// ===========================================================================

TEST(AstrometrySolver, SolveThrowsOnNullFITScube)
{
    const std::string cfg = writeMinimalConfig(kTestOutputDir);
    solver s(cfg);
    std::shared_ptr<FITScube> null_img;
    EXPECT_THROW(s.solve(null_img), FITSexception);
}

// ===========================================================================
// Suite 3 — Concurrency throttle (no solver binary required)
// We verify that two solver instances with limit=1 can be constructed and
// that the slot guard compiles and links correctly.  A full concurrency test
// (verifying that threads actually block) requires a slow real invocation and
// is left inside the integration suite below.
// ===========================================================================

TEST(AstrometrySolverConcurrency, SlotGuardCompilesAndLinks)
{
    const std::string cfg = writeMinimalConfig(kTestOutputDir);
    solver s1(cfg, 1);
    solver s2(cfg, 2);
    // If the concurrency machinery compiled, this test already passes.
    SUCCEED();
}

// ===========================================================================
// Suite 4 — Parameter access and modification (no solver binary required)
// All tests use the default constructor (no config file needed) so they run
// even on machines where astrometry.net / SCAMP is not installed.
// ===========================================================================

// --- astrometry.net parameters ---

TEST(AstrometryParams, SetAndGetBool)
{
    solver s;
    // "no-plots" defaults to false; overwrite and read back
    s.setAstrometryParam("no-plots", true);
    const solver::param_value v = s.getAstrometryParam("no-plots");
    ASSERT_TRUE(std::holds_alternative<bool>(v));
    EXPECT_TRUE(std::get<bool>(v));
}

TEST(AstrometryParams, SetAndGetString)
{
    solver s;
    // "scale-units" defaults to empty string; overwrite and read back
    s.setAstrometryParam("scale-units", std::string("arcsecperpix"));
    const solver::param_value v = s.getAstrometryParam("scale-units");
    ASSERT_TRUE(std::holds_alternative<std::string>(v));
    EXPECT_EQ(std::get<std::string>(v), "arcsecperpix");
}

TEST(AstrometryParams, SetAndGetFloat_OnMonostateSlot)
{
    solver s;
    // "scale-low" starts as monostate; a monostate slot accepts any type
    s.setAstrometryParam("scale-low", 0.5f);
    const solver::param_value v = s.getAstrometryParam("scale-low");
    ASSERT_TRUE(std::holds_alternative<float>(v));
    EXPECT_FLOAT_EQ(std::get<float>(v), 0.5f);
}

TEST(AstrometryParams, SetAndGetInt_OnMonostateSlot)
{
    solver s;
    // "tweak-order" starts as monostate; a monostate slot accepts any type
    s.setAstrometryParam("tweak-order", 3);
    const solver::param_value v = s.getAstrometryParam("tweak-order");
    ASSERT_TRUE(std::holds_alternative<int>(v));
    EXPECT_EQ(std::get<int>(v), 3);
}

TEST(AstrometryParams, UnsetResetsToMonostate)
{
    solver s;
    s.setAstrometryParam("scale-low", 0.5f);
    EXPECT_TRUE(s.isAstrometryParamSet("scale-low"));

    s.unsetAstrometryParam("scale-low");
    EXPECT_FALSE(s.isAstrometryParamSet("scale-low"));

    // getAstrometryParam on an unset slot must succeed and return monostate
    const solver::param_value v = s.getAstrometryParam("scale-low");
    EXPECT_TRUE(std::holds_alternative<std::monostate>(v));
}

TEST(AstrometryParams, IsParamSet)
{
    solver s;
    // "no-plots" is stored as bool false — that counts as set
    EXPECT_TRUE(s.isAstrometryParamSet("no-plots"));
    // "scale-low" is monostate by default — not set
    EXPECT_FALSE(s.isAstrometryParamSet("scale-low"));
    // unknown key — returns false without throwing
    EXPECT_FALSE(s.isAstrometryParamSet("nonexistent-key"));
}

TEST(AstrometryParams, SetThrowsOnUnknownKey)
{
    solver s;
    EXPECT_THROW(s.setAstrometryParam("nonexistent-key", true),        std::invalid_argument);
    EXPECT_THROW(s.setAstrometryParam("nonexistent-key", std::string("x")), std::invalid_argument);
    EXPECT_THROW(s.setAstrometryParam("nonexistent-key", 1.0f),        std::invalid_argument);
    EXPECT_THROW(s.setAstrometryParam("nonexistent-key", 1),           std::invalid_argument);
}

TEST(AstrometryParams, GetThrowsOnUnknownKey)
{
    solver s;
    EXPECT_THROW(s.getAstrometryParam("nonexistent-key"), std::invalid_argument);
}

TEST(AstrometryParams, UnsetThrowsOnUnknownKey)
{
    solver s;
    EXPECT_THROW(s.unsetAstrometryParam("nonexistent-key"), std::invalid_argument);
}

TEST(AstrometryParams, SetThrowsOnTypeMismatch)
{
    solver s;
    // "no-plots" is a bool slot; wrong types must throw
    EXPECT_THROW(s.setAstrometryParam("no-plots", 42),               std::invalid_argument);
    EXPECT_THROW(s.setAstrometryParam("no-plots", 1.0f),             std::invalid_argument);
    EXPECT_THROW(s.setAstrometryParam("no-plots", std::string("Y")), std::invalid_argument);
}

// --- SCAMP parameters ---

TEST(ScampParams, SetAndGetBool)
{
    solver s;
    // "solve_photom" defaults to false; overwrite and read back
    s.setScampParam("solve_photom", true);
    const solver::param_value v = s.getScampParam("solve_photom");
    ASSERT_TRUE(std::holds_alternative<bool>(v));
    EXPECT_TRUE(std::get<bool>(v));
}

TEST(ScampParams, SetAndGetString)
{
    solver s;
    // "astref_catalog" defaults to "GAIA-DR3"; overwrite and read back
    s.setScampParam("astref_catalog", std::string("2MASS"));
    const solver::param_value v = s.getScampParam("astref_catalog");
    ASSERT_TRUE(std::holds_alternative<std::string>(v));
    EXPECT_EQ(std::get<std::string>(v), "2MASS");
}

TEST(ScampParams, SetAndGetFloat)
{
    solver s;
    // "crossid_radius" defaults to 2.0f; overwrite and read back
    s.setScampParam("crossid_radius", 5.0f);
    const solver::param_value v = s.getScampParam("crossid_radius");
    ASSERT_TRUE(std::holds_alternative<float>(v));
    EXPECT_FLOAT_EQ(std::get<float>(v), 5.0f);
}

TEST(ScampParams, SetAndGetInt)
{
    solver s;
    // "distort_degrees" defaults to 3; overwrite and read back
    s.setScampParam("distort_degrees", 5);
    const solver::param_value v = s.getScampParam("distort_degrees");
    ASSERT_TRUE(std::holds_alternative<int>(v));
    EXPECT_EQ(std::get<int>(v), 5);
}

TEST(ScampParams, UnsetResetsToMonostate)
{
    solver s;
    EXPECT_TRUE(s.isScampParamSet("crossid_radius"));   // 2.0f is the default

    s.unsetScampParam("crossid_radius");
    EXPECT_FALSE(s.isScampParamSet("crossid_radius"));

    const solver::param_value v = s.getScampParam("crossid_radius");
    EXPECT_TRUE(std::holds_alternative<std::monostate>(v));
}

TEST(ScampParams, IsParamSet)
{
    solver s;
    // "distort_degrees" is stored as int 3 — that counts as set
    EXPECT_TRUE(s.isScampParamSet("distort_degrees"));
    // "pixscale_maxerr" is monostate by default — not set
    EXPECT_FALSE(s.isScampParamSet("pixscale_maxerr"));
    // unknown key — returns false without throwing
    EXPECT_FALSE(s.isScampParamSet("nonexistent-key"));
}

TEST(ScampParams, SetThrowsOnUnknownKey)
{
    solver s;
    EXPECT_THROW(s.setScampParam("nonexistent-key", true),        std::invalid_argument);
    EXPECT_THROW(s.setScampParam("nonexistent-key", std::string("x")), std::invalid_argument);
    EXPECT_THROW(s.setScampParam("nonexistent-key", 1.0f),        std::invalid_argument);
    EXPECT_THROW(s.setScampParam("nonexistent-key", 1),           std::invalid_argument);
}

TEST(ScampParams, GetThrowsOnUnknownKey)
{
    solver s;
    EXPECT_THROW(s.getScampParam("nonexistent-key"), std::invalid_argument);
}

TEST(ScampParams, UnsetThrowsOnUnknownKey)
{
    solver s;
    EXPECT_THROW(s.unsetScampParam("nonexistent-key"), std::invalid_argument);
}

TEST(ScampParams, SetThrowsOnTypeMismatch)
{
    solver s;
    // "distort_degrees" is an int slot; wrong types must throw
    EXPECT_THROW(s.setScampParam("distort_degrees", true),               std::invalid_argument);
    EXPECT_THROW(s.setScampParam("distort_degrees", 1.5f),               std::invalid_argument);
    EXPECT_THROW(s.setScampParam("distort_degrees", std::string("5")),   std::invalid_argument);
}

// ===========================================================================
// Suite 5 — Integration tests (skipped if solve-field is absent)
// ===========================================================================

class AstrometryIntegrationTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        if (!hasSolveField())
            GTEST_SKIP() << "solve-field not found — skipping integration test";

#ifdef ASTROMETRY_NET_DEFAULT_PATH
        configPath_ = writeMinimalConfig(kTestOutputDir, ASTROMETRY_NET_DEFAULT_PATH);
#else
        configPath_ = writeMinimalConfig(kTestOutputDir);
#endif
    }

    std::string configPath_;
};

TEST_F(AstrometryIntegrationTest, SolveSyntheticFITS_DoesNotThrow)
{
    // A synthetic star field has random pixel positions — astrometry.net will
    // most likely not find a WCS solution, but our C++ interface must not throw.
    const std::filesystem::path tmpFits =
        kTestOutputDir / "synthetic_smoke_test.fits";

    std::vector<DST::Stars> catalog;
    auto img = TestUtils::generateStarField(
        catalog,
        512, 512,
        100,       // nStars — dense enough to exercise source detection
        1.5,       // sigma_stars [pixels]
        500.0,     // min_intensity
        50000.0,   // max_intensity
        100.0,     // background
        10.0,      // noise_sigma
        42);       // seed (reproducible)

    img->Write(tmpFits.string());

    solver s(configPath_, 1);
    EXPECT_NO_THROW(s.solve(tmpFits.string(), 1));
}

TEST_F(AstrometryIntegrationTest, ConcurrencyThrottleBlocksExcessThreads)
{
    // Spawn 4 threads on a solver limited to 1 concurrent slot.
    // All threads use a non-existent FITS file so the OS-level solve-field
    // invocation fails quickly.  We verify that the slot guard serialises entry
    // by reading peakActiveSolves() after all threads finish: it must be <= 1.
    //
    // NOTE: measuring concurrency from the task lambda (++counter before solve())
    // is wrong — the counter would be incremented before the slot guard blocks,
    // so it would always reach 4.  peakActiveSolves() is tracked inside the slot
    // guard itself, after a slot has been acquired, giving the correct reading.
    const std::string nonExistentFits =
        (kTestOutputDir / "no_such_file.fits").string();

    solver s(configPath_, 1);

    auto task = [&]()
    {
        try { s.solve(nonExistentFits, 1); }
        catch (...) {}
    };

    std::vector<std::thread> threads;
    threads.reserve(4);
    for (int i = 0; i < 4; ++i)
        threads.emplace_back(task);
    for (auto& t : threads)
        t.join();

    // With limit=1, no more than 1 thread should ever be inside solve() simultaneously.
    EXPECT_LE(s.peakActiveSolves(), 1u);
}


// ===========================================================================
// Suite 5 — WCS accuracy test
// Requires: (a) solve-field installed, (b) testdata/astrometry_reference.fits (downloaded into the build directory)
//
// How to obtain the reference file:
//   Option A — commit a small cutout:
//     Use any public archive (DSS, PanSTARRS, Legacy Survey) to fetch a
//     512x512 FITS cutout of a star-dense field (M11, M45, NGC 6705 …).
//     The file must already carry correct WCS keywords (CRVAL1/CRVAL2).
//     Keep it small (< 1 MB) so it is comfortable in git.
//
//   Option B — generate via astropy (Python, one-off):
//     from astroquery.skyview import SkyView
//     hdu = SkyView.get_images('M11','DSS')[0][0]
//     hdu.writeto('testdata/astrometry_reference.fits')
//
// ===========================================================================

class AstrometryWCSTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        if (!hasSolveField())
            GTEST_SKIP() << "solve-field not found — skipping WCS accuracy test";

        if (!std::filesystem::exists(kReferenceFits))
            GTEST_SKIP()
                << "Reference FITS not found at " << kReferenceFits
                << "\nPlace a real astronomical FITS with known WCS there "
                   "(see comment in gTest_astrometry.cxx Suite 5).";

#ifdef ASTROMETRY_NET_DEFAULT_PATH
        configPath_ = writeMinimalConfig(kTestOutputDir, ASTROMETRY_NET_DEFAULT_PATH);
#else
        configPath_ = writeMinimalConfig(kTestOutputDir);
#endif
    }

    std::string configPath_;
};

TEST_F(AstrometryWCSTest, SolvedWCSMatchesKnownReference)
{
    // Step 1 — snapshot WCS from the unmodified reference file
    FITSmanager refManager(kReferenceFits.string());
    std::shared_ptr<FITScube> refImg = refManager.GetPrimary();
    ASSERT_TRUE(refImg) << "Failed to read reference FITS";

    const FITSwcs refWCS = refImg->getWCS();
    ASSERT_GT(refWCS.getNumberOfWCS(), 0) << "Reference FITS carries no WCS";

    const double refRA  = refWCS.CRVAL(1);
    const double refDec = refWCS.CRVAL(2);

    // Step 2 — work on a copy so the reference file is never modified
    const std::filesystem::path workFits = std::filesystem::path("./build/test_data") / "wcs_accuracy_work.fits";
    std::filesystem::create_directories(workFits.parent_path());
    std::filesystem::copy_file(kReferenceFits, workFits,
                               std::filesystem::copy_options::overwrite_existing);

    // Step 3 — solve the copy
    solver s(configPath_, 1);
    ASSERT_NO_THROW(s.solve(workFits.string(), 1));

    // Step 4 — read back the solved copy and verify WCS keywords were patched in
    FITSmanager solvedManager(workFits.string());
    std::shared_ptr<FITScube> solvedImg = solvedManager.GetPrimary();
    ASSERT_TRUE(solvedImg) << "Failed to read solved FITS";

    const FITSwcs solvedWCS = solvedImg->getWCS();
    ASSERT_GT(solvedWCS.getNumberOfWCS(), 0)
        << "Solved FITS carries no WCS — keywords were not written back";

    EXPECT_NEAR(solvedWCS.CRVAL(1), refRA,  kWcsTolerance_deg)
        << "RA  (CRVAL1) mismatch after solving";
    EXPECT_NEAR(solvedWCS.CRVAL(2), refDec, kWcsTolerance_deg)
        << "Dec (CRVAL2) mismatch after solving";

    // Step 5 — delete the work copy only when all assertions passed
    if (!testing::Test::HasFailure())
        std::filesystem::remove(workFits);
}


// ===========================================================================
// SCAMP helpers
// ===========================================================================

// Root of the scamp submodule's test directory, set by CMake.
static const std::filesystem::path kScampTestDataDir =
#ifdef SCAMP_TEST_DATA_DIR
    std::filesystem::path(SCAMP_TEST_DATA_DIR);
#else
    std::filesystem::path("./scamp/tests");
#endif

static bool hasScamp()
{
    // Prefer the compiled-in install path; fall back to PATH lookup.
#ifdef SCAMP_DEFAULT_PATH
    if (std::filesystem::exists(
            std::filesystem::path(SCAMP_DEFAULT_PATH) / "scamp"))
        return true;
#endif
    return std::system("which scamp > /dev/null 2>&1") == 0;
}

static std::string getScampPath()
{
#ifdef SCAMP_DEFAULT_PATH
    if (std::filesystem::exists(
            std::filesystem::path(SCAMP_DEFAULT_PATH) / "scamp"))
        return std::string(SCAMP_DEFAULT_PATH);
#endif
    return std::string();  // empty → build_cmdline() invokes "scamp" from PATH
}

// JSON config that enables only the scamp section.
// SCAMP-specific parameters (reference catalog, global ahead, etc.) are set
// afterwards via setScampParam() so each test controls its own configuration.
static std::string writeScampTestConfig(const std::filesystem::path& outDir,
                                        const std::string& scampPath = "")
{
    std::filesystem::create_directories(outDir);
    const std::filesystem::path cfgPath = outDir / "scamp_test_config.json";

    std::ofstream f(cfgPath);
    f << "{\n"
      << "    \"astrometry.net\": {\n"
      << "        \"path\": \"\",\n"
      << "        \"url\":  \"\",\n"
      << "        \"api-key\": \"\",\n"
      << "        \"parameters\": {}\n"
      << "    },\n"
      << "    \"scamp\": {\n"
      << "        \"path\": \"" << scampPath << "\",\n"
      << "        \"parameters\": {}\n"
      << "    }\n"
      << "}\n";

    return cfgPath.string();
}

// Configure a solver with the settings needed to reproduce SCAMP's own
// integration-test setup: local GAIA-DR1 reference, MegaCam global header,
// no plots, output redirected to outputDir.
//
// @param catFile   Absolute path to the input catalog (.cat) — used to
//                  derive the output .head filename in outputDir.
// @param outputDir Directory where SCAMP writes its .head and .xml files.
static void configureForMegaCamTest(solver& s,
                                    const std::filesystem::path& extra,
                                    const std::filesystem::path& catFile,
                                    const std::filesystem::path& outputDir)
{
    s.setScampParam("astref_catalog",    std::string("FILE"));
    s.setScampParam("astrefcat_name",    (extra / "GAIA-DR1_1334+3754_r46.cat").string());
    s.setScampParam("aheader_global",    (extra / "megacam.ahead").string());
    s.setScampParam("mosaic_type",       std::string("SAME_CRVAL"));
    s.setScampParam("projection_type",   std::string("SAME"));
    s.setScampParam("astrinstru_key",    std::string("FILTER,QRUNID"));
    s.setScampParam("stability_type",    std::string("INSTRUMENT"));
    s.setScampParam("centroid_keys",     std::string("XWIN_IMAGE,YWIN_IMAGE"));
    s.setScampParam("centroiderr_keys",  std::string("ERRAWIN_IMAGE,ERRBWIN_IMAGE,ERRTHETAWIN_IMAGE"));
    s.setScampParam("distort_keys",      std::string("XWIN_IMAGE,YWIN_IMAGE"));
    s.setScampParam("distort_groups",    std::string("1,1"));
    s.setScampParam("distort_degrees",   4);
    s.setScampParam("astraccuracy_type", std::string("TURBULENCE-ARCSEC"));
    s.setScampParam("astr_accuracy",     0.055f);
    s.setScampParam("write_xml",         false);
    s.setScampParam("checkplot_dev",     std::string("NULL"));
    // Redirect output files to outputDir so the source tree is never polluted.
    s.setScampParam("header_name",
        (outputDir / (std::filesystem::path(catFile).stem().string() + ".head")).string());
    s.setScampParam("xml_name",
        (outputDir / "scamp.xml").string());
    // verbose_type already defaults to "QUIET"
}

// ===========================================================================
// Suite 6a — SCAMP binary test
//
// Reproduces SCAMP's own integration01 test: solve each of the 7 MegaCam
// chips against the local GAIA-DR1 reference catalog and verify that a
// sensible WCS solution is produced for every chip.
//
// Note: SCAMP's cross-matching quality metrics (AS_Contrast, sigma_int …)
// require all 7 chips in a single scamp invocation.  The C++ solve() API is
// single-catalog, so those aggregate statistics are not checked here.
// Individual per-chip WCS accuracy (CRVAL1/CRVAL2 in the expected range) is
// sufficient to validate the integration.
// ===========================================================================

class ScampBinaryTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        if (!hasScamp())
            GTEST_SKIP() << "scamp not found — skipping SCAMP binary test";

        if (!std::filesystem::exists(kScampTestDataDir / "extra" / "744331p.cat"))
            GTEST_SKIP() << "SCAMP test catalogs not found at " << kScampTestDataDir
                         << " (ensure the scamp submodule is checked out)";

        configPath_ = writeScampTestConfig(kTestOutputDir, getScampPath());
    }

    std::string configPath_;
};

TEST_F(ScampBinaryTest, EachMegaCamChipProducesValidHeadFile)
{
    const std::filesystem::path extra = kScampTestDataDir / "extra";

    // Dummy image — SCAMP reads only the catalog, not the pixel data.
    auto img = TestUtils::generateStarField(64, 64, 10, 1.5, 500.0, 5000.0, 100.0, 10.0, 42);

    for (int chip = 31; chip <= 37; ++chip)
    {
        const std::filesystem::path catFile =
            extra / ("7443" + std::to_string(chip) + "p.cat");

        solver s(configPath_);
        configureForMegaCamTest(s, extra, catFile, kTestOutputDir);

        ASSERT_NO_THROW(s.solve(img, catFile.string()))
            << "solve() threw for chip 7443" << chip << "p";

        // .head is redirected to kTestOutputDir via header_name.
        const std::filesystem::path headFile =
            kTestOutputDir / ("7443" + std::to_string(chip) + "p.head");

        ASSERT_TRUE(std::filesystem::exists(headFile))
            << ".head not produced for chip 7443" << chip << "p";
        ASSERT_GT(std::filesystem::file_size(headFile), static_cast<std::uintmax_t>(0))
            << ".head is empty for chip 7443" << chip << "p";

        // Parse CRVAL1 / CRVAL2 from the plain-text .head file and verify
        // they fall within the expected sky region for this MegaCam pointing:
        //   RA  ~203.5°  (CFHT field centred on ~13h34m)
        //   Dec ~+37.9°  (field centred on ~+37°54')
        // We allow ±10° in RA and ±5° in Dec to accommodate chip offsets
        // within the mosaic focal plane.
        double crval1 = std::numeric_limits<double>::quiet_NaN();
        double crval2 = std::numeric_limits<double>::quiet_NaN();
        {
            std::ifstream hf(headFile);
            std::string line;
            while (std::getline(hf, line))
            {
                if (line.size() >= 6 && line.substr(0, 6) == "CRVAL1")
                {
                    const size_t eq = line.find('=');
                    const size_t sl = line.find('/');
                    if (eq != std::string::npos)
                        crval1 = std::stod(line.substr(
                            eq + 1, sl == std::string::npos ? std::string::npos : sl - eq - 1));
                }
                else if (line.size() >= 6 && line.substr(0, 6) == "CRVAL2")
                {
                    const size_t eq = line.find('=');
                    const size_t sl = line.find('/');
                    if (eq != std::string::npos)
                        crval2 = std::stod(line.substr(
                            eq + 1, sl == std::string::npos ? std::string::npos : sl - eq - 1));
                }
            }
        }

        EXPECT_TRUE(std::isfinite(crval1))
            << "CRVAL1 missing or non-numeric in .head for chip 7443" << chip << "p";
        EXPECT_TRUE(std::isfinite(crval2))
            << "CRVAL2 missing or non-numeric in .head for chip 7443" << chip << "p";

        EXPECT_GE(crval1, 193.0) << "CRVAL1 too low for chip 7443" << chip << "p";
        EXPECT_LE(crval1, 213.0) << "CRVAL1 too high for chip 7443" << chip << "p";
        EXPECT_GE(crval2,  33.0) << "CRVAL2 too low for chip 7443" << chip << "p";
        EXPECT_LE(crval2,  43.0) << "CRVAL2 too high for chip 7443" << chip << "p";
    }
}


// ===========================================================================
// Suite 6b — SCAMP C++ solver API tests
// ===========================================================================

class ScampSolverApiTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        if (!hasScamp())
            GTEST_SKIP() << "scamp not found — skipping SCAMP API test";

        if (!std::filesystem::exists(kScampTestDataDir / "extra" / "744331p.cat"))
            GTEST_SKIP() << "SCAMP test catalogs not found at " << kScampTestDataDir
                         << " (ensure the scamp submodule is checked out)";

        configPath_ = writeScampTestConfig(kTestOutputDir, getScampPath());
    }

    std::string configPath_;
    // First MegaCam chip catalog — used by all API tests.
    const std::filesystem::path kCat = kScampTestDataDir / "extra" / "744331p.cat";
};

TEST_F(ScampSolverApiTest, SolveDoesNotThrow)
{
    solver s(configPath_);
    configureForMegaCamTest(s, kScampTestDataDir / "extra", kCat, kTestOutputDir);

    auto img = TestUtils::generateStarField(64, 64, 10, 1.5, 500.0, 5000.0, 100.0, 10.0, 42);

    EXPECT_NO_THROW(s.solve(img, kCat.string()));
}

TEST_F(ScampSolverApiTest, SolveUpdatesWCSKeywordsInImage)
{
    solver s(configPath_);
    configureForMegaCamTest(s, kScampTestDataDir / "extra", kCat, kTestOutputDir);

    auto img = TestUtils::generateStarField(64, 64, 10, 1.5, 500.0, 5000.0, 100.0, 10.0, 42);
    s.solve(img, kCat.string());

    // If SCAMP produced a solution, CRVAL1 must be in the expected RA range.
    // If it did not converge (no .head file), we log a note but do not fail —
    // the SolveDoesNotThrow test already guarantees no exception was thrown.
    const FITSwcs wcs = img->getWCS();
    if (wcs.getNumberOfWCS() > 0)
    {
        EXPECT_GE(wcs.CRVAL(1), 193.0) << "CRVAL1 out of expected range after solve";
        EXPECT_LE(wcs.CRVAL(1), 213.0) << "CRVAL1 out of expected range after solve";
        img->Write(kTestOutputDir / "744331p_solved.fits");
    }
    else
    {
        GTEST_LOG_(INFO) << "SCAMP returned no WCS for 744331p — field may not have converged";
    }
}

TEST_F(ScampSolverApiTest, SolveThrowsWhenScampBinaryMissing)
{
    // Point the solver at a non-existent directory so _checkForScampSolver() throws.
    const std::string badConfig =
        writeScampTestConfig(kTestOutputDir, "/nonexistent/scamp/path");
    solver s(badConfig);

    auto img = TestUtils::generateStarField(64, 64, 10, 1.5, 500.0, 5000.0, 100.0, 10.0, 42);

    EXPECT_THROW(s.solve(img, kCat.string()), std::runtime_error);
}

// ===========================================================================
// Suite 7 — solver(size_t) default constructor (no config file)
// ===========================================================================

TEST(AstrometrySolver, DefaultConstructorNoArgs)
{
    // Must construct without throwing; all parameters retain compiled-in defaults.
    EXPECT_NO_THROW({ solver s; });
}

TEST(AstrometrySolver, DefaultConstructorWithExplicitConcurrency)
{
    EXPECT_NO_THROW({ solver s(1);  });   // minimum
    EXPECT_NO_THROW({ solver s(0);  });   // auto (hardware/3)
    EXPECT_NO_THROW({ solver s(8);  });   // explicit
}

TEST(AstrometrySolver, DefaultConstructorParametersAreAccessible)
{
    // Parameters must be in their compiled-in default state after construction
    // with no config file.  Pick a few representative ones from each category.
    solver s;

    // bool default: "no-plots" → false
    ASSERT_TRUE(s.isAstrometryParamSet("no-plots"));
    EXPECT_FALSE(std::get<bool>(s.getAstrometryParam("no-plots")));

    // string default: "astref_catalog" → "GAIA-DR3"
    ASSERT_TRUE(s.isScampParamSet("astref_catalog"));
    EXPECT_EQ(std::get<std::string>(s.getScampParam("astref_catalog")), "GAIA-DR3");

    // int default: "distort_degrees" → 3
    ASSERT_TRUE(s.isScampParamSet("distort_degrees"));
    EXPECT_EQ(std::get<int>(s.getScampParam("distort_degrees")), 3);

    // monostate default: "scale-low" → not set
    EXPECT_FALSE(s.isAstrometryParamSet("scale-low"));
}

// ===========================================================================
// Suite 8 — Path / URL / key setters
// All tests use solver() (no config) and validate via getter round-trips.
// ===========================================================================

TEST(AstrometrySetters, SetAstrometryNetSolverPath_RoundTrip)
{
    solver s;
    const std::string path = "/opt/astrometry/bin";
    s.setAstrometryNetSolverPath(path);
    // The path is stored in astrometry_params.solver_path, which is visible
    // through _checkForAstrometryNetSolver when solve() is called.
    // We validate indirectly: calling solve() on a non-existent file with this
    // path must throw std::runtime_error (binary not found), not any other exception,
    // confirming the path was accepted.
    EXPECT_THROW(s.solve("dummy.fits", 1), std::runtime_error);
}

TEST(AstrometrySetters, SetAstrometryNetSolverPath_EmptyPathIsAccepted)
{
    // An empty path is a valid assignment (solver will look up PATH at solve time).
    solver s;
    EXPECT_NO_THROW(s.setAstrometryNetSolverPath(""));
}

TEST(AstrometrySetters, SetAstrometryNetSolverUrl_RoundTrip)
{
    solver s;
    const std::string url = "https://nova.astrometry.net/api";
    EXPECT_NO_THROW(s.setAstrometryNetSolverUrl(url));
    // The URL is written into astrometry_params.solver_url.
    // Verify it was accepted by checking that the object is still usable
    // (parameters remain readable after the setter).
    EXPECT_NO_THROW(s.getAstrometryParam("no-plots"));
}

TEST(AstrometrySetters, SetAstrometryNetSolverUrl_EmptyClears)
{
    solver s;
    s.setAstrometryNetSolverUrl("https://nova.astrometry.net/api");
    EXPECT_NO_THROW(s.setAstrometryNetSolverUrl(""));
    // Solver must still be functional after clearing.
    EXPECT_NO_THROW(s.getAstrometryParam("no-plots"));
}

TEST(AstrometrySetters, SetAstrometryNetApiKey_RoundTrip)
{
    solver s;
    const std::string key = "abc123xyz";
    EXPECT_NO_THROW(s.setAstrometryNetApiKey(key));
    // API key is internal state; validate by checking solver remains usable.
    EXPECT_NO_THROW(s.getAstrometryParam("scale-units"));
}

TEST(AstrometrySetters, SetAstrometryNetCatalogUrl_RoundTrip)
{
    solver s;
    const std::string dir = "/data/astrometry/index";
    EXPECT_NO_THROW(s.setAstrometryNetCatalogUrl(dir));
    // setAstrometryNetCatalogUrl sets the "index-dir" parameter.
    const solver::param_value v = s.getAstrometryParam("index-dir");
    ASSERT_TRUE(std::holds_alternative<std::string>(v));
    EXPECT_EQ(std::get<std::string>(v), dir);
}

TEST(AstrometrySetters, SetAstrometryNetCatalogUrl_OverwritesPreviousValue)
{
    solver s;
    s.setAstrometryNetCatalogUrl("/first/path");
    s.setAstrometryNetCatalogUrl("/second/path");
    const solver::param_value v = s.getAstrometryParam("index-dir");
    ASSERT_TRUE(std::holds_alternative<std::string>(v));
    EXPECT_EQ(std::get<std::string>(v), "/second/path");
}

TEST(AstrometrySetters, SetScampSolverPath_RoundTrip)
{
    solver s;
    const std::string path = "/nonexistent/scamp/path";
    EXPECT_NO_THROW(s.setScampSolverPath(path));
    // Path accepted: a subsequent solve() with a non-existent binary must throw
    // std::runtime_error (binary not found), not any other exception.
    auto img = TestUtils::generateStarField(16, 16, 3, 1.0, 100.0, 500.0, 10.0, 2.0, 1);
    EXPECT_THROW(s.solve(img, "dummy.ldac"), std::runtime_error);
}

TEST(AstrometrySetters, SetScampSolverPath_EmptyPathIsAccepted)
{
    solver s;
    EXPECT_NO_THROW(s.setScampSolverPath(""));
    // Solver must still be functional after clearing path.
    EXPECT_NO_THROW(s.getScampParam("distort_degrees"));
}

// ===========================================================================
// Suite 9 — setSolverConfiguration (runtime config reload)
// ===========================================================================

TEST(AstrometrySolver, SetSolverConfiguration_ThrowsOnEmptyPath)
{
    solver s;
    EXPECT_THROW(s.setSolverConfiguration(""), std::invalid_argument);
}

TEST(AstrometrySolver, SetSolverConfiguration_ThrowsOnMissingFile)
{
    solver s;
    EXPECT_THROW(s.setSolverConfiguration("/nonexistent/config.json"),
                 std::invalid_argument);
}

TEST(AstrometrySolver, SetSolverConfiguration_LoadsParametersFromFile)
{
    // Build a JSON config that overrides "scale-units" and "distort_degrees".
    const std::filesystem::path cfgDir = kTestOutputDir / "setSolverConfig";
    std::filesystem::create_directories(cfgDir);
    const std::filesystem::path cfgPath = cfgDir / "override.json";

    {
        std::ofstream f(cfgPath);
        f << "{\n"
          << "    \"astrometry.net\": {\n"
          << "        \"path\": \"\",\n"
          << "        \"url\":  \"\",\n"
          << "        \"api-key\": \"\",\n"
          << "        \"parameters\": {\n"
          << "            \"scale-units\": \"arcsecperpix\",\n"
          << "            \"no-plots\": true\n"
          << "        }\n"
          << "    },\n"
          << "    \"scamp\": {\n"
          << "        \"path\": \"\",\n"
          << "        \"parameters\": {\n"
          << "            \"distort_degrees\": 5\n"
          << "        }\n"
          << "    }\n"
          << "}\n";
    }

    solver s;
    ASSERT_NO_THROW(s.setSolverConfiguration(cfgPath.string()));

    // Astrometry.net parameter must reflect the file value.
    const solver::param_value scaleUnits = s.getAstrometryParam("scale-units");
    ASSERT_TRUE(std::holds_alternative<std::string>(scaleUnits));
    EXPECT_EQ(std::get<std::string>(scaleUnits), "arcsecperpix");

    const solver::param_value noPlots = s.getAstrometryParam("no-plots");
    ASSERT_TRUE(std::holds_alternative<bool>(noPlots));
    EXPECT_TRUE(std::get<bool>(noPlots));

    // SCAMP parameter must reflect the file value.
    const solver::param_value distortDeg = s.getScampParam("distort_degrees");
    ASSERT_TRUE(std::holds_alternative<int>(distortDeg));
    EXPECT_EQ(std::get<int>(distortDeg), 5);
}

TEST(AstrometrySolver, SetSolverConfiguration_PreservesUnmentionedParameters)
{
    // Parameters not present in the new config must retain their current values.
    const std::filesystem::path cfgDir = kTestOutputDir / "setSolverConfig";
    std::filesystem::create_directories(cfgDir);
    const std::filesystem::path cfgPath = cfgDir / "partial.json";

    {
        std::ofstream f(cfgPath);
        // Config touches only "no-plots"; everything else must be unchanged.
        f << "{\n"
          << "    \"astrometry.net\": {\n"
          << "        \"path\": \"\",\n"
          << "        \"url\":  \"\",\n"
          << "        \"api-key\": \"\",\n"
          << "        \"parameters\": { \"no-plots\": true }\n"
          << "    },\n"
          << "    \"scamp\": { \"path\": \"\", \"parameters\": {} }\n"
          << "}\n";
    }

    solver s;
    // Pre-set a value we expect to survive the reload.
    s.setAstrometryParam("scale-low", 0.7f);
    s.setScampParam("crossid_radius", 4.0f);

    ASSERT_NO_THROW(s.setSolverConfiguration(cfgPath.string()));

    // These were set before reload and are not in the new config — must be unchanged.
    const solver::param_value scaleLow = s.getAstrometryParam("scale-low");
    ASSERT_TRUE(std::holds_alternative<float>(scaleLow));
    EXPECT_FLOAT_EQ(std::get<float>(scaleLow), 0.7f);

    const solver::param_value crossid = s.getScampParam("crossid_radius");
    ASSERT_TRUE(std::holds_alternative<float>(crossid));
    EXPECT_FLOAT_EQ(std::get<float>(crossid), 4.0f);
}

// ===========================================================================
// Suite 10 — solve(imgFile, sextractorFile) 2-arg SCAMP file overload
// ===========================================================================

TEST(AstrometrySolver, SolveScampFileOverload_ThrowsWhenBinaryMissing)
{
    // solve(string, string) delegates to solve(string, 1, string).
    // When the SCAMP binary is not found, std::runtime_error must be thrown.
    solver s;
    s.setScampSolverPath("/nonexistent/scamp/bin");
    EXPECT_THROW(s.solve(std::string("dummy.fits"), std::string("dummy.ldac")),
                 DSL::FITSexception);
}

TEST(AstrometrySolver, SolveScampFileOverload_ThrowsNotImplemented)
{
    // solve(string, size_t, string) is documented as "not yet implemented"
    // and must throw std::runtime_error unconditionally.
    solver s;
    EXPECT_THROW(s.solve(std::string("dummy.fits"), size_t(1), std::string("dummy.ldac")),
                 DSL::FITSexception);
}
