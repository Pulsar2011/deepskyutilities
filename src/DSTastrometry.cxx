/**
 * @file DSTastrometry.cxx
 * @brief Implementation of DST::astrometry::solver.
 * @author GILLARD William — Centre de Physique des Particules de Marseille
 *
 * @copyright Licensed under PolyForm Noncommercial 1.0.0.
 * You may use, share and adapt this code for noncommercial purposes.
 * Commercial use requires a separate licence from the author.
 * Full licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/
 *
 * @details
 * This file contains the out-of-line implementations of all `solver` member
 * functions declared in DSTastrometry.h.  The code is split into four logical
 * sections:
 *
 * 1. **Initialisation** (`init`, `_checkForAstrometryNetSolver`,
 *    `_checkForScampSolver`) — parse JSON config and validate binary paths.
 *
 * 2. **Parameter helpers** (anonymous namespace) — `typeCompatible`, `checkKey`,
 *    `checkTypeAndSet` — shared validation logic used by all typed parameter setters.
 *
 * 3. **Parameter accessors** (`setAstrometryParam`, `getAstrometryParam`, …,
 *    `setScampParam`, `getScampParam`, …) — thin wrappers around the helpers.
 *
 * 4. **Solve overloads** — the actual plate-solving pipeline for each combination
 *    of input type (FITScube / file path) and engine (Astrometry.net / SCAMP /
 *    user star list).
 *
 * ### Thread-safety notes
 * Each `solve(file, hdu)` call makes a **thread-local copy** of
 * `astrometry_params` before building the command line.  This prevents data
 * races when multiple threads call `solve()` concurrently on the same `solver`
 * instance: the output path (`"new-fits"`) is patched only in the local copy,
 * leaving the shared `astrometry_params` untouched.
 *
 * SCAMP's `scamp_params` is read (not written) during `solve()`, so it is safe
 * to read from multiple threads simultaneously as long as no `setScampParam()`
 * call is in progress.
 */

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <chrono>
#include <string>

#include <DSTools/DSTastrometry.h>
#include <nlohmann/json.hpp>

using namespace DSL;

namespace DST
{
    namespace astrometry
    {
        void solver::init(const std::string& configJSON)
        {
            if(configJSON.empty() || configJSON.size() < 1)
            {
                throw std::invalid_argument("Configuration file is empty. ["+std::to_string(__LINE__)+"]");
            }

            if(!std::filesystem::exists(configJSON))
            {
                throw std::invalid_argument("Configuration file not found: " + configJSON + " ["+std::to_string(__LINE__)+"]");
            }

            const char* tmpdir = std::getenv("TMPDIR");
            if(tmpdir != nullptr && std::filesystem::exists(tmpdir))
            {
                astrometry_params.params["temp-dir"] = std::string(tmpdir);
            }
            else
            {
                std::filesystem::path tmpPath = std::filesystem::temp_directory_path();
                astrometry_params.params["temp-dir"] = tmpPath.string();
            }

            std::ifstream configFile(configJSON);

            nlohmann::json config;
            configFile >> config;

            if (config.contains("astrometry.net"))
            {
                const nlohmann::json& an = config["astrometry.net"];

                if (an.contains("path") || an.contains("PATH"))
                    astrometry_params.solver_path  = an["path"].get<std::string>();
                if (an.contains("url") || an.contains("URL"))
                    astrometry_params.solver_url   = an["url"].get<std::string>();
                if (an.contains("api-key") || an.contains("API_KEY"))
                    astrometry_params.api_key      = an["api-key"].get<std::string>();

                if (an.contains("parameters"))
                {
                    const nlohmann::json& an_params = an["parameters"];

                    for (nlohmann::json::const_iterator it = an_params.begin(); it != an_params.end(); ++it)
                    {
                        std::string     key=std::string();
                        key          += it.key();
                        std::transform(key.begin(), key.end(), key.begin(), ::tolower);  // case-insensitive keys

                        const nlohmann::json& jval = it.value();

                        if (astrometry_params.params.find(key) == astrometry_params.params.cend())
                            continue;  // unknown parameter — skip

                        param_value& current = astrometry_params.params[key];

                        if (std::holds_alternative<bool>(current))
                        {
                            current = jval.get<bool>();
                        }
                        else if (std::holds_alternative<std::string>(current))
                        {
                            current = jval.get<std::string>();
                        }
                        else if (std::holds_alternative<float>(current))
                        {
                            current = jval.get<float>();
                        }
                        else if (std::holds_alternative<int>(current))
                        {
                            current = jval.get<int>();
                        }
                        else  // std::monostate: type not yet set, infer from JSON
                        {
                            if (jval.is_boolean())
                                current = jval.get<bool>();
                            else if (jval.is_string())
                                current = jval.get<std::string>();
                            else if (jval.is_number_integer())
                                current = jval.get<int>();
                            else if (jval.is_number_float())
                                current = jval.get<float>();
                        }
                    }
                }
            }

            if (config.contains("scamp"))
            {
                const nlohmann::json& sc = config["scamp"];

                if (sc.contains("path"))
                    scamp_params.solver_path   = sc["path"].get<std::string>();

                if (sc.contains("parameters"))
                {
                    const nlohmann::json& sc_params = sc["parameters"];

                    for (nlohmann::json::const_iterator it = sc_params.begin(); it != sc_params.end(); ++it)
                    {
                        std::string     key=std::string();
                        key          += it.key();
                        std::transform(key.begin(), key.end(), key.begin(), ::tolower);  // case-insensitive keys
                        
                        const nlohmann::json& jval = it.value();

                        if (scamp_params.params.find(key) == scamp_params.params.cend())
                            continue;  // unknown parameter — skip

                        param_value& current = scamp_params.params[key];

                        if (std::holds_alternative<bool>(current))
                        {
                            current = jval.get<bool>();
                        }
                        else if (std::holds_alternative<std::string>(current))
                        {
                            current = jval.get<std::string>();
                        }
                        else if (std::holds_alternative<float>(current))
                        {
                            current = jval.get<float>();
                        }
                        else if (std::holds_alternative<int>(current))
                        {
                            current = jval.get<int>();
                        }
                        else  // std::monostate: type not yet set, infer from JSON
                        {
                            if (jval.is_boolean())
                                current = jval.get<bool>();
                            else if (jval.is_string())
                                current = jval.get<std::string>();
                            else if (jval.is_number_integer())
                                current = jval.get<int>();
                            else if (jval.is_number_float())
                                current = jval.get<float>();
                        }
                    }
                }
            }

        }

        void solver::_checkForAstrometryNetSolver()
        {
            if(astrometry_params.solver_path.empty())
            {
                throw std::runtime_error("Astrometry.net solver not found. Please provide a valid path or URL. ["+std::to_string(__LINE__)+"]");
            }

            bool solverFound = false;
            if(std::filesystem::exists(astrometry_params.solver_path+"/solve-field"))
            {
                solverFound |= true;
            }
            else
            {
                solverFound &= false;
            }

            if(!solverFound)
            {
                throw std::runtime_error("Astrometry.net solver not found. Install astrometry.net command line tools to use solver. ["+std::to_string(__LINE__)+"]");
            }
        }

        void solver::_checkForScampSolver()
        {
            if(scamp_params.solver_path.empty())
            {
                throw std::runtime_error("SCAMP solver not found. Please provide a valid path. ["+std::to_string(__LINE__)+"]");
            }

            bool solverFound = false;
            if(std::filesystem::exists(scamp_params.solver_path+"/scamp"))
            {
                solverFound |= true;
            }
            else
            {
                solverFound &= false;
            }

            if(!solverFound)
            {
                throw std::runtime_error("SCAMP solver not found. Install SCAMP command line tools to use solver. ["+std::to_string(__LINE__)+"]");
            }
        }

        // ---------------------------------------------------------------------------
        // Parameter access helpers
        // ---------------------------------------------------------------------------

        namespace
        {
            // Returns true if the new value's type is compatible with the slot already
            // stored: monostate slots accept any type; otherwise types must match.
            bool typeCompatible(const solver::param_value& existing,
                                const solver::param_value& newVal)
            {
                if (std::holds_alternative<std::monostate>(existing))
                    return true;
                return existing.index() == newVal.index();
            }

            void checkKey(const std::map<std::string, solver::param_value>& params,
                          const std::string& key,
                          const std::string& tool)
            {
                if (params.find(key) == params.cend())
                    throw std::invalid_argument("Unknown " + tool + " parameter: \"" + key + "\"");
            }

            void checkTypeAndSet(std::map<std::string, solver::param_value>& params,
                                 const std::string& key,
                                 const solver::param_value& newVal,
                                 const std::string& tool)
            {
                auto it = params.find(key);
                if (it == params.end())
                    throw std::invalid_argument("Unknown " + tool + " parameter: \"" + key + "\"");
                if (!typeCompatible(it->second, newVal))
                    throw std::invalid_argument("Type mismatch for " + tool + " parameter \"" + key + "\"");
                it->second = newVal;
            }
        }  // anonymous namespace

        // --- astrometry.net ---

        void solver::setAstrometryParam(const std::string& key, bool val)
        { checkTypeAndSet(astrometry_params.params, key, param_value{val}, "astrometry.net"); }

        void solver::setAstrometryParam(const std::string& key, std::string val)
        { checkTypeAndSet(astrometry_params.params, key, param_value{std::move(val)}, "astrometry.net"); }

        void solver::setAstrometryParam(const std::string& key, float val)
        { checkTypeAndSet(astrometry_params.params, key, param_value{val}, "astrometry.net"); }

        void solver::setAstrometryParam(const std::string& key, int val)
        { checkTypeAndSet(astrometry_params.params, key, param_value{val}, "astrometry.net"); }

        void solver::unsetAstrometryParam(const std::string& key)
        {
            checkKey(astrometry_params.params, key, "astrometry.net");
            astrometry_params.params[key] = std::monostate{};
        }

        solver::param_value solver::getAstrometryParam(const std::string& key) const
        {
            checkKey(astrometry_params.params, key, "astrometry.net");
            return astrometry_params.params.at(key);
        }

        bool solver::isAstrometryParamSet(const std::string& key) const
        {
            auto it = astrometry_params.params.find(key);
            return it != astrometry_params.params.cend() &&
                   !std::holds_alternative<std::monostate>(it->second);
        }

        // --- SCAMP ---

        void solver::setScampParam(const std::string& key, bool val)
        { checkTypeAndSet(scamp_params.params, key, param_value{val}, "SCAMP"); }

        void solver::setScampParam(const std::string& key, std::string val)
        { checkTypeAndSet(scamp_params.params, key, param_value{std::move(val)}, "SCAMP"); }

        void solver::setScampParam(const std::string& key, float val)
        { checkTypeAndSet(scamp_params.params, key, param_value{val}, "SCAMP"); }

        void solver::setScampParam(const std::string& key, int val)
        { checkTypeAndSet(scamp_params.params, key, param_value{val}, "SCAMP"); }

        void solver::unsetScampParam(const std::string& key)
        {
            checkKey(scamp_params.params, key, "SCAMP");
            scamp_params.params[key] = std::monostate{};
        }

        solver::param_value solver::getScampParam(const std::string& key) const
        {
            checkKey(scamp_params.params, key, "SCAMP");
            return scamp_params.params.at(key);
        }

        bool solver::isScampParamSet(const std::string& key) const
        {
            auto it = scamp_params.params.find(key);
            return it != scamp_params.params.cend() &&
                   !std::holds_alternative<std::monostate>(it->second);
        }

        // ---------------------------------------------------------------------------

        void solver::solve(const std::shared_ptr<FITScube> img)
        {
            if (!img)
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Input FITScube is null ["+std::to_string(__LINE__)+"]");

            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
            std::chrono::seconds unix_timestamp = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch());

            std::string timeStamp = std::to_string(unix_timestamp.count());

            // 1. Write FITScube to temporary FITS file
            std::string tempDir=std::string();
            tempDir += std::get<std::string>(astrometry_params.params["temp-dir"]);

            if(tempDir.empty())
                tempDir += std::filesystem::temp_directory_path().string();

            tempDir += ((tempDir.back() == '/') ? std::string("") : std::string("/"))
                + std::string("astrometry_temp_")
                + timeStamp
                + std::string(".fits");

            img->Write(tempDir);
            
            // delegate solving to the string-based solve method, which will read the temporary FITS file, solve it, and write the WCS back to the same file
            solve(tempDir,1);

            // Once solver is done, read the WCS from the temporary FITS file and update the FITScube's WCS header
            FITSmanager tempFits(tempDir);
            size_t N = tempFits.NumberOfHeader();  // force reading the file and populating internal structures
            std::shared_ptr<FITScube> tempImg = tempFits.GetImageAtIndex(N);  // get the last HDU, which should be the one we wrote and solved
            
            if(!tempImg)
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Failed to read temporary FITS file "+tempDir+" after solving ["+std::to_string(__LINE__)+"]");

            tempImg->reLoadWCS();  // ensure WCS is up-to-date with the HDU header before reading
            const FITSwcs tmpWCS = tempImg->getWCS();
            FITShdu wcsHDU = tmpWCS.asFITShdu();
            
            for(FITSDictionary::const_iterator it = wcsHDU.begin(); it != wcsHDU.end(); ++it)
            {
                img->HDU().ValueForKey(it->first, it->second.value(),it->second.type(),it->second.comment());
            }

            // Clean up temporary file
            std::filesystem::remove(tempDir);
        }

        void solver::solve(const std::string& imgFile, const size_t hdu)
        {
            _checkForAstrometryNetSolver();

            _solve_slot_guard slot(*this);  // blocks until a concurrency slot is free; released on scope exit

            // Compute output filename before launching the thread so outString is
            // in scope both inside the lambda and in the post-join cleanup below.
            bool autoGeneratedOutput = false;
            std::string outString;

            {
                const param_value& newFitsVal = astrometry_params.params["new-fits"];
                if(std::holds_alternative<std::string>(newFitsVal) &&
                   std::get<std::string>(newFitsVal) == "none")
                {
                    std::filesystem::path p(imgFile);
                    std::filesystem::path stem = p.filename();

                    // Strip compound extension (.fits.gz, .fit.gz) or any single extension
                    if(stem.extension() == ".gz")
                        stem = stem.stem();   // drop .gz → leaves .fits / .fit
                    stem = stem.stem();       // drop remaining extension

                    outString = (p.parent_path() / (stem.string() + "_wcs.fits")).string();
                    autoGeneratedOutput = true;
                }
                else if(std::holds_alternative<std::string>(newFitsVal))
                {
                    outString = std::get<std::string>(newFitsVal);
                }
            }

            std::thread solverThread([this, imgFile, outString]()
            {
                try
                {
                    // Build a per-call local copy of params with the resolved output path so
                    // we never mutate shared solver state (which would cause data races on
                    // concurrent solve() calls on the same solver instance).
                    astrometry_config_parameters local_params = astrometry_params;
                    local_params.params["new-fits"] = outString;
                    std::string cmdArgs = local_params.build_cmdline(imgFile);

                    // Execute the solver command
                    int result = std::system(cmdArgs.c_str());
                    if (result != 0)
                    {
                        throw std::runtime_error("Astrometry solver failed with exit code: " + std::to_string(result) + " ["+std::to_string(__LINE__)+"]");
                    }
                }
                catch (const std::exception& ex)
                {
                    std::cerr << "Error during astrometry solving: " << ex.what() << std::endl;
                }
            });

            solverThread.join();

            // Solver did not produce an output file — no WCS solution found, return silently
            if(!std::filesystem::exists(outString))
                return;

            // Read the WCS solution from the output file, then patch each keyword
            // directly into the existing HDU of the input file.
            // AppendKeyToHeader uses fits_update_key — it replaces the keyword if it
            // already exists and adds it if it doesn't. No image data is read or rewritten.
            FITSmanager outFits(outString);
            if(outFits.Status() != 0 || !outFits.isOpen())
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Failed to read output FITS file "+outString+" after solving ["+std::to_string(__LINE__)+"]");

            std::shared_ptr<FITScube> outImg = outFits.GetPrimary();
            if(!outImg)
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Failed to read primary HDU from output FITS file "+outString+" after solving ["+std::to_string(__LINE__)+"]");

            outImg->reLoadWCS();  // ensure WCS is up-to-date with the HDU header before reading
            const FITSwcs outWCS = outImg->getWCS();
            FITShdu wcsHDU = outWCS.asFITShdu();

            FITSmanager inFits(imgFile,false);
            if(inFits.Status() != 0 || !inFits.isOpen())
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Failed to open input FITS file "+imgFile+" for WCS update ["+std::to_string(__LINE__)+"]");

            for(FITSDictionary::const_iterator it = wcsHDU.begin(); it != wcsHDU.end(); ++it)
            {
                inFits.AppendKeyToHeader(hdu, it->first, it->second);
            }
            inFits.Write();
            inFits.Close();

            // Clean up the output FITS file if it was auto-generated by this method
            if(autoGeneratedOutput)
            {
                std::filesystem::remove(outString);
            }
        }

        void solver::solve(const std::shared_ptr<FITScube> img, const std::string& sextractorFile)
        {
            _checkForScampSolver();

            _solve_slot_guard slot(*this);  // blocks until a concurrency slot is free; released on scope exit

            std::thread solverThread([this, sextractorFile]()
            {
                try
                {

                    // Build command line arguments based on configuration parameters
                    temp_config_handle cmdArgs = scamp_params.build_cmdline(sextractorFile);

                    // Execute the solver command.
                    // On success the temp config file is deleted automatically by
                    // cmdArgs destructor (no release() call).
                    // On failure the file is kept and its path is printed so the
                    // user can inspect the configuration that caused the problem.
                    int result = std::system(cmdArgs.command.c_str());
                    if (result != 0)
                    {
                        std::cerr << "Warning: SCAMP failed — config file retained for investigation: "
                                  << cmdArgs.config_path.string() << std::endl;
                        cmdArgs.release();
                        throw std::runtime_error("SCAMP solver failed with exit code: " + std::to_string(result) + " ["+std::to_string(__LINE__)+"]");
                    }

                }
                catch (const std::exception& ex)
                {
                    std::cerr << "Error during SCAMP solving: " << ex.what() << std::endl;
                }
            });
            
            solverThread.join();

            // Derive where SCAMP wrote its output header file.
            // Priority:
            //   1. header_name — explicit output path set by the caller via
            //      setScampParam("header_name", "/abs/path/to/output.head").
            //      Use this to redirect output away from the catalog directory.
            //   2. Fallback: <catalog_dir>/<catalog_stem><header_suffix>
            //      (SCAMP default, suffix controlled by "header_suffix" param).
            const std::string& header_name =
                std::get<std::string>(scamp_params.params.at("header_name"));
            const std::string& header_suffix =
                std::get<std::string>(scamp_params.params.at("header_suffix"));
            const std::filesystem::path catPath(sextractorFile);
            const std::filesystem::path headPath =
                header_name.empty()
                ? catPath.parent_path() / (catPath.stem().string() + header_suffix)
                : std::filesystem::path(header_name);

            // No .head file → SCAMP found no solution; return silently (same convention
            // as the astrometry.net solve path which returns silently when solve-field
            // produces no output FITS).
            if (!std::filesystem::exists(headPath))
                return;

            // SCAMP's .head file is plain-text FITS header cards, not a proper FITS file.
            // Strategy: inject those cards into a minimal temporary FITS file so the
            // existing FITSmanager / FITSwcs machinery can extract and apply the WCS —
            // the same pattern used in the astrometry.net solve() overload above.
            std::chrono::time_point<std::chrono::system_clock> now2 = std::chrono::system_clock::now();
            const std::string ts2 = std::to_string(
                std::chrono::duration_cast<std::chrono::seconds>(now2.time_since_epoch()).count());
            const std::filesystem::path tmpFits =
                std::filesystem::temp_directory_path() / ("scamp_head_" + ts2 + ".fits");

            {
                fitsfile* fptr = nullptr;
                int status = 0;
                const std::string overwritePath = "!" + tmpFits.string();
                fits_create_file(&fptr, overwritePath.c_str(), &status);

                // A 1×1 image: FITScube::ctor requires NAXIS > 0 to open the HDU.
                long naxes[2] = {1, 1};
                fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);

                std::ifstream headFile(headPath);
                std::string line;
                while (std::getline(headFile, line))
                {
                    if (line.size() < 8) continue;
                    // Skip structural cards already written by fits_create_img
                    // (SIMPLE, BITPIX, NAXIS / NAXIS1 / NAXIS2, EXTEND, END).
                    // Keeping them would produce cfitsio status 212 (illegal NAXIS value).
                    const std::string k = line.substr(0, 8);
                    if (k.compare(0, 6, "SIMPLE") == 0 ||
                        k.compare(0, 6, "BITPIX") == 0 ||
                        k.compare(0, 5, "NAXIS")  == 0 ||
                        k.compare(0, 6, "EXTEND") == 0 ||
                        k.compare(0, 3, "END")    == 0)
                        continue;
                    line.resize(80, ' ');  // FITS cards are exactly 80 characters
                    fits_write_record(fptr, line.c_str(), &status);
                    status = 0;  // tolerate per-card errors
                }
                fits_close_file(fptr, &status);
            }

            // Scope guard: tmpFits is deleted on scope exit for any code path
            // (normal return, FITSexception, or any other exception type).
            struct RemoveOnExit
            {
                const std::filesystem::path& path;
                ~RemoveOnExit() { std::filesystem::remove(path); }
            } tmpFitsGuard{tmpFits};

            try
            {
                FITSmanager headMgr(tmpFits.string());
                std::shared_ptr<FITScube> headImg = headMgr.GetPrimary();
                if (headImg)
                {
                    headImg->reLoadWCS();
                    const FITSwcs headWCS = headImg->getWCS();
                    FITShdu wcsHDU = headWCS.asFITShdu();
                    for (FITSDictionary::const_iterator it = wcsHDU.begin(); it != wcsHDU.end(); ++it)
                        img->HDU().ValueForKey(it->first, it->second.value(),
                                              it->second.type(), it->second.comment());
                    // Refresh the cached FITSwcs object from the updated in-memory HDU.
                    // getWCS() returns the cached member (fwcs); ValueForKey() only updates
                    // the HDU dictionary, not the cache, so reLoadWCS() is required.
                    img->reLoadWCS();
                }
            }
            catch (const FITSexception& ex)
            {
                std::cerr << "Warning: failed to extract WCS from SCAMP .head file: "
                          << ex.what() << std::endl;
            }

        }

        void solver::solve(const std::string& imgFile, const size_t hdu, const std::string& sextractorFile)
        {
            FITSmanager inFits(imgFile,false);
            if(inFits.Status() != 0 || !inFits.isOpen())
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Failed to open input FITS file "+imgFile+" for SCAMP solving ["+std::to_string(__LINE__)+"]");

            std::shared_ptr<FITScube> inImg = inFits.GetImageAtIndex(hdu);
            if(!inImg)
                throw FITSexception(BAD_DIMEN,"DST::astrometry::solver","solve","Failed to read HDU "+std::to_string(hdu)+" from input FITS file "+imgFile+" for SCAMP solving ["+std::to_string(__LINE__)+"]");    
            
            solve(inImg, sextractorFile);
            inImg->reLoadWCS();  // refresh cached FITSwcs after WCS update in solve()

            FITSwcs wcs = inImg->getWCS();
            FITShdu wcsHDU = wcs.asFITShdu();

            for(FITSDictionary::const_iterator it = wcsHDU.begin(); it != wcsHDU.end(); ++it)
            {
                inFits.AppendKeyToHeader(hdu, it->first, it->second);
            }

            inFits.Write();
            inFits.Close();
            inImg.reset();
        }
    }
}