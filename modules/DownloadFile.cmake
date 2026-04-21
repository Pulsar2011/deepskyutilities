# Helper script invoked by download targets via cmake -P (script mode).
# Called as: cmake [-DSKIP_IF_EXISTS=1] [-DOPTIONAL=1] -DURL=<url> -DDEST=<path> -P DownloadFile.cmake
#
#   SKIP_IF_EXISTS  -- silently return when DEST already exists (e.g. downloaded at configure time)
#   OPTIONAL        -- treat a failed download as a warning rather than a fatal error
#
# NOTE: This runs in CMake script mode (-P).  Project commands such as
#       find_program(), find_package(), etc. are NOT available.  We use
#       execute_process(COMMAND which ...) to locate external tools instead.
#
# IMPORTANT: All diagnostic messages use message() WITHOUT a mode keyword
# (not STATUS) so they go to stderr -- the same stream as CMake errors.
# message(STATUS ...) goes to stdout in CMake 3.15+ which CI runners
# may silently discard for failed commands.
cmake_minimum_required(VERSION 3.5)

message("DownloadFile.cmake [START] ========================================")
message("  CMake version : ${CMAKE_VERSION}")
message("  URL           : ${URL}")
message("  DEST          : ${DEST}")
message("  SKIP_IF_EXISTS: ${SKIP_IF_EXISTS}")
message("  OPTIONAL      : ${OPTIONAL}")
message("  System        : ${CMAKE_HOST_SYSTEM_NAME}")
message("================================================================")

if(NOT URL OR NOT DEST)
    message(FATAL_ERROR "DownloadFile.cmake requires -DURL=... -DDEST=...")
endif()

if(SKIP_IF_EXISTS AND EXISTS "${DEST}")
    message("  Already present, skipping download: ${DEST}")
    return()
endif()

# Check destination directory exists
get_filename_component(_dest_dir "${DEST}" DIRECTORY)
if(NOT EXISTS "${_dest_dir}")
    message("  WARNING: Destination directory does not exist: ${_dest_dir}")
else()
    message("  Destination directory OK: ${_dest_dir}")
endif()

message("  Downloading ${URL}")
message("           -> ${DEST}")

# Download strategy: prefer curl/wget (robust on minimal CI images that may
# lack CMake SSL support), then fall back to CMake's file(DOWNLOAD).
set(_download_ok FALSE)

# --- Attempt 1: curl (available on macOS and virtually all Linux CI images) ---
# find_program() is not available in script mode; use 'which' / 'where' instead.
message("")
message("  --- Attempt 1: curl ---")
if(WIN32)
    execute_process(COMMAND where curl OUTPUT_VARIABLE _CURL_PATH ERROR_QUIET RESULT_VARIABLE _which_curl_rc)
else()
    execute_process(COMMAND which curl OUTPUT_VARIABLE _CURL_PATH ERROR_QUIET RESULT_VARIABLE _which_curl_rc)
endif()

if(_which_curl_rc EQUAL 0)
    # Trim trailing newline from 'which' output
    string(STRIP "${_CURL_PATH}" _CURL_PATH)
    message("  Found curl: ${_CURL_PATH}")

    # Print curl version for diagnostics
    execute_process(
        COMMAND "${_CURL_PATH}" --version
        OUTPUT_VARIABLE _curl_version
        ERROR_QUIET
    )
    string(REGEX REPLACE "\n.*" "" _curl_version_line1 "${_curl_version}")
    message("  curl version: ${_curl_version_line1}")

    message("  Running: curl -fSL -v --retry 3 --retry-delay 5 -o <DEST> <URL>")
    execute_process(
        COMMAND "${_CURL_PATH}" -fSL -v --retry 3 --retry-delay 5 -o "${DEST}" "${URL}"
        RESULT_VARIABLE _curl_rc
        ERROR_VARIABLE  _curl_stderr
        OUTPUT_VARIABLE _curl_stdout
    )
    message("  curl exit code: ${_curl_rc}")
    if(_curl_stderr)
        message("  curl stderr:\n${_curl_stderr}")
    endif()
    if(_curl_rc EQUAL 0 AND EXISTS "${DEST}")
        message("  curl SUCCESS - file downloaded: ${DEST}")
        set(_download_ok TRUE)
    else()
        if(EXISTS "${DEST}")
            message("  curl FAILED but partial file exists - removing")
        else()
            message("  curl FAILED - no output file created")
        endif()
        file(REMOVE "${DEST}")
        message("  curl FAILED (rc=${_curl_rc})")
    endif()
else()
    message("  curl not found on this system (which rc=${_which_curl_rc})")
endif()

# --- Attempt 2: wget ---
if(NOT _download_ok)
    message("")
    message("  --- Attempt 2: wget ---")
    if(WIN32)
        execute_process(COMMAND where wget OUTPUT_VARIABLE _WGET_PATH ERROR_QUIET RESULT_VARIABLE _which_wget_rc)
    else()
        execute_process(COMMAND which wget OUTPUT_VARIABLE _WGET_PATH ERROR_QUIET RESULT_VARIABLE _which_wget_rc)
    endif()

    if(_which_wget_rc EQUAL 0)
        string(STRIP "${_WGET_PATH}" _WGET_PATH)
        message("  Found wget: ${_WGET_PATH}")

        # Print wget version for diagnostics
        execute_process(
            COMMAND "${_WGET_PATH}" --version
            OUTPUT_VARIABLE _wget_version
            ERROR_QUIET
        )
        string(REGEX REPLACE "\n.*" "" _wget_version_line1 "${_wget_version}")
        message("  wget version: ${_wget_version_line1}")

        message("  Running: wget -v --tries=3 -O <DEST> <URL>")
        execute_process(
            COMMAND "${_WGET_PATH}" -v --tries=3 -O "${DEST}" "${URL}"
            RESULT_VARIABLE _wget_rc
            ERROR_VARIABLE  _wget_stderr
            OUTPUT_VARIABLE _wget_stdout
        )
        message("  wget exit code: ${_wget_rc}")
        if(_wget_stderr)
            message("  wget stderr:\n${_wget_stderr}")
        endif()
        if(_wget_rc EQUAL 0 AND EXISTS "${DEST}")
            message("  wget SUCCESS - file downloaded: ${DEST}")
            set(_download_ok TRUE)
        else()
            if(EXISTS "${DEST}")
                message("  wget FAILED but partial file exists - removing")
            else()
                message("  wget FAILED - no output file created")
            endif()
            file(REMOVE "${DEST}")
            message("  wget FAILED (rc=${_wget_rc})")
        endif()
    else()
        message("  wget not found on this system (which rc=${_which_wget_rc})")
    endif()
endif()

# --- Attempt 3: CMake file(DOWNLOAD) ---
# Requires CMake built with SSL/TLS support + CA certificates on the system.
if(NOT _download_ok)
    message("")
    message("  --- Attempt 3: CMake file(DOWNLOAD) ---")
    message("    URL       : ${URL}")
    message("    DEST      : ${DEST}")
    message("    TLS_VERIFY: ON")
    file(DOWNLOAD
        "${URL}"
        "${DEST}"
        STATUS _status
        SHOW_PROGRESS
        TLS_VERIFY ON
        LOG _download_log
    )
    list(GET _status 0 _code)
    list(GET _status 1 _msg)
    message("  file(DOWNLOAD) status code: ${_code}")
    message("  file(DOWNLOAD) status msg : ${_msg}")
    if(_download_log)
        message("  file(DOWNLOAD) log:\n${_download_log}")
    endif()
    if(_code EQUAL 0)
        message("  file(DOWNLOAD) SUCCESS: ${DEST}")
        set(_download_ok TRUE)
    else()
        if(EXISTS "${DEST}")
            message("  file(DOWNLOAD) FAILED - partial file exists, removing")
        endif()
        message("  file(DOWNLOAD) FAILED: ${_msg}")
        file(REMOVE "${DEST}")
    endif()
endif()

# --- Final verdict ---
message("")
message("DownloadFile.cmake [RESULT] ======================================")
if(_download_ok)
    message("  DOWNLOAD OK: ${DEST}")
else()
    message("  DOWNLOAD FAILED - all 3 methods exhausted")
    message("    URL : ${URL}")
    message("    DEST: ${DEST}")
    if(OPTIONAL)
        message(WARNING "Download failed (non-fatal): all methods exhausted\n  URL: ${URL}")
    else()
        message(FATAL_ERROR "Download failed: all methods exhausted\n  URL: ${URL}")
    endif()
endif()
message("================================================================")
