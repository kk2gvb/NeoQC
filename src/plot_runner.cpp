#include "../include/plot_runner.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

#ifndef NEOQC_PLOT_SCRIPT
#define NEOQC_PLOT_SCRIPT "scripts/plot_results.py"
#endif

namespace {

std::string shellQuote(const std::string& value) {
    std::string quoted = "'";
    for (const char character : value) {
        if (character == '\'') {
            quoted += "'\"'\"'";
        } else {
            quoted += character;
        }
    }
    quoted += "'";
    return quoted;
}

}  // namespace

PlotRunResult PlotRunner::runAll(const std::string& resultDir,
                                 const std::string& plotDir,
                                 const PlotOptions& options) {
    PlotRunResult result;
    result.manifestPath = (fs::path(plotDir) / "plots_manifest.json").string();
    result.reportPath = (fs::path(resultDir) / "neoqc_qc_report.html").string();

    std::error_code error;
    fs::create_directories(plotDir, error);
    if (error) {
        std::cerr << "[warning] Cannot create plot directory '" << plotDir
                  << "': " << error.message() << "\n";
        return result;
    }

    const fs::path scriptPath = NEOQC_PLOT_SCRIPT;
    if (!fs::is_regular_file(scriptPath)) {
        std::cerr << "[warning] Plot script not found: " << scriptPath << "\n";
        return result;
    }
    if (!fs::is_directory(resultDir)) {
        std::cerr << "[warning] NeoQC result directory not found: " << resultDir << "\n";
        return result;
    }
    if (!options.generateSvg && !options.generatePng) {
        std::cerr << "[warning] Plot generation requested without an output format.\n";
        return result;
    }

    std::string command = "python3 " + shellQuote(scriptPath.string()) + " "
                        + shellQuote(resultDir) + " " + shellQuote(plotDir)
                        + " --formats";
    if (options.generateSvg) command += " svg";
    if (options.generatePng) command += " png";
    if (!options.includeAdapters) command += " --skip-adapters";
    if (options.strict) command += " --strict";

    const int returnCode = std::system(command.c_str());
    if (returnCode != 0) {
        std::cerr << "[warning] Plot generation failed (rc=" << returnCode
                  << "). TSV and summary results remain valid.\n";
        return result;
    }
    if (!fs::is_regular_file(result.manifestPath)) {
        std::cerr << "[warning] Plot generator did not create "
                  << result.manifestPath << "\n";
        return result;
    }
    if (!fs::is_regular_file(result.reportPath)) {
        std::cerr << "[warning] Plot generator did not create "
                  << result.reportPath << "\n";
        return result;
    }

    result.success = true;
    return result;
}
