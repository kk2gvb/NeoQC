#pragma once
#include <string>

struct PlotOptions {
    bool includeAdapters = true;
    bool generateSvg = true;
    bool generatePng = true;
    bool strict = true;
    std::string runManifestPath;
};

struct PlotRunResult {
    bool success = false;
    std::string manifestPath;
    std::string reportPath;
};

class PlotRunner {
public:
    // Checks the Python interpreter used by plot generation before any FASTQ
    // processing starts.  On failure, writes setup instructions to stderr.
    static bool checkDependencies();

    // Generates every available chart in one Python process.  A trimming
    // result root also produces a combined before/after HTML report. Plot
    // failures are reported but never invalidate the underlying QC results.
    static PlotRunResult runAll(const std::string& resultDir,
                                const std::string& plotDir,
                                const PlotOptions& options = {});
};
