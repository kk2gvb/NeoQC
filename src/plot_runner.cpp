#include "../include/plot_runner.h"
#include <iostream>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

bool PlotRunner::run(const std::string& tsvPath, const std::string& outDir) {
    // Создаём каталог для графиков, если его нет
    std::error_code ec;
    fs::create_directories(outDir, ec);
    if (ec) {
        std::cerr << "[warning] Cannot create plot directory '" 
                  << outDir << "': " << ec.message() << "\n";
        return false;
    }

    // Путь к скрипту (относительно корня проекта)
    const std::string scriptPath = "scripts/plot_results.py";

    if (!fs::exists(scriptPath)) {
        std::cerr << "[warning] Plot script not found: " << scriptPath << "\n";
        return false;
    }

    if (!fs::exists(tsvPath)) {
        std::cerr << "[warning] TSV file not found: " << tsvPath << "\n";
        return false;
    }

    // Формируем команду с экранированием путей
    std::string cmd = "python3 \"" + scriptPath + "\" \"" + tsvPath + "\" \"" + outDir + "\"";

    int rc = std::system(cmd.c_str());
    if (rc != 0) {
        std::cerr << "[warning] plot_results.py failed (rc=" << rc 
                  << "). TSV/JSON results are still valid.\n";
        return false;
    }

    return true;
}