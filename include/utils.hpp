#pragma once
#include <string>
#include <iostream>

namespace Utils {
    // Убираем const&, чтобы работать с копией строки внутри функции
    inline std::string trim_path(std::string filepath, std::string folder) {

        // Исправлено: используем filepath вместо str
        size_t pos1 = filepath.find("../data/" + folder +"/"); 
        if (pos1 != std::string::npos) {
            filepath.erase(pos1, 9 + folder.length());
        }

        size_t pos2 = filepath.find(".fq.gz");
        if (pos2 != std::string::npos) {
            filepath.erase(pos2, 6); // Внимание: "data/" — это 5 символов, а не 6
        }

        size_t pos3 = filepath.find(".fastq.gz");
        if (pos3 != std::string::npos) {
            filepath.erase(pos3, 8);
        }

        return filepath; 
    }
}
