#pragma once
#include <string>

class PlotRunner {
public:
    // Запускает Python-скрипт для построения графиков из TSV.
    // Возвращает true, если графики успешно построены.
    // При ошибке выводит предупреждение, но не прерывает работу программы.
    static bool run(const std::string& tsvPath, const std::string& outDir);
};