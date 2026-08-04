#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include "fastq_reader.h"

// ---------------------------------------------------------------------------
// Структура результатов анализа
// ---------------------------------------------------------------------------
struct QualityStats {
    uint64_t totalReads = 0;
    uint64_t totalBases = 0;
    uint64_t minLength  = 0;
    uint64_t maxLength  = 0;
    double   avgLength  = 0.0;

    uint64_t countA = 0;
    uint64_t countC = 0;
    uint64_t countG = 0;
    uint64_t countT = 0;
    uint64_t countN = 0;

    double avgGC      = 0.0;
    double percentN   = 0.0;
    double percentQ20 = 0.0;
    double percentQ30 = 0.0;
    double percentWithAdapter = 0.0;

    // Качество по позициям (среднее Phred-значение на каждой позиции)
    std::vector<double> meanQualityPerPosition;

    // Распределение среднего качества прочтений (гистограмма)
    std::vector<uint64_t> qualityDistribution; // индекс = Phred score

    // Распределение среднего Phred-качества по прочтениям.
    // Индекс = среднее качество прочтения, округлённое до ближайшего целого.
    std::vector<uint64_t> perSequenceQualityDistribution;
};

enum class ReadDirection {
    R1,
    R2
};

// ---------------------------------------------------------------------------
// Анализатор качества
// ---------------------------------------------------------------------------
class QualityAnalyzer {
public:
    explicit QualityAnalyzer(ReadDirection direction = ReadDirection::R1);

    // Обработка одной FASTQ-записи
    void processRecord(const FastqRecord& record);

    // Поиск адаптеров в записи
    void analyzeAdapters(const FastqRecord& record);

    // Получить итоговую статистику
    QualityStats getStats() const;

    // -----------------------------------------------------------------------
    // Публичные поля (нужны для вывода в TSV)
    // -----------------------------------------------------------------------
    // Список адаптеров (название, последовательность)
    struct Adapter {
        std::string name;
        std::string sequence;
        // Начальный k-mer, используемый для обнаружения adapter в риде.
        std::string detectionSequence;
    };

    std::vector<Adapter> adapters;

    // Кумулятивные счётчики: после обнаружения adapter на позиции значение
    // увеличивается до конца рида, как в FastQC Adapter Content.
    std::vector<std::vector<uint64_t>> adapterPosCounts;

private:
    // Базовые счётчики
    uint64_t totalGC   = 0;
    uint64_t totalBases = 0;
    uint64_t totalReads = 0;

    uint64_t countA = 0;
    uint64_t countC = 0;
    uint64_t countG = 0;
    uint64_t countT = 0;
    uint64_t countN = 0;

    uint64_t minLength = UINT64_MAX;
    uint64_t maxLength = 0;
    uint64_t totalLength = 0; // для avgLength

    uint64_t q20Count = 0;
    uint64_t q30Count = 0;
    uint64_t readsWithAdapter = 0;

    // Качество по позициям
    std::vector<uint64_t> qualitySum;
    std::vector<uint64_t> qualityCount;

    // Распределение среднего качества прочтений
    std::vector<uint64_t> qualityDistribution;
    std::vector<uint64_t> perSequenceQualityDistribution;
};
