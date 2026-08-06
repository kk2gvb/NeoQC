#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <array>
#include <unordered_map>
#include "fastq_reader.h"

constexpr size_t DUPLICATION_PREFIX_LENGTH = 50;
constexpr double OVERREPRESENTED_SEQUENCE_THRESHOLD = 0.1;

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

    std::vector<uint64_t> gcDistribution;

    std::vector<uint64_t> lengthDistribution;

    double avgGC      = 0.0;
    double percentN   = 0.0;
    double percentQ20 = 0.0;
    double percentQ30 = 0.0;
    double percentWithAdapter = 0.0;

    // Качество по позициям (среднее Phred-значение на каждой позиции)
    std::vector<double> meanQualityPerPosition;
    std::vector<double> lowerQuartileQualityPerPosition;
    std::vector<double> medianQualityPerPosition;

    // Распределение среднего Phred-качества по прочтениям.
    // Индекс = среднее качество прочтения, округлённое до ближайшего целого.
    std::vector<uint64_t> perSequenceQualityDistribution;

    std::vector<uint64_t> baseCountA;
    std::vector<uint64_t> baseCountC;
    std::vector<uint64_t> baseCountG;
    std::vector<uint64_t> baseCountT;
    std::vector<uint64_t> baseCountN;

    std::vector<uint64_t> readsPerPosition;
    
    // Количество каждой уникальной последовательности
    std::unordered_map<std::string, uint64_t> sequenceCounts;

    struct OverrepresentedSequence
    {
        std::string sequence;
        uint64_t count;
        double percent;
    };

    std::vector<OverrepresentedSequence> overrepresentedSequences;
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

    // Частота каждой уникальной последовательности
    std::unordered_map<std::string, uint64_t> sequenceCounts;

    std::vector<uint64_t> baseCountA;
    std::vector<uint64_t> baseCountC;
    std::vector<uint64_t> baseCountG;
    std::vector<uint64_t> baseCountT;
    std::vector<uint64_t> baseCountN;

    std::vector<uint64_t> readsPerPosition;

    std::vector<uint64_t> gcDistribution = std::vector<uint64_t>(101, 0);
    std::vector<uint64_t> lengthDistribution;

    uint64_t minLength = UINT64_MAX;
    uint64_t maxLength = 0;
    uint64_t totalLength = 0; // для avgLength

    uint64_t q20Count = 0;
    uint64_t q30Count = 0;
    uint64_t readsWithAdapter = 0;

    // Качество по позициям
    std::vector<uint64_t> qualitySum;
    std::vector<uint64_t> qualityCount;
    std::vector<std::array<uint64_t, 94>> qualityHistogram;

    std::vector<uint64_t> perSequenceQualityDistribution;
};
