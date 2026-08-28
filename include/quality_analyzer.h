#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <array>
#include <unordered_map>
#include <memory>
#include "gc_model.h"
#include "fastq_reader.h"

constexpr std::size_t DUPLICATION_PREFIX_LENGTH = 50;
constexpr double OVERREPRESENTED_SEQUENCE_THRESHOLD = 0.1;

struct DuplicationKey {
    // 50 symbols from A/C/G/T/N encoded with three bits each (150 bits).
    std::array<uint64_t, 3> words{};
    bool operator==(const DuplicationKey&) const = default;
};

struct DuplicationEntry {
    DuplicationKey key;
    uint64_t count;
};

// Перенесно из fastq_reader.h, чтобы выполнять проверку сразу в QualityAnalyzer::processRecord, а не в FastqReader::readNext и не создавать лишние циклы
struct BaseValidationError {
    bool found = false;
    std::size_t recordNumber = 0;
    std::size_t position = 0;
    char base = '\0';
};

struct DuplicationKeyHash {
    std::size_t operator()(const DuplicationKey& key) const noexcept;
};

struct DuplicationLevelRow {
    std::string label;
    double totalSequencesPercent = 0.0;
    double deduplicatedSequencesPercent = 0.0;
};

struct OverrepresentedSequence {
    std::string sequence;
    uint64_t count = 0;
    double percent = 0.0;
};

struct DuplicationStats {
    std::vector<DuplicationLevelRow> levels;
    std::vector<OverrepresentedSequence> overrepresentedSequences;
    uint64_t totalReads = 0;
    uint64_t uniqueSequences = 0;
    double deduplicatedRemainingPercent = 100.0;
};

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

    // Реальное распределение NeoQC:
    // один read -> один GC bin
    std::vector<uint64_t> gcDistribution;

    // FastQC-compatible observed distribution:
    // один read может быть распределён между несколькими bins
    std::vector<double> gcDistributionFastQC;

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
    std::vector<uint64_t> perSequenceQualityDistributionTruncate;


    std::vector<uint64_t> baseCountA;
    std::vector<uint64_t> baseCountC;
    std::vector<uint64_t> baseCountG;
    std::vector<uint64_t> baseCountT;
    std::vector<uint64_t> baseCountN;

    std::vector<uint64_t> readsPerPosition;
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
    BaseValidationError processRecord(const FastqRecord& record);

    // Поиск адаптеров в записи
    void analyzeAdapters(const FastqRecord& record);

    // Получить итоговую статистику
    QualityStats getStats() const;

    // Получить точную статистику по всем уникальным 50-nt префиксам.
    DuplicationStats getDuplicationStats() const;

    DuplicationStats getDuplicationStats(const std::vector<DuplicationEntry>& entries) const;

    uint64_t getTotalReads() const { return totalReads; }

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

    // -----------------------------------------------------------------------
    // Объединение статистики с другим анализатором (для параллельной обработки)
    // -----------------------------------------------------------------------
    void merge(const QualityAnalyzer& other);

    std::vector<DuplicationEntry> getDuplicationEntries() const;
    

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

    // Все уникальные 50-нуклеотидные префиксы отслеживаются до конца файла.
    // Компактное 3-битное представление уменьшает стоимость точной карты.
    std::unordered_map<DuplicationKey, uint64_t, DuplicationKeyHash> sequenceCounts;

    std::vector<uint64_t> baseCountA;
    std::vector<uint64_t> baseCountC;
    std::vector<uint64_t> baseCountG;
    std::vector<uint64_t> baseCountT;
    std::vector<uint64_t> baseCountN;

    std::vector<uint64_t> readsPerPosition;

    std::vector<uint64_t> gcDistribution = std::vector<uint64_t>(101, 0);
    std::vector<double> gcDistributionFastQC = std::vector<double>(101, 0.0);
    std::vector<uint64_t> lengthDistribution;

    std::unordered_map<size_t, std::unique_ptr<GCModel>> gcModels;

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
    std::vector<uint64_t> perSequenceQualityDistributionTruncate;

};
