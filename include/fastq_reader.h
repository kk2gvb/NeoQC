#pragma once

#include <string>
#include <cstdint>
#include <chrono>
#include <zlib.h>
#include <vector>

// ---------------------------------------------------------------------------
// FASTQ-запись
// ---------------------------------------------------------------------------
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string separator;
    std::string quality;
    std::uint64_t recordNumber = 0;
};


// ---------------------------------------------------------------------------
// Читатель FASTQ (поддерживает plain и .gz)
// ---------------------------------------------------------------------------
class FastqReader {
public:
    FastqReader(const std::string& filename);
    ~FastqReader();

    // Читает следующую запись. Возвращает false, если файл закончился.
    // Бросает исключение при ошибках формата или распаковки.
    bool readNext(FastqRecord& record);

    bool readBatch(std::vector<FastqRecord>& batch,
               std::size_t batchSize);

    // Количество прочитанных записей
    std::size_t getReadCount() const;

    // Имя файла (для сообщений об ошибках)
    const std::string& getFilename() const;


private:
    // Читает одну строку произвольной длины (без \n и \r)
    // Возвращает false, если достигнут конец файла
    bool readLine(std::string& line);

    // Удаляет \r и \n в конце строки
    static void trimNewlines(std::string& s);

    // Проверяет, что символ — допустимое основание (A/C/G/T/N, регистр не важен)
    static bool isValidBase(char c);

    gzFile fileHandle = nullptr;
    std::string filename;
    std::uint64_t readCount = 0;
};
