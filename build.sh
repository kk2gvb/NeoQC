#!/bin/bash
set -e

# Переходим в корень проекта, откуда бы ни был запущен скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# Обработка флагов самого build.sh
# ---------------------------------------------------------------------------
CLEAN=0
CMAKE_ARGS=()

for arg in "$@"; do
    case "$arg" in
        --clean)
            CLEAN=1
            ;;
        --help|-h)
            echo "Usage: $0 [--clean] [cmake-args...]"
            echo ""
            echo "Options:"
            echo "  --clean          Remove build/ directory before configuring"
            echo "  cmake-args       Any extra arguments are passed to CMake"
            echo "                   (e.g. -DCMAKE_BUILD_TYPE=Debug)"
            echo ""
            echo "Examples:"
            echo "  $0"
            echo "  $0 --clean"
            echo "  $0 -DCMAKE_BUILD_TYPE=Debug"
            exit 0
            ;;
        *)
            CMAKE_ARGS+=("$arg")
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Очистка (опционально)
# ---------------------------------------------------------------------------
if [[ $CLEAN -eq 1 ]]; then
    echo "=== Cleaning build directory ==="
    rm -rf build
fi

# ---------------------------------------------------------------------------
# Конфигурация
# ---------------------------------------------------------------------------
mkdir -p build
cd build

echo "=== Configuring with CMake ==="
cmake .. -DCMAKE_BUILD_TYPE=Release "${CMAKE_ARGS[@]}"

# ---------------------------------------------------------------------------
# Сборка
# ---------------------------------------------------------------------------
echo "=== Building ==="
# Определяем число ядер: nproc (Linux) или sysctl (macOS)
if command -v nproc >/dev/null 2>&1; then
    JOBS=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
    JOBS=$(sysctl -n hw.ncpu)
else
    JOBS=2
fi

cmake --build . --config Release -j"$JOBS"

echo ""
echo "=== Build complete ==="
echo "Binary: $SCRIPT_DIR/build/neoqc"
echo ""
echo "Run examples:"
echo "  Single-end:"
echo "    ./build/neoqc --r1 data/sample.fastq.gz --sample-id sample01 --out results/sample01"
echo ""
echo "  Paired-end:"
echo "    ./build/neoqc --r1 data/sample_R1.fastq.gz --r2 data/sample_R2.fastq.gz \\"
echo "                  --sample-id sample01 --out results/sample01"
echo ""
echo "  With plots (requires Python + scripts/plot_results.py):"
echo "    ./build/neoqc --r1 data/sample.fastq.gz --sample-id sample01 --out results/sample01 --plot"