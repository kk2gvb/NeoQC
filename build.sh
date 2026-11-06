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
CUDA_MODE=""
CUDA_MODE_SET_BY_CMAKE=0

for arg in "$@"; do
    case "$arg" in
        --clean)
            CLEAN=1
            ;;
        --cuda)
            CUDA_MODE="ON"
            ;;
        --no-cuda)
            CUDA_MODE="OFF"
            ;;
        -DNEOQC_ENABLE_CUDA=*)
            CUDA_MODE_SET_BY_CMAKE=1
            CMAKE_ARGS+=("$arg")
            ;;
        --help|-h)
            echo "Usage: $0 [--clean] [--cuda|--no-cuda] [cmake-args...]"
            echo ""
            echo "Options:"
            echo "  --clean          Remove build/ directory before configuring"
            echo "  --cuda           Build with the optional CUDA QC backend"
            echo "  --no-cuda        Build CPU-only (default)"
            echo "  cmake-args       Any extra arguments are passed to CMake"
            echo "                   (e.g. -DCMAKE_BUILD_TYPE=Debug or"
            echo "                   -DCMAKE_CUDA_ARCHITECTURES=86)"
            echo ""
            echo "Examples:"
            echo "  $0"
            echo "  $0 --clean"
            echo "  $0 -DCMAKE_BUILD_TYPE=Debug"
            echo "  $0 --cuda -DCMAKE_CUDA_ARCHITECTURES=86"
            exit 0
            ;;
        *)
            CMAKE_ARGS+=("$arg")
            ;;
    esac
done

if [[ -n "$CUDA_MODE" && $CUDA_MODE_SET_BY_CMAKE -eq 1 ]]; then
    echo "Error: use either --cuda/--no-cuda or -DNEOQC_ENABLE_CUDA=..., not both." >&2
    exit 2
fi

# Always set the default explicitly.  This prevents a previous --cuda build
# in the same build/ directory from silently persisting into a CPU-only build.
if [[ -n "$CUDA_MODE" ]]; then
    CMAKE_ARGS+=("-DNEOQC_ENABLE_CUDA=$CUDA_MODE")
elif [[ $CUDA_MODE_SET_BY_CMAKE -eq 0 ]]; then
    CMAKE_ARGS+=("-DNEOQC_ENABLE_CUDA=OFF")
fi

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
echo ""
echo "  CUDA-enabled build (GPU per-base QC; remaining QC and trimming use OpenMP):"
echo "    ./build.sh --cuda -DCMAKE_CUDA_ARCHITECTURES=86"
