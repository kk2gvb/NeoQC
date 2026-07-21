#!/bin/bash
set -e

if [ -n "$1" ]; then
    # Записываем первый аргумент в понятную переменную
    text="$1"
    echo "Передан путь к файлам: $text"
else
    echo "Вы не передали путь к файлам при запуске скрипта!"
    exit 1
fi

echo "Cleaning build..."
rm -rf build/* results/* plots/*
mkdir -p build results plots

cd build

echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

echo "Building..."
make -j$(nproc)

echo "Build complete."

echo "Setting up Python environment..."
cd ..
python -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
echo "Python environment ready."

cd build
./qc $text