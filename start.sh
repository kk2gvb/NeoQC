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

source venv/bin/activate

cd build
./qc $text