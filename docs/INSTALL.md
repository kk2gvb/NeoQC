# Установка NeoQC

## 1. Требования

Для сборки NeoQC необходимы:

- CMake `3.22+`;
- компилятор с поддержкой C++20;
- zlib;
- OpenMP;
- Python 3;
- Matplotlib `3.7+` (для графиков и HTML-отчёта).

Основное вычислительное ядро собирается на C++20. Python используется для построения графиков, формирования отчётов и QC-оценки.

Шрифты IBM Plex для графиков и HTML-отчёта входят в репозиторий (`assets/fonts/ibm-plex/`, лицензия SIL OFL 1.1) — устанавливать их в систему не нужно. Сам отчёт собирается без Matplotlib (`scripts/generate_qc_report.py`), если графики уже построены.

## 2. Проверка окружения

```bash
cmake --version
g++ --version
python3 --version
```

Для сборки также требуется поддержка OpenMP со стороны компилятора.

## 3. Python-окружение

При первой конфигурации CMake автоматически создаёт `venv/` в корне
репозитория и устанавливает зависимости из `requirements.txt`. Активировать
виртуальное окружение не нужно: NeoQC и CTest используют его автоматически.

При необходимости проверить установленную версию Matplotlib напрямую:

```bash
./venv/bin/python3 -c "import matplotlib; print(matplotlib.__version__)"
```

## 4. Сборка

Из корня репозитория:

```bash
./build.sh
```

Результат:

```text
build/neoqc
```

Ручная сборка через CMake:

```bash
cmake -S . -B build
cmake --build build -j$(nproc)
```

## 5. Проверка установки

```bash
ctest --test-dir build --output-on-failure
```

## 6. Проверка CLI

```bash
./build/neoqc --help
```

## 7. Генерация тестовых данных

Если сборочная конфигурация предоставляет соответствующую цель:

```bash
cmake --build build --target generate_test_data
```

Проверить созданные файлы:

```bash
find tests/data -maxdepth 1 -type f | sort
```

## 8. Типовая установка

```bash
git clone https://github.com/kk2gvb/NeoQC.git
cd NeoQC

./build.sh

ctest --test-dir build --output-on-failure

./build/neoqc --help
```

Файл `requirements.txt` обязателен: CMake использует его для автоматической
установки Python-зависимостей.

## 9. Обновление сборки

После изменения исходного кода:

```bash
cmake --build build -j$(nproc)
```

После изменения CMake-конфигурации:

```bash
cmake -S . -B build
cmake --build build -j$(nproc)
```

Затем:

```bash
ctest --test-dir build --output-on-failure
```

## 10. Диагностика

### Не найден CMake

Установите CMake средствами пакетного менеджера ОС.

### Не найден zlib

Установите пакет разработки zlib и повторите конфигурацию CMake.

### Не найден Python, venv или Matplotlib

Проверьте:

```bash
python3 --version
cmake -S . -B build
```

На этапе конфигурации CMake покажет причину, если не удалось создать `venv/`
или установить зависимости. После успешной конфигурации проверка выполняется
без активации окружения: `./venv/bin/python3 -c "import matplotlib"`.

## 11. Релизная проверка

Перед публикацией релиза:

```bash
cmake -S . -B build
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
./build/neoqc --help
```

После этого необходимо проверить версию, документацию, конфигурационные файлы и работоспособность single-end и paired-end запуска.
