# Установка NeoQC

## 1. Требования

Для сборки NeoQC необходимы:

- CMake `3.22+`;
- компилятор с поддержкой C++20;
- zlib;
- OpenMP;
- Python 3;
- Matplotlib.

Основное вычислительное ядро собирается на C++20. Python используется для построения графиков, формирования отчётов и QC-оценки.

## 2. Проверка окружения

```bash
cmake --version
g++ --version
python3 --version
```

Для сборки также требуется поддержка OpenMP со стороны компилятора.

## 3. Python-окружение

Рекомендуется использовать отдельное виртуальное окружение:

```bash
python3 -m venv venv
source venv/bin/activate
```

Если в репозитории присутствует `requirements.txt`:

```bash
pip install -r requirements.txt
```

Проверка Matplotlib:

```bash
python -c "import matplotlib; print(matplotlib.__version__)"
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

python3 -m venv venv
source venv/bin/activate

pip install -r requirements.txt

./build.sh

ctest --test-dir build --output-on-failure

./build/neoqc --help
```

Если `requirements.txt` отсутствует в конкретной ревизии, необходимо установить Python-зависимости, используемые скриптами проекта.

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

### Не найден Python или Matplotlib

Проверьте:

```bash
python3 --version
python3 -c "import matplotlib"
```

При необходимости:

```bash
source venv/bin/activate
```

## 11. Релизная проверка

Перед публикацией релиза:

```bash
cmake -S . -B build
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
./build/neoqc --help
```

После этого необходимо проверить версию, документацию, конфигурационные файлы и работоспособность single-end и paired-end запуска.
