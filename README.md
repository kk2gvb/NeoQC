# neo-mrna-vax
## Как запускать и собирать

Для сборки должны быть установлены:
- CMake (версии 3.22+)
- Cuda
- C++ (Стандарт 20+)

```shell
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release 
    cmake --build . -- -j$(nproc)
    ./qc ../data/SRR27872625_1.fastq.gz  
```