#pragma once

#include <string>

// This header deliberately has no CUDA Runtime includes.  It is consumed by
// regular C++ code as well as the CUDA translation unit.
std::string runCudaBackendSelfTest();
