#ifndef CANARY_TESTS_DEVICE_CUDA_BACKEND_HPP
#define CANARY_TESTS_DEVICE_CUDA_BACKEND_HPP

// Execution of the checks on the current CUDA device (see host_backend.hpp for
// the common functions). Every runtime call is checked; a failure throws
// std::runtime_error naming the call and the CUDA error. release does not throw:
// an error it leaves is reported by the next checked call.

#include <cstddef>
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace canary::test::backend
{

inline void check(cudaError_t error, const char* operation)
{
    if (error != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + " failed: " + cudaGetErrorName(error) +
                                 ": " + cudaGetErrorString(error));
    }
}

inline void* allocate_bytes(std::size_t bytes)
{
    void* data = nullptr;
    check(cudaMalloc(&data, bytes), "cudaMalloc");
    return data;
}

inline void release(void* data) noexcept { (void)cudaFree(data); }

inline void copy_to_device(void* destination, const void* source, std::size_t bytes)
{ check(cudaMemcpy(destination, source, bytes, cudaMemcpyHostToDevice), "cudaMemcpy to device"); }

inline void copy_to_host(void* destination, const void* source, std::size_t bytes)
{ check(cudaMemcpy(destination, source, bytes, cudaMemcpyDeviceToHost), "cudaMemcpy to host"); }

// Launches Entry with one block of one thread and waits for it to complete.
template <auto Entry, typename... Args> void launch(const Args&... args)
{
    Entry<<<1, 1>>>(args...);
    check(cudaGetLastError(), "kernel launch");
    check(cudaDeviceSynchronize(), "kernel execution");
}

// The current device, its compute capability and the runtime and driver versions.
inline std::string description()
{
    int device = 0;
    check(cudaGetDevice(&device), "cudaGetDevice");
    cudaDeviceProp properties {};
    check(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
    int runtime = 0;
    check(cudaRuntimeGetVersion(&runtime), "cudaRuntimeGetVersion");
    int driver = 0;
    check(cudaDriverGetVersion(&driver), "cudaDriverGetVersion");
    return "CUDA device " + std::to_string(device) + " (" + properties.name +
           ", compute capability " + std::to_string(properties.major) + "." +
           std::to_string(properties.minor) + "), runtime " + std::to_string(runtime) +
           ", driver " + std::to_string(driver);
}

}  // namespace canary::test::backend

#endif  // CANARY_TESTS_DEVICE_CUDA_BACKEND_HPP
