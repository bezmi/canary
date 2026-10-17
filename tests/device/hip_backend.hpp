#ifndef CANARY_TESTS_DEVICE_HIP_BACKEND_HPP
#define CANARY_TESTS_DEVICE_HIP_BACKEND_HPP

// Execution of the checks on the current HIP device (see host_backend.hpp for
// the common functions). Every runtime call is checked; a failure throws
// std::runtime_error naming the call and the HIP error. release does not throw:
// an error it leaves is reported by the next checked call.

#include <cstddef>
#include <hip/hip_runtime.h>
#include <stdexcept>
#include <string>

namespace canary::test::backend
{

inline void check(hipError_t error, const char* operation)
{
    if (error != hipSuccess) {
        throw std::runtime_error(std::string(operation) + " failed: " + hipGetErrorName(error) +
                                 ": " + hipGetErrorString(error));
    }
}

inline void* allocate_bytes(std::size_t bytes)
{
    void* data = nullptr;
    check(hipMalloc(&data, bytes), "hipMalloc");
    return data;
}

inline void release(void* data) noexcept { (void)hipFree(data); }

inline void copy_to_device(void* destination, const void* source, std::size_t bytes)
{ check(hipMemcpy(destination, source, bytes, hipMemcpyHostToDevice), "hipMemcpy to device"); }

inline void copy_to_host(void* destination, const void* source, std::size_t bytes)
{ check(hipMemcpy(destination, source, bytes, hipMemcpyDeviceToHost), "hipMemcpy to host"); }

// Launches Entry with one block of one thread and waits for it to complete.
template <auto Entry, typename... Args> void launch(const Args&... args)
{
    Entry<<<1, 1>>>(args...);
    check(hipGetLastError(), "kernel launch");
    check(hipDeviceSynchronize(), "kernel execution");
}

// The current device, its architecture and the runtime and driver versions.
inline std::string description()
{
    int device = 0;
    check(hipGetDevice(&device), "hipGetDevice");
    hipDeviceProp_t properties {};
    check(hipGetDeviceProperties(&properties, device), "hipGetDeviceProperties");
    int runtime = 0;
    check(hipRuntimeGetVersion(&runtime), "hipRuntimeGetVersion");
    int driver = 0;
    check(hipDriverGetVersion(&driver), "hipDriverGetVersion");
    return "HIP device " + std::to_string(device) + " (" + properties.name + ", " +
           properties.gcnArchName + "), runtime " + std::to_string(runtime) + ", driver " +
           std::to_string(driver);
}

}  // namespace canary::test::backend

#endif  // CANARY_TESTS_DEVICE_HIP_BACKEND_HPP
