#ifndef CANARY_TESTS_DEVICE_HOST_BACKEND_HPP
#define CANARY_TESTS_DEVICE_HOST_BACKEND_HPP

// Execution of the checks in the host companion. Storage is host memory and a
// launch is an ordinary call of the entry. The CUDA and HIP backends
// (cuda_backend.hpp, hip_backend.hpp) provide the same functions, and
// description() for the runtime executable.

#include <cstddef>
#include <cstring>
#include <new>

namespace canary::test::backend
{

// Bytes of storage aligned for any type whose alignment does not exceed
// __STDCPP_DEFAULT_NEW_ALIGNMENT__.
inline void* allocate_bytes(std::size_t bytes) { return ::operator new(bytes); }

inline void release(void* data) noexcept { ::operator delete(data); }

inline void copy_to_device(void* destination, const void* source, std::size_t bytes)
{ std::memcpy(destination, source, bytes); }

inline void copy_to_host(void* destination, const void* source, std::size_t bytes)
{ std::memcpy(destination, source, bytes); }

// Calls Entry once with args. Returns after the call completes.
template <auto Entry, typename... Args> void launch(const Args&... args) { Entry(args...); }

}  // namespace canary::test::backend

#endif  // CANARY_TESTS_DEVICE_HOST_BACKEND_HPP
