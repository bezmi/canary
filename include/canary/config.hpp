#ifndef CANARY_CONFIG_HPP
#define CANARY_CONFIG_HPP

// Build settings a consumer may define before including any canary header.

// Annotation for every function that per-cell code calls, so that it can run
// on a GPU. Define it as the consumer's host-and-device annotation, for
// example AMREX_GPU_HOST_DEVICE or KOKKOS_FUNCTION. Empty by default.
#ifndef CANARY_HOST_DEVICE
    #define CANARY_HOST_DEVICE
#endif

// Declares a function inline and requires the compiler to inline every call
// to it, independent of the compiler's size limits. It is for per-cell
// functions that are called in hot loops. Place it after CANARY_HOST_DEVICE
// and in place of the inline keyword. A consumer may define it, for example
// as AMREX_FORCE_INLINE. With a compiler that has no such attribute it is the
// inline keyword only.
#ifndef CANARY_FORCE_INLINE
    #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        #define CANARY_FORCE_INLINE __forceinline__
    #elif defined(__GNUC__)
        #define CANARY_FORCE_INLINE inline __attribute__((always_inline))
    #elif defined(_MSC_VER)
        #define CANARY_FORCE_INLINE __forceinline
    #else
        #define CANARY_FORCE_INLINE inline
    #endif
#endif

// Check of a precondition in per-cell code. Define it when assert is not
// usable in device code. Defaults to assert, so it is disabled by NDEBUG.
#ifndef CANARY_ASSERT
    #include <cassert>
    #define CANARY_ASSERT(cond) assert(cond)
#endif

// Species capacity of generic mixture models (the default Capacity template
// argument). Named models use their exact species count.
#ifndef CANARY_MAX_SPECIES
    #define CANARY_MAX_SPECIES 16
#endif

#endif  // CANARY_CONFIG_HPP
