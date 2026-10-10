#pragma once

// Build settings a consumer may define before including any gaslib header.

// Annotation for every function that per-cell code calls, so that it can run
// on a GPU. Define it as the consumer's host-and-device annotation, for
// example AMREX_GPU_HOST_DEVICE or KOKKOS_FUNCTION. Empty by default.
#ifndef GASLIB_HOST_DEVICE
    #define GASLIB_HOST_DEVICE
#endif

// Check of a precondition in per-cell code. Define it when assert is not
// usable in device code. Defaults to assert, so it is disabled by NDEBUG.
#ifndef GASLIB_ASSERT
    #include <cassert>
    #define GASLIB_ASSERT(cond) assert(cond)
#endif
