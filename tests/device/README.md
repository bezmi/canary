# Numerical device checks

These targets compile Canary's public numerical operations as CUDA or HIP
kernels, and the same sources as ordinary host functions. A successful build
shows that the operations compile for the selected device toolchain. It does not
run a device and does not check numerical results. Installed consumers never
build these targets and need no device toolchain.

## Selecting a backend

| Option | Effect |
| --- | --- |
| `CANARY_BUILD_CUDA_CHECKS=ON` | Enables the CUDA language and builds `canary_cuda_check` and `canary_device_host_check`. |
| `CANARY_BUILD_HIP_CHECKS=ON` | Enables the HIP language and builds `canary_hip_check` and `canary_device_host_check`. |

Both default to `OFF` and are independent of `CANARY_BUILD_TESTS`,
`CANARY_BUILD_DATABASE` and `CANARY_BUILD_BENCHMARKS`. A selected backend is
required: if its compiler or toolkit is missing, configuration fails instead of
skipping the checks.

Set the device architecture explicitly for a reproducible or headless build
(`CMAKE_CUDA_ARCHITECTURES` or `CMAKE_HIP_ARCHITECTURES`). Otherwise CMake and
the toolchain choose their default. Example for CUDA only:

```sh
cmake -S . -B build/cuda-checks -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCANARY_BUILD_TESTS=OFF -DCANARY_BUILD_DATABASE=OFF -DCANARY_BUILD_BENCHMARKS=OFF \
  -DCANARY_BUILD_CUDA_CHECKS=ON -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build build/cuda-checks --parallel 16
```

The `cuda-compile` and `hip-compile` workflow presets add the backend checks to
the full release host validation (`cmake --workflow --preset cuda-compile`).
The checks are split into one translation unit per required instance so that
they compile in parallel. Each CUDA unit needs about 1.2 GB of memory; with 16
jobs the thermodynamic checks build in about two minutes on a 24-core machine.
Lower the job count if memory is short.

## What is compiled

Each check reads one cell from caller storage, builds the views the public
operation accepts and writes the complete result to a one-element output view.
`thermodynamics_cases.hpp` holds the shared input preparation, the check entries
and the required-instance lists. `CMakeLists.txt` lists the required instances
once and generates one source per instance from `thermodynamics.in` for each
language (host `.cpp`, CUDA `.cu`, HIP `.hip`), covering:

| Parameter type P | State scalar T | Units |
| --- | --- | --- |
| float | float, Dual\<float\> | SI, scaled |
| double | double, Dual\<double\> | SI, scaled |

Thermally perfect energy closures are also compiled with P and T of different
precision. For every model and P/T/unit combination, every closure pair
((rho, e), (rho, rho_e), (rho, p), (rho, T), (p, T)) is compiled with contiguous
mass fractions. Each composition input (mass fractions, prepared mixture,
partial densities), each storage (contiguous, the middle cell of a layout_right
or layout_left grid), each mode input form (quantity mdspan, ModeIndex view) and
each mode input kind (temperatures, energy densities) is compiled at least once,
on (rho, T). Properties, species operations, electron pressure and NASA9 cp/h and
entropy are compiled for the models that provide them.

| Models |
| --- |
| Calorically perfect gas; perfect mixtures of 5 and 11 species at and below capacity |
| Thermally perfect: one-segment controlled species, nitrogen, five-species air |
| Two-temperature five-species air |
| Multi-temperature: four-temperature air5, three-temperature air11, controlled two-species model |
| NASA9 curves of one, two and four segments |
