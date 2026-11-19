# Numerical device checks

These targets compile Canary's public numerical operations as CUDA or HIP
kernels, and the same sources as ordinary host functions. A successful build
shows that the operations compile for the selected device toolchain. The
optional runtime tests also run every kernel on the current device and check
the results (see [Running the checks](#running-the-checks)). Installed
consumers never build these targets and need no device toolchain.

## Selecting a backend

| Option | Effect |
| --- | --- |
| `CANARY_BUILD_CUDA_CHECKS=ON` | Enables the CUDA language and builds `canary_cuda_check` and `canary_device_host_check`. |
| `CANARY_BUILD_HIP_CHECKS=ON` | Enables the HIP language and builds `canary_hip_check` and `canary_device_host_check`. |
| `CANARY_BUILD_CUDA_RUNTIME_TESTS=ON` | Also builds `canary_cuda_runtime_tests` and registers its CTest tests. Requires `CANARY_BUILD_CUDA_CHECKS=ON`. |
| `CANARY_BUILD_HIP_RUNTIME_TESTS=ON` | Also builds `canary_hip_runtime_tests` and registers its CTest tests. Requires `CANARY_BUILD_HIP_CHECKS=ON`. |

All default to `OFF`. A runtime option without its check option is a
configuration error. The options are independent of `CANARY_BUILD_TESTS`,
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
cmake --build build/cuda-checks --parallel 20
```

The `cuda-compile` and `hip-compile` workflow presets add the backend checks to
the full release host validation (`cmake --workflow --preset cuda-compile`).
The checks are split into one translation unit per required instance so that
they compile in parallel. Each CUDA unit needs about 1.1 GB of memory; with 20
jobs the 64 CUDA units (about 5000 kernels) build in about 3.5 minutes on a
24-core machine. Lower the job count if memory is short.

## What is compiled

Each check reads one cell from caller storage, builds the views the public
operation accepts and writes the complete result to a one-element output view.
`cell_inputs.hpp` holds the shared input preparation. `thermodynamics_cases.hpp`,
`reactions_cases.hpp`, `exchange_cases.hpp` and `transport_cases.hpp` hold the
check entries of each subsystem and, in one group function per subsystem part,
the fixtures they are evaluated at: each group builds its models and inputs on
the host, runs every entry through `native_runtime.hpp` and checks the result.
An entry is a host-and-device function; a backend calls it from a kernel.
A group compiles exactly the entries it runs. `CMakeLists.txt`
lists the required instances once and generates one source per instance from
`check.in` for each language (host `.cpp`, CUDA `.cu`, HIP `.hip`), covering:

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

| Subsystem | Models and operations |
| --- | --- |
| Thermodynamics | Calorically perfect gas; perfect mixtures of 5 and 11 species at and below capacity; thermally perfect one-segment controlled species, nitrogen and five-species air; two-temperature air5; multi-temperature four-temperature air5, three-temperature air11 and a controlled two-species model; NASA9 curves of one, two and four segments |
| Reactions | `rates_of_progress` and `production_rates` for 13 reactions at one temperature, 5 reactions at two temperatures (scalar T_ve, owning modes, single temperature), 40 reactions at three temperatures (owning modes, each borrowed mode form and grid storage once, single temperature) and 12 interleaved reactions with four internal modes; multi-temperature rows also with P and T of different precision |
| Exchange | Every standalone term and the combined sources, for two-temperature air5 (full, empty, vibration-translation only, chemistry-vibration only), three-temperature air11 (full, empty, and each of the five processes alone; borrowed state view for the full exchange) and four-temperature air5 without free electrons |
| Transport | Constant, Sutherland and power-law viscosity, alone and with a constant Prandtl number on a perfect gas and mixture; Wilke and Gupta-Yos transport, full binary and mixture-averaged diffusion on one- and two-temperature air5 and on 1, 2 and 16 species (also with P and T of different precision) |

Exchange and transport cases read compositions from contiguous storage, except
where noted; the storage forms themselves are compiled by the other groups.

## Fixtures and derived inputs

`{thermodynamics,reactions,exchange,transport}_inputs.hpp` declare the fixture
models and states as binary64 values in SI units; `host_inputs.hpp` holds the
shared states, compositions, unit scales and seeding. A state input is rounded
to the value type of the state scalar, then converted to the unit system.

Inputs that a forward closure derives from a declared state (specific energy,
pressure, mode energies, and the density of a state declared by pressure) and
the rates of progress given to the standalone chemistry-coupled exchange terms
are generated offline into `derived_inputs.hpp` by `derive_inputs.cpp`, with
Canary's host closures in binary64 and SI units. They are inputs only, never
expected values. After changing a fixture or a model the fixtures use,
regenerate them in a tree configured with a check option and review the diff:

```sh
scripts/generate_device_inputs.sh build/cuda-compile
```

## Running the checks

With a runtime option selected, the check objects of the backend are linked
with a generated `main` into `canary_cuda_runtime_tests` or
`canary_hip_runtime_tests`. The executable takes one subsystem name
(`thermodynamics`, `reactions`, `exchange` or `transport`), runs every group
instance of that subsystem on the current device (one thread per evaluation)
and prints a report. CTest registers one test per subsystem,
`canary.device.<cuda|hip>.<subsystem>`, with the labels `canary.device` and the
backend name, and the resource lock `canary_device_runtime` so that device
tests do not run concurrently within one `ctest` invocation.

```sh
cmake --workflow --preset cuda-runtime
ctest --test-dir build/cuda-runtime -L canary.device --output-on-failure
```

Each fixture state is evaluated once without a derivative seed and, for Dual
state scalars, once per seeded input (see `host_inputs.hpp`). Every evaluation
runs on the device and, with the same inputs, as an ordinary function on the
host. Each number of the device result (a value, or the derivative of a Dual)
is compared with the same number of the host result. The two agree when they
are equal, or both NaN, or both finite with

```
|device - host| <= 1e-9 * max(|device|, |host|, scale)
```

where `scale` is 0 for a value and, for a derivative, the magnitude of the
value of the same field. The second form lets a derivative whose exact value
is zero differ by the rounding of the value. Two fixtures state wider limits
where the result is a difference of nearly equal numbers; the reason is given
at each (`nasa9_curve` in `thermodynamics_cases.hpp`, `two_temperature_terms`
in `exchange_cases.hpp`). The host results themselves are checked by the host
test suite, not here.

The report lists:

| Entry | Meaning | Test result |
| --- | --- | --- |
| error | A device call failed or a group could not run. | Fail |
| mismatch | An exact field (status, reason, diagnostic mode) of the device result differs from the fixture's declared value. | Fail |
| disagreement with the host | A number of the device result does not agree with the host result. | Fail |

A test passes only when at least one evaluation ran and nothing is listed. A
missing device or a device error is a failure, never a skip.
