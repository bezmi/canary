# Canary

Canary is a header-only C++20 library of gas physics for compressible flow
solvers. Given the state of one cell, it returns thermodynamic properties,
transport coefficients, chemical production rates and energy-exchange source
terms. It has no mesh, no time integration and no I/O in its per-cell
functions; a flow solver calls it from its own loops, on a CPU or inside a GPU
kernel.

- [What it provides](#what-it-provides)
- [Requirements](#requirements)
- [Example](#example)
- [Using Canary in a solver](#using-canary-in-a-solver)
- [Precision and devices](#precision-and-devices)
- [Building and testing](#building-and-testing)
- [Known limitations](#known-limitations)
- [Licence](#licence)

## What it provides

| Part | Models |
| --- | --- |
| Thermodynamics | Calorically perfect gas and mixture; thermally perfect mixture with NASA 9-coefficient fits; two-temperature mixture (translation-rotation, vibration-electronic); multi-temperature mixture with separate vibrational, electronic and electron modes |
| Transport | Constant, Sutherland and power-law viscosity, with a constant Prandtl number; Wilke and Gupta-Yos mixture rules; binary and mixture-averaged diffusion coefficients |
| Reactions | Arrhenius rates with third bodies, equilibrium constants from the thermodynamic fits, and controlling temperatures for multi-temperature models |
| Energy exchange | Vibration-translation (Landau-Teller with Millikan-White times and high-temperature corrections), chemistry-vibration coupling, electron-heavy elastic exchange, vibration-electron exchange, electron-impact ionisation |

Shipped data, in `include/canary/data/`, generated from `database/data/`:

| Data | Content |
| --- | --- |
| `air5` | N2, O2, N, O, NO: thermodynamic fits, electronic levels, viscosity and conductivity fits, collision integrals |
| `air11` | air5 and N+, O+, N2+, O2+, NO+, e- |
| `air5_park`, `air5_kim_jo`, `air11_kim_jo` | Reaction mechanisms |
| `air5_park_exchange`, `air11_kim_jo_exchange` | Energy-exchange models |

Each data file cites the publications of its values in its `reference`
fields. Other species and mechanisms are read from JSON files with the
optional database library (see [database/README.md](database/README.md)).

Every physical value has a type that states its kind and unit system
(`Density<double, si>`, `Temperature<double, si>`), so a pressure cannot be
passed where a temperature is expected. Models can also work in a scaled unit
system defined by reference values.

## Requirements

- A C++20 compiler. Tested with GCC 16.2 and Clang 23.1.
- [mdspan](https://github.com/kokkos/mdspan) (the Kokkos reference
  implementation).
- CMake 3.29 or newer to build the tests or to install.
- The optional database library needs a C++23 compiler and glaze to build.

## Example

The program builds the five-species air model and the Park mechanism, then
evaluates one cell from its density, temperature and mass fractions.

```cpp
#include <array>
#include <canary/data/air5.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/reactions.hpp>
#include <canary/thermo.hpp>
#include <cstddef>
#include <cstdio>
#include <mdspan/mdspan.hpp>

int main()
{
    using namespace canary;

    // Host: build the model and the mechanism once.
    const UnitScales<double, si> scales {};
    const auto gas = make_air5(scales).model;
    const auto mechanism = make_air5_park(scales);

    // Per cell: mass fractions of N2, O2, N, O, NO, and two state variables.
    const std::array<double, 5> y {0.70, 0.05, 0.05, 0.15, 0.05};
    const Density<double, si> rho {0.01};
    const Temperature<double, si> temp {6000.0};

    const Kokkos::mdspan<const double, Kokkos::extents<std::size_t, 5>> y_view(y.data());
    const auto mix = mixture(gas, as_kind<kind::mass_fraction>(y_view));
    const auto closure = thermo_state(gas, mix, rho, temp);
    if (closure.status != Status::Valid) { return 1; }
    const auto& state = closure.state;
    const auto props = properties(gas, mix, state);
    const auto rates = production_rates(mechanism, gas, mix, rho, temp);

    std::printf("p = %.6g Pa, e = %.6g J/kg, a = %.6g m/s\n",
                state.p.value(),
                state.e.value(),
                props.acoustics.sound_speed.value());
    std::printf("production rate of NO = %.6g kg/(m^3 s)\n", rates[4].value());
    return 0;
}
```

Output:

```
p = 20534.4 Pa, e = 9.55087e+06 J/kg, a = 1649.74 m/s
production rate of NO = 26.2878 kg/(m^3 s)
```

`thermo_state` accepts any supported pair of state variables, for example
density and specific energy, or pressure and temperature. When the pair
requires a temperature solve, check `closure.status`: a solve that did not
converge returns a status other than `Status::Valid` and a diagnostic, and
never throws.

## Using Canary in a solver

**With CMake**, after installing Canary (see [cmake/README.md](cmake/README.md)):

```cmake
find_package(canary 2026.11 CONFIG REQUIRED)
target_link_libraries(my_solver PRIVATE canary::canary)
```

A requested version must equal the installed version. **Without CMake**, add
`include/` and the mdspan headers to the include path.

**Host and per-cell code.** Build models, mechanisms and exchange models once
on the host with the `make_...` functions. They validate their input and throw
`ParameterError` on invalid data. The built objects are trivially copyable,
so a solver can copy them into a GPU kernel. Per-cell
functions do not allocate and do not throw.

**Which headers to include.** `<canary/canary.hpp>` includes everything,
including every generated data set and the host builders, and is convenient
for host code. In a file that defines a GPU kernel or a hot loop, include only
what the kernel calls, which shortens compilation:

| Need | Header |
| --- | --- |
| `thermo_state`, `properties` | `<canary/thermo.hpp>` and the model header, for example `<canary/models/thermally_perfect.hpp>` |
| Production rates | `<canary/reactions.hpp>` |
| Exchange source terms | `<canary/energy_exchange.hpp>` |
| Transport coefficients | `<canary/transport/mixture.hpp>` for the Wilke and Gupta-Yos rules; `<canary/transport/sutherland.hpp>`, `power_law.hpp`, `constant.hpp`, `prandtl.hpp` for the simple laws |
| A shipped data set and its builder | `<canary/data/air5.hpp>` and the others in `canary/data/` |

**Consumer settings.** Define these before including any Canary header
(see `include/canary/config.hpp`):

| Macro | Purpose | Default |
| --- | --- | --- |
| `CANARY_HOST_DEVICE` | Annotation of every per-cell function, for example `AMREX_GPU_HOST_DEVICE` or `KOKKOS_FUNCTION` | empty |
| `CANARY_FORCE_INLINE` | Forced-inline specifier of hot per-cell functions | compiler attribute |
| `CANARY_ASSERT` | Precondition check in per-cell code | `assert` |
| `CANARY_MAX_SPECIES` | Species capacity of generic mixture models | 16 |

## Precision and devices

- **Double precision** (`double` state values, and dual numbers of `double`
  for automatic differentiation) is supported and is what the tests require.
- **Single precision** (`float` state values) compiles and is experimental.
  Temperature solves can end with `Status::NotConverged` or a diagnostic with
  `ClosureReason::PrecisionLimit`, and results that subtract nearly equal numbers lose most of
  their digits. Check the status of every closure.
- **CUDA**: every public per-cell function is compiled as a CUDA kernel and
  run on a device by the device tests (CUDA 13.4, C++20,
  `--expt-relaxed-constexpr`). The device results in double precision agree
  with the host results to a relative difference of 1e-9; see
  [tests/device/README.md](tests/device/README.md).
- **HIP**: the same checks exist for HIP and have not been compiled or run.
- With NVIDIA's compiler, building a model in a `.cu` file prints the warning
  `'long double' is treated as 'double' in device code`. It comes from the
  host-only construction code, which uses `long double`; that code does not
  run on the device and the warning does not affect results.

## Building and testing

```sh
cmake --workflow --preset release      # configure, build and test with GCC
cmake --workflow --preset debug-clang  # the same with Clang, debug
cmake --workflow --preset cuda-runtime # host tests and the CUDA device tests
```

`CMakePresets.json` lists every preset. The host tests compare Canary with
reference values from the cited publications, with independent high-precision
calculations and with the GDTk gas package.

## Known limitations

- Single precision is experimental (see above).
- HIP is unverified.
- The NASA9 fits of adjacent temperature intervals are joined by a smooth
  blend whose widths (400 K to 1000 K) are a choice of this library, not
  published data.
- The electron-vibration relaxation time is a fit of Lee's 1985 expression;
  Lee recommends a modified rate form that Canary does not implement.
- Releases make no compatibility promise to each other: an interface can
  change between calendar versions.

## Licence

Copyright (C) 2026 Shahzeb Imran and the canary contributors.

Canary is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version, along with the Canary Cerberus Exception. See
[LICENSE](LICENSE) and [LICENSE.CERBERUS-EXCEPTION](LICENSE.CERBERUS-EXCEPTION).

Canary is distributed in the hope that it will be useful, but without any
warranty; without even the implied warranty of merchantability or fitness for
a particular purpose.

The exception concerns [Cerberus](https://github.com/PlasmaSimUQ/cerberus)
only. Cerberus, and works based on it, may be combined with Canary without
the AGPL extending to their own code, on condition that they state that they
use Canary and make the source of the Canary version they use available to
everyone who receives the program or uses it over a network. Canary itself,
modified or not, stays under the AGPL. Any other software that uses Canary is
subject to the AGPL as a whole. The licence files are authoritative; this
paragraph is a summary.
