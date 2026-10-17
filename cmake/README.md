# Installing and consuming Canary

Canary requires CMake 3.29 or newer and a C++20 compiler. Install mdspan
separately, then configure Canary with that installation in `CMAKE_PREFIX_PATH`:

```sh
cmake -S . -B build/install -G Ninja \
  -DCMAKE_INSTALL_PREFIX=/path/to/canary-install \
  -DCMAKE_PREFIX_PATH=/path/to/mdspan-install \
  -DCANARY_FIND_MDSPAN=ON \
  -DCANARY_BUILD_TESTS=OFF \
  -DCANARY_BUILD_BENCHMARKS=OFF \
  -DCANARY_BUILD_DATABASE=OFF
cmake --install build/install
```

Set `CANARY_BUILD_DATABASE=ON` to include the static database library. This
requires a C++23 compiler and glaze while building Canary. Build before installing:

```sh
cmake --build build/install --target canary_database --parallel 2
cmake --install build/install
```

The installed database interface supports C++20 consumers and does not require
glaze. Canary does not install mdspan or other build dependencies.

In the consuming project's `CMakeLists.txt`:

```cmake
find_package(canary CONFIG REQUIRED)
target_link_libraries(my_solver PRIVATE canary::canary)
# Available when Canary was installed with CANARY_BUILD_DATABASE=ON:
target_link_libraries(my_database_tool PRIVATE canary::database)
```

Configure that project with both installation prefixes:

```sh
cmake -S /path/to/consumer -B /path/to/consumer-build \
  '-DCMAKE_PREFIX_PATH=/path/to/canary-install;/path/to/mdspan-install'
cmake --build /path/to/consumer-build --parallel 2
```

`find_package(canary)` requires the separately installed `mdspan::mdspan` target.
The Canary installation can be moved to another prefix; update the consumer's
`CMAKE_PREFIX_PATH` to the new location. No package version is declared yet.

## Numerical device checks

CUDA and HIP compilation checks of the numerical operations are opt-in test
targets; see [tests/device/README.md](../tests/device/README.md). They are not
installed and add no requirement to consumers.
