include(CPM)

# Both branches provide the target mdspan::mdspan. A project that also uses
# mdspan through CPM must pin the same commit: CPM reuses the first package
# added under a name and does not compare GIT_TAG.
if(CANARY_FIND_MDSPAN)
  find_package(mdspan CONFIG REQUIRED)
else()
  CPMAddPackage(
      NAME mdspan
      GITHUB_REPOSITORY kokkos/mdspan
      GIT_TAG c8290f0f1977135450cef6cd460dcd925b23a38f
      SYSTEM YES
      EXCLUDE_FROM_ALL YES
      OPTIONS
          "MDSPAN_ENABLE_TESTS OFF"
          "MDSPAN_ENABLE_EXAMPLES OFF"
          "MDSPAN_ENABLE_BENCHMARKS OFF"
  )
endif()

if(NOT TARGET mdspan::mdspan)
  message(FATAL_ERROR "canary requires the CMake target mdspan::mdspan, which the mdspan package did not define")
endif()
