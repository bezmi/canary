// A view of partial densities is not accepted where mass fractions are.
#include <array>
#include <cstddef>
#include <canary/species.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
using View = Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>;
constexpr canary::SpeciesData<double, canary::scaled, 1> species {
  .n = 1,
  .mass = {canary::ParticleMass<double, canary::scaled> {1.0}},
  .charge = {canary::Charge<double, canary::scaled> {0.0}},
};
const std::array<double, 1> values {1.0};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::mean_particle_mass(
  species, canary::as_kind<canary::kind::density, canary::scaled>(View(values.data(), 1)));
#endif
auto good = canary::mean_particle_mass(
  species, canary::as_kind<canary::kind::mass_fraction>(View(values.data(), 1)));
