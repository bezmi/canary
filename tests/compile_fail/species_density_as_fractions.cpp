// A view of partial densities is not accepted where mass fractions are.
#include <array>
#include <cstddef>
#include <gaslib/species.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
using View = Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>;
constexpr gaslib::SpeciesData<double, gaslib::scaled, 1> species {
  .n = 1,
  .mass = {gaslib::ParticleMass<double, gaslib::scaled> {1.0}},
  .charge = {gaslib::Charge<double, gaslib::scaled> {0.0}},
};
const std::array<double, 1> values {1.0};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::mean_particle_mass(
  species, gaslib::as_kind<gaslib::kind::density, gaslib::scaled>(View(values.data(), 1)));
#endif
auto good = gaslib::mean_particle_mass(
  species, gaslib::as_kind<gaslib::kind::mass_fraction>(View(values.data(), 1)));
