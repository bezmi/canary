// Closure throughput from (rho, rho_e) over 10^6 cells: the one-species
// model through gaslib and on plain doubles, mixtures of 5 and 11 species
// with the exact capacity and with capacity 16, and the thermally perfect
// model for N2 alone and for air (N2 and O2 in the five-species model),
// without and with a temperature guess equal to the previous cell's
// temperature, the two-temperature model for air from (rho, rho_e,
// rho_e_ve), and the four-temperature model of five-species air and the
// three-temperature model of weakly ionised eleven-species air from
// (rho, rho_e, rho_e_m). The last cases are the mass production rates of the Kim-Jo
// five-species air mechanism from (rho, T), and the mass production rates
// and vibrational-electronic energy source of the Park five-species air
// mechanism and energy exchange on the two-temperature model from
// (rho, T, T_ve), with mass fractions that vary from cell to cell. Then the
// transport coefficients and the mixture-averaged diffusion coefficients of
// five-species air by the Wilke model on the thermally perfect model and by
// the Gupta-Yos model on the two-temperature model, from states with mass
// fractions that vary from cell to cell. Each case runs several times; the
// fastest run is reported in nanoseconds per cell.

#include <array>
#include <chrono>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <gaslib/gaslib.hpp>
#include <mdspan/mdspan.hpp>
#include <ratio>
#include <tuple>
#include <vector>

namespace
{

using gaslib::si;
constexpr std::size_t cells = 1'000'000;
constexpr int repeats = 7;

// Fastest of the repeated runs of body(), in nanoseconds per cell. The sum of
// the results is printed so the work is not removed.
template <typename Body> double time_per_cell(const char* name, Body body)
{
    double best = 1.0e300;
    double sink = 0;
    for (int r = 0; r < repeats; ++r) {
        const auto start = std::chrono::steady_clock::now();
        sink += body();
        const auto stop = std::chrono::steady_clock::now();
        const double ns = std::chrono::duration<double, std::nano>(stop - start).count() / cells;
        best = ns < best ? ns : best;
    }
    std::printf("%-40s %8.3f ns/cell   (checksum %.6e)\n", name, best, sink);
    return best;
}

struct Cells {
    std::vector<double> rho;
    std::vector<double> rho_e;
};

Cells make_cells()
{
    Cells c {.rho = std::vector<double>(cells), .rho_e = std::vector<double>(cells)};
    for (std::size_t i = 0; i < cells; ++i) {
        c.rho[i] = 1.0 + (0.5 * static_cast<double>(i % 97) / 97.0);
        c.rho_e[i] = 2.5e5 * (1.0 + (0.3 * static_cast<double>(i % 89) / 89.0));
    }
    return c;
}

// Mass fractions of n species per cell, contiguous per cell.
std::vector<double> make_fractions(std::size_t n)
{
    std::vector<double> y(cells * n);
    for (std::size_t i = 0; i < cells; ++i) {
        double sum = 0;
        for (std::size_t s = 0; s < n; ++s) {
            y[(i * n) + s] = 1.0 + static_cast<double>((i + (3 * s)) % 7);
            sum += y[(i * n) + s];
        }
        for (std::size_t s = 0; s < n; ++s) { y[(i * n) + s] /= sum; }
    }
    return y;
}

template <std::size_t Capacity>
gaslib::CaloricallyPerfectMixture<double, si, Capacity> make_mixture(std::size_t n)
{
    std::vector<gaslib::ParticleMass<double, si>> m;
    std::vector<gaslib::Charge<double, si>> q;
    std::vector<gaslib::HeatCapacityRatio<double>> g;
    for (std::size_t s = 0; s < n; ++s) {
        m.emplace_back(4.65e-26 * (1.0 + (0.1 * static_cast<double>(s))));
        q.emplace_back(0.0);
        g.emplace_back(s % 2 == 0 ? 1.4 : 5.0 / 3.0);
    }
    return gaslib::make_calorically_perfect_mixture<Capacity>(
             m,
             q,
             gaslib::boltzmann_constant(gaslib::UnitScales<double, si> {}),
             g)
      .model;
}

template <std::size_t Capacity> double mixture_case(const char* name, const Cells& c, std::size_t n)
{
    const auto gas = make_mixture<Capacity>(n);
    const std::vector<double> y = make_fractions(n);
    return time_per_cell(name, [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data() + (i * n),
                                                                             n));
            const auto closure =
              gaslib::thermo_state(gas,
                                   fractions,
                                   gaslib::Density<double, si> {c.rho[i]},
                                   gaslib::EnergyDensity<double, si> {c.rho_e[i]});
            sum += closure.state.p.value() + closure.state.temp.value();
        }
        return sum;
    });
}

// Thermally perfect model G with the same mass fractions y in every cell.
// With guess, each cell starts from the temperature of the previous cell
// (300 K for the first).
template <typename G, std::size_t N>
double thermally_perfect_case(
  const char* name, const Cells& c, const G& gas, const std::array<double, N>& y, bool guess)
{
    const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
      Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), N));
    return time_per_cell(name, [&] {
        double sum = 0;
        gaslib::Temperature<double, si> previous {300.0};
        for (std::size_t i = 0; i < cells; ++i) {
            const gaslib::Density<double, si> rho {c.rho[i]};
            const gaslib::EnergyDensity<double, si> rho_e {c.rho_e[i]};
            const auto closure =
              guess ? gaslib::thermo_state(gas,
                                           fractions,
                                           rho,
                                           rho_e,
                                           gaslib::TemperatureGuess<double, si> {previous})
                    : gaslib::thermo_state(gas, fractions, rho, rho_e);
            previous = closure.state.temp;
            sum += closure.state.p.value() + closure.state.temp.value();
        }
        return sum;
    });
}

// Two-temperature model G with the same mass fractions y in every cell, from
// (rho, rho_e, rho_e_ve). The inputs of each cell are those of the state with
// its density, T_ve = 0.8 T, and T the temperature of the single-temperature
// model of the same species at the cell's (rho, rho_e).
template <typename G, std::size_t N>
double two_temperature_case(const char* name,
                            const Cells& c,
                            const G& gas,
                            const std::array<double, N>& y)
{
    const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
      Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), N));
    const auto mix = gaslib::mixture(gas, fractions);
    std::vector<double> rho_e(cells);
    std::vector<double> rho_e_ve(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        const gaslib::Density<double, si> rho {c.rho[i]};
        const auto temp = gaslib::thermo_state(gas.thermally_perfect,
                                               mix,
                                               rho,
                                               gaslib::EnergyDensity<double, si> {c.rho_e[i]})
                            .state.temp;
        const std::array<double, 1> t_ve {0.8 * temp.value()};
        const auto state =
          gaslib::thermo_state(
            gas,
            mix,
            rho,
            temp,
            gaslib::as_kind<gaslib::kind::mode_temperature, si>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(t_ve.data(), 1)))
            .state;
        rho_e[i] = c.rho[i] * state.e.value();
        rho_e_ve[i] = c.rho[i] * state.mode_e[0].value();
    }
    return time_per_cell(name, [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto modes = gaslib::as_kind<gaslib::kind::mode_energy_density, si>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(&rho_e_ve[i], 1));
            const auto closure = gaslib::thermo_state(gas,
                                                      fractions,
                                                      gaslib::Density<double, si> {c.rho[i]},
                                                      gaslib::EnergyDensity<double, si> {rho_e[i]},
                                                      modes);
            sum += closure.state.p.value() + closure.state.temp.value() +
                   closure.state.mode_temp[0].value();
        }
        return sum;
    });
}

// Multi-temperature model G with the same mass fractions y in every cell,
// from (rho, rho_e, rho_e_m). The inputs of each cell are those of the state
// with its density, the temperature temp[i] and every T_m = 0.8 T.
template <typename G, std::size_t N>
double multi_temperature_case(const char* name,
                              const Cells& c,
                              const G& gas,
                              const std::array<double, N>& y,
                              const std::vector<double>& temp)
{
    constexpr std::size_t modes = std::tuple_size_v<decltype(gaslib::state_t<G, double> {}.mode_e)>;
    const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
      Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), N));
    const auto mix = gaslib::mixture(gas, fractions);
    std::vector<double> rho_e(cells);
    std::vector<double> rho_e_m(cells * modes);
    for (std::size_t i = 0; i < cells; ++i) {
        std::array<double, modes> t_m {};
        t_m.fill(0.8 * temp[i]);
        const auto state =
          gaslib::thermo_state(
            gas,
            mix,
            gaslib::Density<double, si> {c.rho[i]},
            gaslib::Temperature<double, si> {temp[i]},
            gaslib::as_kind<gaslib::kind::mode_temperature, si>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(t_m.data(), modes)))
            .state;
        rho_e[i] = c.rho[i] * state.e.value();
        for (std::size_t m = 0; m < modes; ++m) {
            rho_e_m[(i * modes) + m] = c.rho[i] * state.mode_e[m].value();
        }
    }
    return time_per_cell(name, [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto energies = gaslib::as_kind<gaslib::kind::mode_energy_density, si>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(&rho_e_m[i * modes],
                                                                             modes));
            const auto closure = gaslib::thermo_state(gas,
                                                      fractions,
                                                      gaslib::Density<double, si> {c.rho[i]},
                                                      gaslib::EnergyDensity<double, si> {rho_e[i]},
                                                      energies);
            sum += closure.state.p.value() + closure.state.temp.value() +
                   closure.state.mode_temp[modes - 1].value();
        }
        return sum;
    });
}

// Mass production rates of mechanism mech with thermo model gas, per cell
// from the cell's density, a temperature between 3000 K and 15000 K and the
// mass fractions of make_fractions.
template <typename M, typename G>
double production_case(const char* name, const Cells& c, const M& mech, const G& gas)
{
    constexpr std::size_t n = 5;
    const std::vector<double> y = make_fractions(n);
    std::vector<double> temp(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        temp[i] = 3000.0 + (12000.0 * static_cast<double>(i % 101) / 101.0);
    }
    return time_per_cell(name, [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data() + (i * n),
                                                                             n));
            const auto omega =
              gaslib::production_rates(mech,
                                       gas,
                                       gaslib::mixture(gas, fractions),
                                       gaslib::Density<double, si> {1.0e-3 * c.rho[i]},
                                       gaslib::Temperature<double, si> {temp[i]});
            sum += omega[2].value() + omega[3].value();
        }
        return sum;
    });
}

// two_temperature_sources of mechanism mech and exchange with the
// two-temperature model gas, per cell from the cell's density, the
// temperatures T of production_case and T_ve = T / 2 + 150 K, and the mass
// fractions of make_fractions.
template <typename M, typename X, typename G>
double
  sources_case(const char* name, const Cells& c, const M& mech, const X& exchange, const G& gas)
{
    constexpr std::size_t n = 5;
    const std::vector<double> y = make_fractions(n);
    std::vector<double> temp(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        temp[i] = 3000.0 + (12000.0 * static_cast<double>(i % 101) / 101.0);
    }
    return time_per_cell(name, [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data() + (i * n),
                                                                             n));
            const auto sources = gaslib::two_temperature_sources(
              mech,
              exchange,
              gas,
              gaslib::mixture(gas, fractions),
              gaslib::Density<double, si> {1.0e-3 * c.rho[i]},
              gaslib::Temperature<double, si> {temp[i]},
              gaslib::Temperature<double, si> {(0.5 * temp[i]) + 150.0});
            sum += sources.species[2].value() + sources.species[3].value() +
                   sources.vibrational_electronic.value();
        }
        return sum;
    });
}

// States of gas per cell from (p, T) and, for a model with modes,
// T_ve = T / 2 + 150 K, with the mass fractions y of make_fractions, the
// temperatures of production_case and p = 1e4 Pa times the cell's density.
template <typename G>
std::vector<gaslib::state_t<G, double>>
  make_states(const Cells& c, const G& gas, const std::vector<double>& y)
{
    constexpr std::size_t n = 5;
    std::vector<gaslib::state_t<G, double>> states(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
          Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data() + (i * n), n));
        const auto mix = gaslib::mixture(gas, fractions);
        const double temp = 3000.0 + (12000.0 * static_cast<double>(i % 101) / 101.0);
        const gaslib::Pressure<double, si> p {1.0e4 * c.rho[i]};
        if constexpr (std::same_as<gaslib::state_t<G, double>, gaslib::ThermoState<double, si>>) {
            states[i] =
              gaslib::thermo_state(gas, mix, p, gaslib::Temperature<double, si> {temp}).state;
        } else {
            const std::array<double, 1> t_ve {(0.5 * temp) + 150.0};
            states[i] =
              gaslib::thermo_state(
                gas,
                mix,
                p,
                gaslib::Temperature<double, si> {temp},
                gaslib::as_kind<gaslib::kind::mode_temperature, si>(
                  Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(t_ve.data(), 1)))
                .state;
        }
    }
    return states;
}

// Transport coefficients of model with gas per cell, at the states of
// make_states, with the mixture value computed from the cell's mass
// fractions. With diffusion, the mixture-averaged diffusion coefficients
// instead.
template <typename M, typename G>
double
  transport_case(const char* name, const Cells& c, const M& model, const G& gas, bool diffusion)
{
    constexpr std::size_t n = 5;
    const std::vector<double> y = make_fractions(n);
    const auto states = make_states(c, gas, y);
    return time_per_cell(name, [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto fractions = gaslib::as_kind<gaslib::kind::mass_fraction>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data() + (i * n),
                                                                             n));
            const auto mix = gaslib::mixture(gas, fractions);
            if (diffusion) {
                const auto d = gaslib::mixture_diffusion_coefficients(model, gas, mix, states[i]);
                sum += d[0].value() + d[2].value();
            } else {
                const auto coefficients = gaslib::transport(model, gas, mix, states[i]);
                sum += coefficients.viscosity.value() + coefficients.conductivity.value();
            }
        }
        return sum;
    });
}

}  // namespace

int main()
{
    const Cells c = make_cells();
    const auto air =
      gaslib::make_calorically_perfect(gaslib::SpecificGasConstant<double, si> {287.0},
                                       gaslib::HeatCapacityRatio {1.4});

    time_per_cell("one species, plain doubles", [&] {
        const double gamma = air.gamma.value();
        const double r = air.gas_constant.value();
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const double p = (gamma - 1.0) * c.rho_e[i];
            const double temp = p / (c.rho[i] * r);
            sum += p + temp;
        }
        return sum;
    });
    time_per_cell("one species, gaslib", [&] {
        double sum = 0;
        for (std::size_t i = 0; i < cells; ++i) {
            const auto closure =
              gaslib::thermo_state(air,
                                   gaslib::Density<double, si> {c.rho[i]},
                                   gaslib::EnergyDensity<double, si> {c.rho_e[i]});
            sum += closure.state.p.value() + closure.state.temp.value();
        }
        return sum;
    });
    mixture_case<5>("5 species, capacity 5", c, 5);
    mixture_case<16>("5 species, capacity 16", c, 5);
    mixture_case<11>("11 species, capacity 11", c, 11);
    mixture_case<16>("11 species, capacity 16", c, 11);

    const gaslib::UnitScales<double, si> si_scales {};
    const auto n2 =
      gaslib::make_thermally_perfect_mixture(std::array {gaslib::air5_species[0]}, si_scales).model;
    const std::array<double, 1> n2_y {1.0};
    thermally_perfect_case("thermally perfect N2", c, n2, n2_y, false);
    thermally_perfect_case("thermally perfect N2, guess", c, n2, n2_y, true);
    const auto air5 = gaslib::make_air5<double>(si_scales).model;
    const std::array<double, 5> air_y {0.767, 0.233, 0.0, 0.0, 0.0};
    thermally_perfect_case("thermally perfect air (5 species)", c, air5, air_y, false);
    thermally_perfect_case("thermally perfect air (5 species), guess", c, air5, air_y, true);
    const auto air5_two_temperature = gaslib::make_air5_two_temperature<double>(si_scales).model;
    two_temperature_case("two-temperature air (5 species)", c, air5_two_temperature, air_y);
    // Temperatures of the two-temperature case: those of the
    // single-temperature model at the cells' (rho, rho_e).
    std::vector<double> air_temp(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        air_temp[i] =
          gaslib::thermo_state(
            air5,
            gaslib::as_kind<gaslib::kind::mass_fraction>(
              Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(air_y.data(), 5)),
            gaslib::Density<double, si> {c.rho[i]},
            gaslib::EnergyDensity<double, si> {c.rho_e[i]})
            .state.temp.value();
    }
    multi_temperature_case("four-temperature air (5 species)",
                           c,
                           gaslib::make_air5_four_temperature<double>(si_scales).model,
                           air_y,
                           air_temp);
    // Weakly ionised air, charge neutral, at 3000 K to 15000 K.
    const std::array<double, 11> ionised_y {
      0.4999995398459613,
      0.049999953984596134,
      0.19999981593838453,
      0.19999981593838453,
      0.01999998159383845,
      0.009999990796919225,
      0.009999990796919225,
      0.003999996318767691,
      0.0019999981593838453,
      0.003999996318767691,
      9.203080774345904e-07,
    };
    std::vector<double> ionised_temp(cells);
    for (std::size_t i = 0; i < cells; ++i) {
        ionised_temp[i] = 3000.0 + (12000.0 * static_cast<double>(i % 101) / 101.0);
    }
    multi_temperature_case("three-temperature air (11 species)",
                           c,
                           gaslib::make_air11_three_temperature<double>(si_scales).model,
                           ionised_y,
                           ionised_temp);
    const auto air5_built = gaslib::make_air5<double>(si_scales);
    const auto kim_jo = gaslib::make_air5_kim_jo(air5_built.species, si_scales);
    production_case("Kim-Jo air (13 reactions) production", c, kim_jo, air5);
    const auto air5_two_temperature_built = gaslib::make_air5_two_temperature<double>(si_scales);
    const auto park = gaslib::make_air5_park(air5_two_temperature_built.species, si_scales);
    const auto park_exchange =
      gaslib::make_air5_park_exchange(air5_two_temperature_built.species, park, si_scales);
    sources_case("Park air (5 reactions) 2T sources",
                 c,
                 park,
                 park_exchange,
                 air5_two_temperature_built.model);
    const auto wilke = gaslib::make_air5_wilke_transport(air5_built.species, si_scales);
    const auto gupta_yos =
      gaslib::make_air5_gupta_yos_transport(air5_two_temperature_built.species, si_scales);
    transport_case("Wilke air (5 species) transport", c, wilke, air5, false);
    transport_case("Gupta-Yos air (5 species) 2T transport",
                   c,
                   gupta_yos,
                   air5_two_temperature_built.model,
                   false);
    transport_case("Wilke air (5 species) diffusion", c, wilke, air5, true);
    transport_case("Gupta-Yos air (5 species) 2T diffusion",
                   c,
                   gupta_yos,
                   air5_two_temperature_built.model,
                   true);
}
