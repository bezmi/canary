// Closure throughput from (rho, rho_e) over 10^6 cells: the one-species
// model through gaslib and on plain doubles, and mixtures of 5 and 11
// species with the exact capacity and with capacity 16. Each case runs
// several times; the fastest run is reported in nanoseconds per cell.

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <gaslib/gaslib.hpp>
#include <mdspan/mdspan.hpp>
#include <ratio>
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
}
