#ifndef CANARY_MODELS_MULTI_TEMPERATURE_HPP
#define CANARY_MODELS_MULTI_TEMPERATURE_HPP

#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/detail/solve.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <string_view>

namespace canary
{

// Host input: the electronic levels of one species, named by name. Level k
// (0 <= k < levels) has degeneracy g_k = degeneracies[k] and level
// temperature theta_k = E_k / k_B = temperatures[k] in K. Level 0 is the
// ground level, with theta_0 = 0. Entries past levels are ignored.
template <std::size_t L>
    requires(L > 0)
struct ElectronicLevelsInput {
    const char* name;
    int levels;
    std::array<int, L> degeneracies;
    std::array<double, L> temperatures;
};

// The kinds of energy that an energy mode can hold for one species:
// vibration and electronic excitation of a heavy species, and the
// translational energy of a species with the electron role.
enum class ModeComponent : std::uint8_t { vibration, electronic, electron };

// Host input: one component of an energy mode, the species by name and the
// kind of energy.
struct ModeComponentInput {
    const char* species;
    ModeComponent component;
};

// Host input: one energy mode, made of components[0 .. count). Entries past
// count are ignored.
template <std::size_t C>
    requires(C > 0)
struct EnergyModeInput {
    int count;
    std::array<ModeComponentInput, C> components;
};

// Electronic levels of one species in the unit system U: the first levels
// entries of degeneracies (g_k) and temperatures (theta_k, in the
// temperature unit of U), with theta_0 = 0.
template <std::floating_point P, UnitSystem U, std::size_t L>
    requires(L > 0)
struct ElectronicLevels {
    int levels;
    std::array<P, L> degeneracies;
    std::array<Temperature<P, U>, L> temperatures;
};

// The components of one species that one energy mode holds.
struct ModeComponentSet {
    bool vibration;
    bool electronic;
    bool electron;
};

// Ideal gas mixture of N species with a translational-rotational
// temperature T (State::temp) and M energy modes, mode m with its own
// temperature T_m (State::mode_temp[m]) and energy per unit mass e_m
// (State::mode_e[m]). SpeciesMetadataView supplies the heavy/electron role
// during construction. For species i, with gas constant R_i,
// cp_i(T) and h_i(T) per unit mass from its NASA 9-coefficient curve, T_ref
// the reference temperature and d_i its rotational degrees of freedom
// (0 for an atom and for free electrons, 2 for a linear molecule, 3 for a nonlinear
// molecule):
// - cp_tr,i = (5/2 + d_i/2) R_i, cv_tr,i = cp_tr,i - R_i and
//   dh_f,i = h_i(T_ref);
// - eps_i(T) = R_i sum_k g_k theta_k exp(-theta_k / T) /
//   sum_k g_k exp(-theta_k / T) is the electronic energy of its levels (0
//   without levels), and c_el,i(T) = d eps_i / dT.
// The components and their energies per unit mass at the mode temperature
// T_m are:
// - vibration of a heavy species with d_i > 0:
//   e_vib,i = h_i(T_m) - cp_tr,i (T_m - T_ref) - dh_f,i
//             - (eps_i(T_m) - eps_i(T_ref)),
//   cv_vib,i = cp_i(T_m) - cp_tr,i - c_el,i(T_m);
// - electronic of a heavy species: e_el,i = eps_i(T_m) - eps_i(T_ref),
//   cv_el,i = c_el,i(T_m);
// - electron of e-: e_e = h_e(T_m) - R_e T_m, cv_e = cp_e(T_m) - R_e.
// Each component is in at most one mode; a component in no mode is frozen
// and contributes no energy. With mass fractions y_i:
// - e_tr = sum over heavy i of y_i (cp_tr,i (T - T_ref) + dh_f,i - R_i T);
// - e_m = sum over the components (i, c) of mode m of y_i e_c,i(T_m);
// - e = e_tr + sum_m e_m;
// - p = rho (R_h T + y_e R_e T_e), with R_h = sum over heavy i of y_i R_i
//   and T_e the temperature of the mode that holds e-;
// - cv = sum over heavy i of y_i cv_tr,i + sum_m cv_m, with cv_m the sum of
//   y_i cv_c,i(T_m) over the components of mode m, and cp = cv + R_h.
// A heavy species whose vibration and electronic components are both in
// modes has its NASA 9-coefficient energy when all temperatures are equal.
// All N species are used. The fields are in unit system U and are derived
// by make_multi_temperature_mixture:
// - thermally_perfect: the single-temperature model of the same species.
//   Its gas constants, curves and reference temperature define the species
//   above; its temperature limits, bracket width, tolerance and iteration
//   limit are those of every temperature solve; its default initial guess
//   starts the solve at thermal equilibrium;
// - translational_rotational_cp: cp_tr,i (used for heavy species);
// - formation_enthalpy: dh_f,i (used for heavy species);
// - electronic_levels: the levels of each species, and
//   reference_electronic_energy: eps_i(T_ref);
// - mode_components[m][i]: the components of species i in mode m;
// - electron: the optional free-electron species index; electron_mode: its
//   energy mode, or -1 without free electrons.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S = 4,
          std::size_t L = 8>
    requires(N > 0 && M > 0 && S > 0 && L > 0)
struct MultiTemperatureMixture {
    using parameter_type = P;
    using units = U;
    template <Scalar T> using State = MultiTemperatureState<T, U, M>;
    template <Scalar T> using Mixture = NormalisedComposition<T, U, N>;
    static constexpr std::size_t mode_count = M;

    ThermallyPerfectMixture<P, U, N, S> thermally_perfect;
    std::array<SpecificHeat<P, U>, N> translational_rotational_cp;
    std::array<SpecificEnthalpy<P, U>, N> formation_enthalpy;
    std::array<ElectronicLevels<P, U, L>, N> electronic_levels;
    std::array<SpecificEnergy<P, U>, N> reference_electronic_energy;
    std::array<std::array<ModeComponentSet, N>, M> mode_components;
    OptionalSpeciesIndex electron;
    int electron_mode;
};

namespace detail
{

// e = eps(t) and cv = c_el(t) of a species with the given levels and gas
// constant, at temperature t (model units), with plain values of type W.
// With theta_0 = 0 every weight w_k = g_k exp(-theta_k / t) is at most g_k
// and their sum at least g_0, so no term overflows for any t > 0. The
// derivative is computed as the variance of theta over the weights,
// c_el = R sum_k w_k (theta_k - eps / R)^2 / (t^2 sum_k w_k), which has no
// cancellation.
template <Scalar W, std::floating_point P, UnitSystem U, std::size_t L>
    requires ScalarWith<W, P> && HasExp<W>
CANARY_HOST_DEVICE EnergyAndCv<W> electronic_energy(const ElectronicLevels<P, U, L>& levels,
                                                    SpecificGasConstant<P, U> gas_constant,
                                                    W t)
{
    using V = ScalarTraits<W>::value_type;
    if (levels.levels <= 0) { return {.e = W(V(0)), .cv = W(V(0))}; }
    const auto n = static_cast<std::size_t>(levels.levels);
    const W inv_t = W(V(1)) / t;
    std::array<W, L> weight {};
    weight[0] = W(static_cast<V>(levels.degeneracies[0]));
    W sum = weight[0];
    W first = W(V(0));
    for (std::size_t k = 1; k < n; ++k) {
        const P theta = levels.temperatures[k].value();
        weight[k] = W(levels.degeneracies[k] * math::exp(W(-(theta * inv_t))));
        sum = sum + weight[k];
        first = first + W(weight[k] * theta);
    }
    const W mean = W(first / sum);
    W second = W(V(0));
    for (std::size_t k = 0; k < n; ++k) {
        const W d = W(levels.temperatures[k].value() - mean);
        second = second + W(weight[k] * (d * d));
    }
    const P r = gas_constant.value();
    return {.e = W(r * mean), .cv = W(r * ((second / sum) * (inv_t * inv_t)))};
}

// True when c holds at least one component.
CANARY_HOST_DEVICE constexpr bool holds_any(ModeComponentSet c)
{ return c.vibration || c.electronic || c.electron; }

// The components of species i over all modes.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
CANARY_HOST_DEVICE constexpr ModeComponentSet
  species_components(const MultiTemperatureMixture<P, U, N, M, S, L>& gas, SpeciesIndex index)
{
    const auto i = index.value();
    ModeComponentSet all {.vibration = false, .electronic = false, .electron = false};
    for (std::size_t m = 0; m < M; ++m) {
        const ModeComponentSet c = gas.mode_components[m][i];
        all.vibration = all.vibration || c.vibration;
        all.electronic = all.electronic || c.electronic;
        all.electron = all.electron || c.electron;
    }
    return all;
}

// Energy per unit mass and cv of the components c of species i at
// temperature t (model units), with ln_t = ln t. With both vibration and
// electronic, the eps terms cancel and only the curve is evaluated.
template <Scalar W,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<W, P> && HasExp<W> && HasLog<W>
CANARY_HOST_DEVICE EnergyAndCv<W>
  component_energy(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                   SpeciesIndex index,
                   ModeComponentSet c,
                   W t,
                   W ln_t)
{
    const auto i = index.value();
    using V = ScalarTraits<W>::value_type;
    const auto& single = gas.thermally_perfect;
    const P r = single.gas_constant[i].value();
    if (c.electron) {
        const Nasa9CpH<W> h = nasa9_cp_h(single.curve[i], t, ln_t);
        return {.e = W(r * (h.h_over_r - t)), .cv = W(r * (h.cp_over_r - P(1)))};
    }
    W e = W(V(0));
    W cv = W(V(0));
    if (c.vibration) {
        const P t_ref = single.reference_temperature.value();
        const P cp_tr = gas.translational_rotational_cp[i].value();
        const P h_f = gas.formation_enthalpy[i].value();
        const Nasa9CpH<W> h = nasa9_cp_h(single.curve[i], t, ln_t);
        e = W((r * h.h_over_r) - (cp_tr * (t - t_ref)) - h_f);
        cv = W((r * h.cp_over_r) - cp_tr);
    }
    if (c.vibration != c.electronic) {
        const EnergyAndCv<W> el =
          electronic_energy(gas.electronic_levels[i], single.gas_constant[i], t);
        const W e_el = W(el.e - gas.reference_electronic_energy[i].value());
        if (c.electronic) {
            e = W(e + e_el);
            cv = W(cv + el.cv);
        } else {
            e = W(e - e_el);
            cv = W(cv - el.cv);
        }
    }
    return {.e = e, .cv = cv};
}

// e_m and cv_m of mode m at its temperature t (model units), with plain
// values of type W and one logarithm.
template <Scalar W,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<W, P> && HasLog<W> && HasExp<W>
CANARY_HOST_DEVICE EnergyAndCv<W> mode_energy(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                              const std::array<MassFraction<W>, N>& y,
                                              std::size_t m,
                                              W t)
{
    using V = ScalarTraits<W>::value_type;
    const W ln_t = math::log(t);
    W e = W(V(0));
    W cv = W(V(0));
    for (const SpeciesIndex i : species_indices(0, N)) {
        const ModeComponentSet c = gas.mode_components[m][i.value()];
        if (!holds_any(c)) { continue; }
        const EnergyAndCv<W> r = component_energy(gas, i, c, t, ln_t);
        e = e + W(y[i.value()].value() * r.e);
        cv = cv + W(y[i.value()].value() * r.cv);
    }
    return {.e = e, .cv = cv};
}

// Sums over the heavy species of y_i cp_tr,i, y_i dh_f,i and y_i R_i (R_h),
// and y_e R_e (0 without e-), with plain values of type W.
template <Scalar W> struct HeavySpeciesSums {
    W cp;
    W formation_enthalpy;
    W gas_constant;
    W electron_gas_constant;
};

template <Scalar W,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<W, P>
CANARY_HOST_DEVICE HeavySpeciesSums<W>
  heavy_species_sums(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                     const std::array<MassFraction<W>, N>& y)
{
    using V = ScalarTraits<W>::value_type;
    HeavySpeciesSums<W> sums {
      .cp = W(V(0)),
      .formation_enthalpy = W(V(0)),
      .gas_constant = W(V(0)),
      .electron_gas_constant = W(V(0)),
    };
    for (const SpeciesIndex i : species_indices(0, N)) {
        const P r = gas.thermally_perfect.gas_constant[i.value()].value();
        if (gas.electron == OptionalSpeciesIndex {i}) {
            sums.electron_gas_constant = W(y[i.value()].value() * r);
            continue;
        }
        sums.cp =
          sums.cp + W(y[i.value()].value() * gas.translational_rotational_cp[i.value()].value());
        sums.formation_enthalpy =
          sums.formation_enthalpy +
          W(y[i.value()].value() * gas.formation_enthalpy[i.value()].value());
        sums.gas_constant = sums.gas_constant + W(y[i.value()].value() * r);
    }
    return sums;
}

// e_tr = sum(y cp_tr) (t - T_ref) + sum(y dh_f) - R_h t over the heavy
// species.
template <Scalar W,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<W, P>
CANARY_HOST_DEVICE W translational_energy(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                          const HeavySpeciesSums<W>& sums,
                                          W t)
{
    const P t_ref = gas.thermally_perfect.reference_temperature.value();
    return W((sums.cp * W(t - t_ref)) + sums.formation_enthalpy - (sums.gas_constant * t));
}

// The inverse of translational_energy for t:
// t = (e_tr + T_ref sum(y cp_tr) - sum(y dh_f)) / (sum(y cp_tr) - R_h).
template <Scalar W,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<W, P>
CANARY_HOST_DEVICE W translational_temperature(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                               const HeavySpeciesSums<W>& sums,
                                               W e_tr)
{
    const P t_ref = gas.thermally_perfect.reference_temperature.value();
    return W((e_tr + W(t_ref * sums.cp) - sums.formation_enthalpy) / (sums.cp - sums.gas_constant));
}

// e and cv of the mixture with every mode at temperature t: e_tr(t) plus
// the energy of every component held by a mode, with one logarithm.
template <Scalar W,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<W, P> && HasLog<W> && HasExp<W>
CANARY_HOST_DEVICE EnergyAndCv<W>
  equilibrium_energy(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                     const std::array<MassFraction<W>, N>& y,
                     const HeavySpeciesSums<W>& sums,
                     W t)
{
    const W ln_t = math::log(t);
    W e = translational_energy(gas, sums, t);
    W cv = W(sums.cp - sums.gas_constant);
    for (const SpeciesIndex i : species_indices(0, N)) {
        const ModeComponentSet c = species_components(gas, i);
        if (!holds_any(c)) { continue; }
        const EnergyAndCv<W> r = component_energy(gas, i, c, t, ln_t);
        e = e + W(y[i.value()].value() * r.e);
        cv = cv + W(y[i.value()].value() * r.cv);
    }
    return {.e = e, .cv = cv};
}

// True when some species of mode m has a positive mass fraction (value
// parts).
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          std::floating_point V>
CANARY_HOST_DEVICE bool mode_present(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                     const std::array<MassFraction<V>, N>& y,
                                     std::size_t m)
{
    V sum = V(0);
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (holds_any(gas.mode_components[m][i.value()])) { sum += y[i.value()].value(); }
    }
    return sum > V(0);
}

// The temperature of the mode that holds e-, or temp without e-.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
CANARY_HOST_DEVICE T electron_temperature(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                          const std::array<ModeTemperature<T, U>, M>& mode_temp,
                                          T temp)
{
    return gas.electron_mode < 0 ? temp
                                 : mode_temp[static_cast<std::size_t>(gas.electron_mode)].value();
}

// The elements of a mode range of M modes.
template <std::size_t M, SpeciesRange V>
CANARY_HOST_DEVICE std::array<range_quantity_t<V>, M> mode_values(const V& modes)
{
    // NOLINTNEXTLINE(misc-static-assert): run-time extent
    CANARY_ASSERT(static_cast<std::size_t>(modes.extent(0)) == M);
    std::array<range_quantity_t<V>, M> out {};
    for (std::size_t m = 0; m < M; ++m) { out[m] = modes[m]; }
    return out;
}

// e_m = (rho e_m) / rho of every mode.
template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE std::array<SpecificEnergy<T, U>, M>
  mode_specific_energies(const std::array<ModeEnergyDensity<T, U>, M>& rho_e_m, Density<T, U> rho)
{
    std::array<SpecificEnergy<T, U>, M> out {};
    for (std::size_t m = 0; m < M; ++m) {
        out[m] = SpecificEnergy<T, U> {T(rho_e_m[m].value() / rho.value())};
    }
    return out;
}

template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE T energy_sum(const std::array<SpecificEnergy<T, U>, M>& mode_e)
{
    T sum = mode_e[0].value();
    for (std::size_t m = 1; m < M; ++m) { sum = sum + mode_e[m].value(); }
    return sum;
}

// Closure of the state with the given members. The status combines the
// status of (rho, e, T, p), of (rho, e_m, T_m) for every mode (see
// thermally_perfect_status) and status.
template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  multi_temperature_closure(Density<T, U> rho,
                            SpecificEnergy<T, U> e,
                            Pressure<T, U> p,
                            Temperature<T, U> temp,
                            const std::array<SpecificEnergy<T, U>, M>& mode_e,
                            const std::array<ModeTemperature<T, U>, M>& mode_temp,
                            Status status)
{
    Status s = combine(thermally_perfect_status(rho, e, temp, p), status);
    for (std::size_t m = 0; m < M; ++m) {
        s = combine(
          s,
          thermally_perfect_status(rho, mode_e[m], Temperature<T, U> {mode_temp[m].value()}));
    }
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = p,
          .temp = temp,
          .mode_e = mode_e,
          .mode_temp = mode_temp,
        },
      .status = s,
    };
}

// p = rho (R_h T + y_e R_e T_e).
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Pressure<T, U>
  multi_temperature_pressure(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                             const HeavySpeciesSums<T>& sums,
                             Density<T, U> rho,
                             Temperature<T, U> temp,
                             const std::array<ModeTemperature<T, U>, M>& mode_temp)
{
    const T t_e = electron_temperature(gas, mode_temp, temp.value());
    return Pressure<T, U> {
      T(rho.value() * ((sums.gas_constant * temp.value()) + (sums.electron_gas_constant * t_e))),
    };
}

// State from rho, T and the mode temperatures: e_m(T_m), e = e_tr(T) +
// sum(e_m), p = rho (R_h T + y_e R_e T_e). The status combines
// multi_temperature_closure with status.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  multi_temperature_from_temperatures(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                      const NormalisedComposition<T, U, N>& mix,
                                      const HeavySpeciesSums<T>& sums,
                                      Density<T, U> rho,
                                      Temperature<T, U> temp,
                                      const std::array<ModeTemperature<T, U>, M>& mode_temp,
                                      Status status,
                                      Pressure<T, U> p)
{
    std::array<SpecificEnergy<T, U>, M> mode_e {};
    for (std::size_t m = 0; m < M; ++m) {
        mode_e[m] =
          SpecificEnergy<T, U> {mode_energy(gas, mix.fractions, m, mode_temp[m].value()).e};
    }
    const SpecificEnergy<T, U> e {
      T(translational_energy(gas, sums, temp.value()) + energy_sum(mode_e)),
    };
    return multi_temperature_closure(rho,
                                     e,
                                     p,
                                     temp,
                                     mode_e,
                                     mode_temp,
                                     combine(status, mix.status));
}

template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  multi_temperature_from_temperatures(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                      const NormalisedComposition<T, U, N>& mix,
                                      const HeavySpeciesSums<T>& sums,
                                      Density<T, U> rho,
                                      Temperature<T, U> temp,
                                      const std::array<ModeTemperature<T, U>, M>& mode_temp,
                                      Status status)
{
    return multi_temperature_from_temperatures(
      gas,
      mix,
      sums,
      rho,
      temp,
      mode_temp,
      status,
      multi_temperature_pressure(gas, sums, rho, temp, mode_temp));
}

// Every mode temperature at temp.
template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE std::array<ModeTemperature<T, U>, M>
  equal_mode_temperatures(Temperature<T, U> temp)
{
    std::array<ModeTemperature<T, U>, M> out {};
    for (std::size_t m = 0; m < M; ++m) { out[m] = ModeTemperature<T, U> {temp.value()}; }
    return out;
}

// T_m of mode m from e_m by solve_temperature of the single-temperature
// model on e_m(T_m) - e_m = 0, with the derivative cv_m(T_m), starting from
// guess. y holds the value parts of the mass fractions.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE SolveResult<T> mode_temperature_from_energy(
  const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
  const NormalisedComposition<T, U, N>& mix,
  const std::array<MassFraction<typename ScalarTraits<T>::value_type>, N>& y,
  std::size_t m,
  SpecificEnergy<T, U> e_m,
  typename ScalarTraits<T>::value_type guess)
{
    using V = ScalarTraits<T>::value_type;
    const V target = ScalarTraits<T>::value(e_m.value());
    const auto f_value = [&gas, &y, m, target](V t) {
        const EnergyAndCv<V> r = mode_energy(gas, y, m, t);
        return Residual<V> {.value = r.e - target, .derivative = r.cv};
    };
    const auto f_full = [&gas, &mix, m, e_m](T t) {
        const EnergyAndCv<T> r = mode_energy(gas, mix.fractions, m, t);
        return Residual<T> {.value = T(r.e - e_m.value()), .derivative = r.cv};
    };
    return solve_temperature<T>(gas.thermally_perfect, f_value, f_full, guess);
}

// Mode temperatures from mode energies: the temperatures, which modes have
// species present, and the combined solve status.
template <Scalar T, UnitSystem U, std::size_t M> struct ModeSolve {
    std::array<ModeTemperature<T, U>, M> temp;
    std::array<bool, M> present;
    Status status;
};

// Solves every mode except skip (pass M to solve all) with
// mode_temperature_from_energy from guess. A mode without species present
// is not solved: its temperature is guess and present[m] is false; the
// caller sets it to T with fill_absent_modes.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE void solve_mode_temperatures(
  const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
  const NormalisedComposition<T, U, N>& mix,
  const std::array<MassFraction<typename ScalarTraits<T>::value_type>, N>& y,
  const std::array<SpecificEnergy<T, U>, M>& mode_e,
  typename ScalarTraits<T>::value_type guess,
  std::size_t skip,
  ModeSolve<T, U, M>& out)
{
    for (std::size_t m = 0; m < M; ++m) {
        if (m == skip) { continue; }
        out.present[m] = mode_present(gas, y, m);
        if (!out.present[m]) {
            out.temp[m] = ModeTemperature<T, U> {T(guess)};
            continue;
        }
        const SolveResult<T> r = mode_temperature_from_energy(gas, mix, y, m, mode_e[m], guess);
        out.temp[m] = ModeTemperature<T, U> {r.root};
        out.status = combine(out.status, r.status);
    }
}

// Sets the temperature of every mode without species present to temp.
template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE void fill_absent_modes(ModeSolve<T, U, M>& solved, Temperature<T, U> temp)
{
    for (std::size_t m = 0; m < M; ++m) {
        if (!solved.present[m]) { solved.temp[m] = ModeTemperature<T, U> {temp.value()}; }
    }
}

// A ModeSolve with every mode present and status Valid.
template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE ModeSolve<T, U, M> empty_mode_solve()
{
    ModeSolve<T, U, M> out {};
    out.present.fill(true);
    out.status = Status::Valid;
    return out;
}

// T_e of the mode that holds e- at given p and T, from that mode's energy
// density rho_e_m. With rho = p / (R_h T + y_e R_e T_e), it is the root of
// e_m(T_e) - rho_e_m (R_h T + y_e R_e T_e) / p, with the derivative
// cv_m(T_e) - rho_e_m y_e R_e / p, by solve_temperature starting from T.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE SolveResult<T> electron_temperature_at_pressure(
  const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
  const NormalisedComposition<T, U, N>& mix,
  const std::array<MassFraction<typename ScalarTraits<T>::value_type>, N>& y,
  const HeavySpeciesSums<T>& sums,
  Pressure<T, U> p,
  Temperature<T, U> temp,
  ModeEnergyDensity<T, U> rho_e_m)
{
    using V = ScalarTraits<T>::value_type;
    const auto m = static_cast<std::size_t>(gas.electron_mode);
    const V k = ScalarTraits<T>::value(rho_e_m.value()) / ScalarTraits<T>::value(p.value());
    const V heavy =
      ScalarTraits<T>::value(sums.gas_constant) * ScalarTraits<T>::value(temp.value());
    const V electron = ScalarTraits<T>::value(sums.electron_gas_constant);
    const auto f_value = [&gas, &y, m, k, heavy, electron](V t) {
        const EnergyAndCv<V> r = mode_energy(gas, y, m, t);
        return Residual<V> {
          .value = r.e - (k * (heavy + (electron * t))),
          .derivative = r.cv - (k * electron),
        };
    };
    const auto f_full = [&gas, &mix, &sums, m, p, temp, rho_e_m](T t) {
        const EnergyAndCv<T> r = mode_energy(gas, mix.fractions, m, t);
        const T k_full = T(rho_e_m.value() / p.value());
        const T heavy_full = T(sums.gas_constant * temp.value());
        return Residual<T> {
          .value = T(r.e - (k_full * T(heavy_full + (sums.electron_gas_constant * t)))),
          .derivative = T(r.cv - (k_full * sums.electron_gas_constant)),
        };
    };
    return solve_temperature<T>(gas.thermally_perfect,
                                f_value,
                                f_full,
                                ScalarTraits<T>::value(temp.value()));
}

}  // namespace detail

namespace detail
{

// Name of a component in error messages.
inline const char* mode_component_name(ModeComponent c)
{
    switch (c) {
    case ModeComponent::vibration: return "vibration";
    case ModeComponent::electronic: return "electronic";
    case ModeComponent::electron: return "electron";
    }
    return "unknown";
}

// Species checks of make_multi_temperature_mixture; throws ParameterError.
template <std::size_t N, std::size_t S>
void check_multi_temperature_species(const std::array<Nasa9Species<S>, N>& species,
                                     const std::array<std::string_view, N>& names,
                                     const std::array<int, N>& rotational_degrees)
{
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (species_index(names, names[i.value()]) != OptionalSpeciesIndex {i}) {
            throw ParameterError(std::format("species {}: name given twice", names[i.value()]));
        }
        const int d = rotational_degrees[i.value()];
        if (species[i.value()].metadata.role == SpeciesRole::electron) {
            if (d != 0) {
                throw ParameterError(
                  std::format("species {}: rotational degrees of freedom must be 0, got {}",
                              names[i.value()],
                              d));
            }
        }
        if (d != 0 && d != 2 && d != 3) {
            throw ParameterError(
              std::format("species {}: rotational degrees of freedom must be 0, 2 or 3, got {}",
                          names[i.value()],
                          d));
        }
    }
}

// Levels of one species in model units; throws ParameterError.
template <std::floating_point P, UnitSystem U, std::size_t L>
ElectronicLevels<P, U, L> make_electronic_levels(const ElectronicLevelsInput<L>& input,
                                                 const UnitScales<P, U>& scales)
{
    if (input.levels < 0 || static_cast<std::size_t>(input.levels) > L) {
        throw ParameterError(
          std::format("level count must be between 0 and {}, got {}", L, input.levels));
    }
    const auto tau = scale_factor<kind::temperature>(scales);
    ElectronicLevels<P, U, L> out {};
    out.levels = input.levels;
    const auto n = static_cast<std::size_t>(input.levels);
    for (std::size_t k = 0; k < n; ++k) {
        const double theta = input.temperatures[k];
        const int g = input.degeneracies[k];
        if (!std::isfinite(theta) || theta < 0.0) {
            throw ParameterError(
              std::format("level {}: temperature must be finite and non-negative, got {}",
                          k,
                          theta));
        }
        if (k == 0 && theta != 0.0) {
            throw ParameterError(
              std::format("level 0: the ground level temperature must be 0, got {}", theta));
        }
        if (g <= 0) {
            throw ParameterError(
              std::format("level {}: degeneracy must be positive, got {}", k, g));
        }
        const long double converted =
          static_cast<long double>(theta) / static_cast<long double>(tau);
        if (!std::isfinite(converted) || converted > std::numeric_limits<P>::max()) {
            throw ParameterError(
              std::format("level {}: converted temperature is outside the finite parameter range",
                          k));
        }
        const Temperature<P, U> stored {static_cast<P>(converted)};
        if (theta > 0.0 && !(stored.value() > P(0))) {
            throw ParameterError(
              std::format("level {}: positive temperature rounds to zero in model units", k));
        }
        out.degeneracies[k] = static_cast<P>(g);
        out.temperatures[k] = stored;
    }
    return out;
}

// Mode components of make_multi_temperature_mixture; throws ParameterError.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          std::size_t C>
void fill_mode_components(MultiTemperatureMixture<P, U, N, M, S, L>& model,
                          const std::array<std::string_view, N>& names,
                          const std::array<Nasa9Species<S>, N>& species,
                          const std::array<int, N>& rotational_degrees,
                          const std::array<EnergyModeInput<C>, M>& modes)
{
    std::array<std::array<int, 3>, N> owner {};
    for (auto& o : owner) { o.fill(-1); }
    for (std::size_t m = 0; m < M; ++m) {
        const int count = modes[m].count;
        if (count < 1 || static_cast<std::size_t>(count) > C) {
            throw ParameterError(
              std::format("mode {}: the component count must be between 1 and {}, got {}",
                          m,
                          C,
                          count));
        }
        const auto n = static_cast<std::size_t>(count);
        for (std::size_t k = 0; k < n; ++k) {
            const ModeComponentInput& input = modes[m].components[k];
            const std::string_view name =
              input.species == nullptr ? std::string_view {} : std::string_view {input.species};
            const auto fail = [&](std::string_view reason) {
                return ParameterError(std::format("mode {}, component {} ({} {}): {}",
                                                  m,
                                                  k,
                                                  name,
                                                  mode_component_name(input.component),
                                                  reason));
            };
            const auto index = species_index(names, name);
            if (!index) { throw fail("unknown species"); }
            const SpeciesIndex species_id = *index;
            const auto i = species_id.value();
            const bool electron = species[i].metadata.role == SpeciesRole::electron;
            switch (input.component) {
            case ModeComponent::vibration:
                if (electron) { throw fail("free electrons have only the electron component"); }
                if (rotational_degrees[i] == 0) { throw fail("an atom has no vibration"); }
                break;
            case ModeComponent::electronic:
                if (electron) { throw fail("free electrons have only the electron component"); }
                break;
            case ModeComponent::electron:
                if (!electron) {
                    throw fail("only an electron-role species has the electron component");
                }
                break;
            default: throw fail("unknown component");
            }
            const auto c = static_cast<std::size_t>(input.component);
            if (owner[i][c] >= 0) { throw fail(std::format("already in mode {}", owner[i][c])); }
            owner[i][c] = static_cast<int>(m);
            ModeComponentSet& set = model.mode_components[m][i];
            set.vibration = set.vibration || input.component == ModeComponent::vibration;
            set.electronic = set.electronic || input.component == ModeComponent::electronic;
            set.electron = set.electron || input.component == ModeComponent::electron;
        }
    }
    model.electron = {};
    model.electron_mode = -1;
    for (const SpeciesIndex e : species_indices(0, N)) {
        if (species[e.value()].metadata.role != SpeciesRole::electron) { continue; }
        const int mode = owner[e.value()][static_cast<std::size_t>(ModeComponent::electron)];
        if (mode < 0) {
            throw ParameterError(std::format(
              "species {}: free electrons must be in a mode as their electron component",
              names[e.value()]));
        }
        model.electron = OptionalSpeciesIndex {e};
        model.electron_mode = mode;
    }
}

}  // namespace detail

// Multi-temperature model of the given species, with d_i in
// rotational_degrees, electronic levels looked up by species name in levels
// (a species not in levels has none) and the energy modes of modes, in the
// unit system of the scales, with its species data. The single-temperature
// model and the species data are those of make_thermally_perfect_mixture,
// with its temperature limits and solve settings. Level temperatures are
// converted to the temperature unit of U in long double precision and then
// rounded to P. Throws ParameterError, naming the species, the levels entry
// or the mode and component, when:
// - a species name is null or empty (the error gives the species index);
// - two species have the same name; an electron-role species does not have
//   charge -1 or has rotational degrees of freedom other than 0; a number
//   of rotational degrees of freedom is not 0, 2 or 3;
// - make_thermally_perfect_mixture rejects a species;
// - a levels entry names an unknown species or a species already given; its
//   level count is not in [0, L]; a level temperature is not finite and
//   non-negative; the first level temperature is not 0; a degeneracy is not
//   positive; a converted temperature is nonfinite or a positive input rounds
//   to zero; a stored reference electronic energy is nonfinite or negative;
// - a stored translational/rotational cp is nonfinite or nonpositive, or
//   a stored reference formation enthalpy is nonfinite;
// - a mode has a component count not in [1, C]; a component names an
//   unknown species; it is vibration of an atom (d_i = 0) or a free electron,
//   electronic of a free electron, electron of a heavy species, or not one of
//   the three components; a component appears twice (in one mode or in two);
// - an electron-role species is present but its electron component is in no mode.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          std::size_t L,
          std::size_t K,
          std::size_t C,
          std::size_t M>
    requires(N > 0 && S > 0 && L > 0 && C > 0 && M > 0)
Built<MultiTemperatureMixture<P, U, N, M, S, L>, SpeciesData<P, U, N>>
  make_multi_temperature_mixture(const std::array<Nasa9Species<S>, N>& species,
                                 const std::array<int, N>& rotational_degrees,
                                 const std::array<ElectronicLevelsInput<L>, K>& levels,
                                 const std::array<EnergyModeInput<C>, M>& modes,
                                 const UnitScales<P, U>& scales)
{
    std::array<std::string_view, N> names {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (species[i.value()].name == nullptr || species[i.value()].name[0] == '\0') {
            throw ParameterError(
              std::format("species {}: name must not be null or empty", i.value()));
        }
        names[i.value()] = species[i.value()].name;
    }
    detail::check_multi_temperature_species(species, names, rotational_degrees);
    const auto single = make_thermally_perfect_mixture(species, scales);

    Built<MultiTemperatureMixture<P, U, N, M, S, L>, SpeciesData<P, U, N>> built {
      .model = {.thermally_perfect = single.model,
                .translational_rotational_cp = {},
                .formation_enthalpy = {},
                .electronic_levels = {},
                .reference_electronic_energy = {},
                .mode_components = {},
                .electron = {},
                .electron_mode = {}},
      .species = single.species};
    auto& model = built.model;
    const P t_ref = single.model.reference_temperature.value();
    const P ln_t_ref = std::log(t_ref);
    for (const SpeciesIndex i : species_indices(0, N)) {
        const P r = single.model.gas_constant[i.value()].value();
        const auto d = static_cast<P>(rotational_degrees[i.value()]);
        try {
            model.translational_rotational_cp[i.value()] =
              detail::positive_thermodynamic_parameter<P>(
                SpecificHeat<P, U> {((P(5) + d) / P(2)) * r},
                "translational/rotational specific heat");
            model.formation_enthalpy[i.value()] = detail::thermodynamic_parameter<P>(
              SpecificEnthalpy<P, U> {
                r * nasa9_cp_h(single.model.curve[i.value()], t_ref, ln_t_ref).h_over_r},
              "reference formation enthalpy");
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", names[i.value()], error.what()));
        }
    }

    std::array<bool, N> has_levels {};
    for (std::size_t k = 0; k < K; ++k) {
        const ElectronicLevelsInput<L>& input = levels[k];
        const std::string_view name =
          input.name == nullptr ? std::string_view {} : std::string_view {input.name};
        const auto index = detail::species_index(names, name);
        if (!index) {
            throw ParameterError(std::format("electronic levels of {}: unknown species", name));
        }
        const SpeciesIndex species_id = *index;
        const auto i = species_id.value();
        if (has_levels[i]) {
            throw ParameterError(std::format("electronic levels of {}: given twice", name));
        }
        has_levels[i] = true;
        try {
            model.electronic_levels[i] = detail::make_electronic_levels(input, scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("electronic levels of {}: {}", name, error.what()));
        }
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        model.reference_electronic_energy[i.value()] = SpecificEnergy<P, U> {
          detail::electronic_energy(model.electronic_levels[i.value()],
                                    single.model.gas_constant[i.value()],
                                    t_ref)
            .e,
        };
        const P reference_energy = model.reference_electronic_energy[i.value()].value();
        if (!std::isfinite(reference_energy) || reference_energy < P(0)) {
            throw ParameterError(std::format("electronic levels of {}: reference electronic energy "
                                             "must be finite and non-negative, got {}",
                                             names[i.value()],
                                             reference_energy));
        }
    }

    detail::fill_mode_components(model, names, species, rotational_degrees, modes);
    return built;
}

// Per-cell composition from mass fractions, as for the single-temperature
// model.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P> &&
             ScalarWith<typename range_quantity_t<V>::value_type, P>
CANARY_HOST_DEVICE auto mixture(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                const V& fractions)
{ return mixture(gas.thermally_perfect, fractions); }

// Per-cell composition from partial densities and the density, as for the
// single-temperature model.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V,
          Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T> && ScalarWith<T, P>
CANARY_HOST_DEVICE auto mixture(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                const V& partial_densities,
                                Density<T, U> rho)
{ return mixture(gas.thermally_perfect, partial_densities, rho); }

// Closures at thermal equilibrium (no mode input): every T_m = T, and the
// state is that of the closure from (rho, T, T_m = T) with the given inputs
// kept. The status combines detail::multi_temperature_closure, the solve
// status where T is solved for, and the composition status.

// T by solve_temperature of the single-temperature model on e(T) - e = 0,
// with every mode at T and the derivative cv(T), from the single-temperature
// model's default initial guess; e is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        SpecificEnergy<T, U> e)
{
    using V = ScalarTraits<T>::value_type;
    const std::array<MassFraction<V>, N> y = detail::fraction_values(mix);
    const detail::HeavySpeciesSums<V> sums_v = detail::heavy_species_sums(gas, y);
    const detail::HeavySpeciesSums<T> sums = detail::heavy_species_sums(gas, mix.fractions);
    const V target = ScalarTraits<T>::value(e.value());
    const auto f_value = [&gas, &y, &sums_v, target](V t) {
        const detail::EnergyAndCv<V> r = detail::equilibrium_energy(gas, y, sums_v, t);
        return detail::Residual<V> {.value = r.e - target, .derivative = r.cv};
    };
    const auto f_full = [&gas, &mix, &sums, e](T t) {
        const detail::EnergyAndCv<T> r = detail::equilibrium_energy(gas, mix.fractions, sums, t);
        return detail::Residual<T> {.value = T(r.e - e.value()), .derivative = r.cv};
    };
    const V guess = detail::default_temperature_guess(gas.thermally_perfect, y, target);
    const detail::SolveResult<T> r =
      detail::solve_temperature<T>(gas.thermally_perfect, f_value, f_full, guess);
    const Temperature<T, U> temp {r.root};
    auto c =
      detail::multi_temperature_from_temperatures(gas,
                                                  mix,
                                                  sums,
                                                  rho,
                                                  temp,
                                                  detail::equal_mode_temperatures<T, U, M>(temp),
                                                  r.status);
    c.state.e = e;
    return c;
}

// T = p / (rho R), with R = sum(y_i R_i) over all species; p is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Pressure<T, U> p)
{
    const Temperature<T, U> temp {temperature(p, rho, mix.gas_constant)};
    auto c =
      detail::multi_temperature_from_temperatures(gas,
                                                  mix,
                                                  detail::heavy_species_sums(gas, mix.fractions),
                                                  rho,
                                                  temp,
                                                  detail::equal_mode_temperatures<T, U, M>(temp),
                                                  Status::Valid,
                                                  p);
    return c;
}

template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Temperature<T, U> temp)
{
    return detail::multi_temperature_from_temperatures(
      gas,
      mix,
      detail::heavy_species_sums(gas, mix.fractions),
      rho,
      temp,
      detail::equal_mode_temperatures<T, U, M>(temp),
      Status::Valid);
}

// rho = p / (R T), with R = sum(y_i R_i) over all species; p is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Pressure<T, U> p,
        Temperature<T, U> temp)
{
    const Density<T, U> rho {density(p, mix.gas_constant, temp)};
    auto c =
      detail::multi_temperature_from_temperatures(gas,
                                                  mix,
                                                  detail::heavy_species_sums(gas, mix.fractions),
                                                  rho,
                                                  temp,
                                                  detail::equal_mode_temperatures<T, U, M>(temp),
                                                  Status::Valid,
                                                  p);
    return c;
}

// Closures with the mode temperatures given as a mode temperature range of
// M elements. e_m = e_m(T_m). The status combines
// detail::multi_temperature_closure and the composition status.

// T in closed form from e_tr = e - sum(e_m)
// (detail::translational_temperature); e is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_temperature, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        SpecificEnergy<T, U> e,
        const V& mode_temperatures)
{
    const auto mode_temp = detail::mode_values<M>(mode_temperatures);
    std::array<SpecificEnergy<T, U>, M> mode_e {};
    for (std::size_t m = 0; m < M; ++m) {
        mode_e[m] =
          SpecificEnergy<T, U> {detail::mode_energy(gas, mix.fractions, m, mode_temp[m].value()).e};
    }
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    const Temperature<T, U> temp {
      detail::translational_temperature(gas, sums, T(e.value() - detail::energy_sum(mode_e))),
    };
    return detail::multi_temperature_closure(
      rho,
      e,
      detail::multi_temperature_pressure(gas, sums, rho, temp, mode_temp),
      temp,
      mode_e,
      mode_temp,
      mix.status);
}

// T = (p / rho - y_e R_e T_e) / R_h; p is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_temperature, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Pressure<T, U> p,
        const V& mode_temperatures)
{
    const auto mode_temp = detail::mode_values<M>(mode_temperatures);
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    const T t_e = detail::electron_temperature(gas, mode_temp, T(p.value() / rho.value()));
    const Temperature<T, U> temp {
      T((T(p.value() / rho.value()) - T(sums.electron_gas_constant * t_e)) / sums.gas_constant),
    };
    auto c = detail::multi_temperature_from_temperatures(gas,
                                                         mix,
                                                         sums,
                                                         rho,
                                                         temp,
                                                         mode_temp,
                                                         Status::Valid,
                                                         p);
    return c;
}

template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_temperature, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Temperature<T, U> temp,
        const V& mode_temperatures)
{
    return detail::multi_temperature_from_temperatures(
      gas,
      mix,
      detail::heavy_species_sums(gas, mix.fractions),
      rho,
      temp,
      detail::mode_values<M>(mode_temperatures),
      Status::Valid);
}

// rho = p / (R_h T + y_e R_e T_e); p is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_temperature, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Pressure<T, U> p,
        Temperature<T, U> temp,
        const V& mode_temperatures)
{
    const auto mode_temp = detail::mode_values<M>(mode_temperatures);
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    const T t_e = detail::electron_temperature(gas, mode_temp, temp.value());
    const Density<T, U> rho {
      T(p.value() / T((sums.gas_constant * temp.value()) + (sums.electron_gas_constant * t_e))),
    };
    auto c = detail::multi_temperature_from_temperatures(gas,
                                                         mix,
                                                         sums,
                                                         rho,
                                                         temp,
                                                         mode_temp,
                                                         Status::Valid,
                                                         p);
    return c;
}

// Closures with the mode energies given as a mode energy density range of M
// elements, rho e_m; e_m = (rho e_m) / rho. Each T_m is found by
// detail::mode_temperature_from_energy, starting from T. A mode whose
// species are all absent (sum of their mass fractions 0) is not solved; its
// temperature is T. The status combines detail::multi_temperature_closure,
// the solve statuses (NotConverged, or Invalid when no bracket has a sign
// change) and the composition status.

// T in closed form from e_tr = e - sum(e_m)
// (detail::translational_temperature); e is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_energy_density, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        SpecificEnergy<T, U> e,
        const V& mode_energy_densities)
{
    const auto mode_e =
      detail::mode_specific_energies(detail::mode_values<M>(mode_energy_densities), rho);
    const auto y = detail::fraction_values(mix);
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    const Temperature<T, U> temp {
      detail::translational_temperature(gas, sums, T(e.value() - detail::energy_sum(mode_e))),
    };
    auto solved = detail::empty_mode_solve<T, U, M>();
    detail::solve_mode_temperatures(gas,
                                    mix,
                                    y,
                                    mode_e,
                                    ScalarTraits<T>::value(temp.value()),
                                    M,
                                    solved);
    detail::fill_absent_modes(solved, temp);
    return detail::multi_temperature_closure(
      rho,
      e,
      detail::multi_temperature_pressure(gas, sums, rho, temp, solved.temp),
      temp,
      mode_e,
      solved.temp,
      combine(solved.status, mix.status));
}

// The mode temperatures are solved from p / (rho R), with R = sum(y_i R_i)
// over all species; then T = (p / rho - y_e R_e T_e) / R_h. p is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_energy_density, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Pressure<T, U> p,
        const V& mode_energy_densities)
{
    const auto mode_e =
      detail::mode_specific_energies(detail::mode_values<M>(mode_energy_densities), rho);
    const auto y = detail::fraction_values(mix);
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    const Temperature<T, U> equilibrium {temperature(p, rho, mix.gas_constant)};
    auto solved = detail::empty_mode_solve<T, U, M>();
    detail::solve_mode_temperatures(gas,
                                    mix,
                                    y,
                                    mode_e,
                                    ScalarTraits<T>::value(equilibrium.value()),
                                    M,
                                    solved);
    const T p_over_rho = T(p.value() / rho.value());
    const T t_e = detail::electron_temperature(gas, solved.temp, p_over_rho);
    const Temperature<T, U> temp {
      T((p_over_rho - T(sums.electron_gas_constant * t_e)) / sums.gas_constant),
    };
    detail::fill_absent_modes(solved, temp);
    const SpecificEnergy<T, U> e {
      T(detail::translational_energy(gas, sums, temp.value()) + detail::energy_sum(mode_e)),
    };
    return detail::multi_temperature_closure(rho,
                                             e,
                                             p,
                                             temp,
                                             mode_e,
                                             solved.temp,
                                             combine(solved.status, mix.status));
}

// e = e_tr(T) + sum(e_m), p = rho (R_h T + y_e R_e T_e).
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_energy_density, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Temperature<T, U> temp,
        const V& mode_energy_densities)
{
    const auto mode_e =
      detail::mode_specific_energies(detail::mode_values<M>(mode_energy_densities), rho);
    const auto y = detail::fraction_values(mix);
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    auto solved = detail::empty_mode_solve<T, U, M>();
    detail::solve_mode_temperatures(gas,
                                    mix,
                                    y,
                                    mode_e,
                                    ScalarTraits<T>::value(temp.value()),
                                    M,
                                    solved);
    detail::fill_absent_modes(solved, temp);
    const SpecificEnergy<T, U> e {
      T(detail::translational_energy(gas, sums, temp.value()) + detail::energy_sum(mode_e)),
    };
    return detail::multi_temperature_closure(
      rho,
      e,
      detail::multi_temperature_pressure(gas, sums, rho, temp, solved.temp),
      temp,
      mode_e,
      solved.temp,
      combine(solved.status, mix.status));
}

// rho = p / (R_h T + y_e R_e T_e), then e_m = (rho e_m) / rho. With e-
// present, T_e and rho are solved together by
// detail::electron_temperature_at_pressure, starting from T; the other mode
// temperatures then as above. p is kept.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> &&
             ModeRangeOf<V, kind::mode_energy_density, T, U>
CANARY_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, M>>
  close(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Pressure<T, U> p,
        Temperature<T, U> temp,
        const V& mode_energy_densities)
{
    using Value = ScalarTraits<T>::value_type;
    const auto rho_e_m = detail::mode_values<M>(mode_energy_densities);
    const auto y = detail::fraction_values(mix);
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    auto solved = detail::empty_mode_solve<T, U, M>();
    std::size_t skip = M;
    T t_e = temp.value();
    if (gas.electron.has_value() && y[(*gas.electron).value()].value() > Value(0)) {
        skip = static_cast<std::size_t>(gas.electron_mode);
        const detail::SolveResult<T> r =
          detail::electron_temperature_at_pressure(gas, mix, y, sums, p, temp, rho_e_m[skip]);
        t_e = r.root;
        solved.temp[skip] = ModeTemperature<T, U> {r.root};
        solved.status = r.status;
    }
    const Density<T, U> rho {
      T(p.value() / T((sums.gas_constant * temp.value()) + (sums.electron_gas_constant * t_e))),
    };
    const auto mode_e = detail::mode_specific_energies(rho_e_m, rho);
    detail::solve_mode_temperatures(gas,
                                    mix,
                                    y,
                                    mode_e,
                                    ScalarTraits<T>::value(temp.value()),
                                    skip,
                                    solved);
    detail::fill_absent_modes(solved, temp);
    const SpecificEnergy<T, U> e {
      T(detail::translational_energy(gas, sums, temp.value()) + detail::energy_sum(mode_e)),
    };
    return detail::multi_temperature_closure(rho,
                                             e,
                                             p,
                                             temp,
                                             mode_e,
                                             solved.temp,
                                             combine(solved.status, mix.status));
}

// cv = sum over heavy i of y_i cv_tr,i + sum_m cv_m(T_m), cp = cv + R_h,
// gamma = cp / cv, the frozen sound speed a = sqrt(gamma p / rho), and the
// mixture gas constant R = sum(y_i R_i) over all species.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE Properties<T, U> properties(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                               const NormalisedComposition<T, U, N>& mix,
                                               const MultiTemperatureState<T, U, M>& state)
{
    const auto sums = detail::heavy_species_sums(gas, mix.fractions);
    T cv_value = T(sums.cp - sums.gas_constant);
    for (std::size_t m = 0; m < M; ++m) {
        cv_value =
          cv_value + detail::mode_energy(gas, mix.fractions, m, state.mode_temp[m].value()).cv;
    }
    const SpecificHeat<T, U> cv {cv_value};
    const SpecificHeat<T, U> cp {T(cv_value + sums.gas_constant)};
    const HeatCapacityRatio<T> gamma {heat_capacity_ratio(cp, cv)};
    return {
      .cv = cv,
      .cp = cp,
      .gamma = gamma,
      .sound_speed = SoundSpeed<T, U> {sound_speed(gamma, state.p, state.rho)},
      .gas_constant = mix.gas_constant,
    };
}

// Electron pressure p_e = rho y_e R_e T_e; 0 without e-.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Pressure<T, U>
  electron_pressure(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    const MultiTemperatureState<T, U, M>& state)
{
    using V = ScalarTraits<T>::value_type;
    if (!gas.electron) { return Pressure<T, U> {T(V(0))}; }
    const auto e = (*gas.electron).value();
    const auto m = static_cast<std::size_t>(gas.electron_mode);
    const T y_r = T(mix.fractions[e].value() * gas.thermally_perfect.gas_constant[e].value());
    return Pressure<T, U> {T(state.rho.value() * y_r * state.mode_temp[m].value())};
}

// Per-species properties of the single-temperature model at temperature
// temp (thermal equilibrium); see SpeciesProperties.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE SpeciesProperties<T, U, N>
  species_properties(const MultiTemperatureMixture<P, U, N, M, S, L>& gas, Temperature<T, U> temp)
{ return species_properties(gas.thermally_perfect, temp); }

// Isobaric specific heat cp_i of every species of the single-temperature
// model at temperature temp (thermal equilibrium).
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE std::array<SpecificHeat<T, U>, N>
  species_isobaric_specific_heat(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                                 Temperature<T, U> temp)
{ return species_isobaric_specific_heat(gas.thermally_perfect, temp); }

// The mode that holds the vibration of species i, or -1 when its vibration
// is in no mode (an atom, a free electron, or a frozen vibration).
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
CANARY_HOST_DEVICE constexpr int
  species_vibration_mode(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                         SpeciesIndex index) noexcept
{
    const auto i = index.value();
    if (i >= N) { return -1; }
    for (std::size_t m = 0; m < M; ++m) {
        if (gas.mode_components[m][i].vibration) { return static_cast<int>(m); }
    }
    return -1;
}

// The free-electron species index resolved from host metadata, or absence.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
CANARY_HOST_DEVICE constexpr OptionalSpeciesIndex
  electron_species_index(const MultiTemperatureMixture<P, U, N, M, S, L>& gas) noexcept
{ return gas.electron; }

// The mode that holds the energy of the free electrons, or -1 without free electrons.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
CANARY_HOST_DEVICE constexpr int
  electron_energy_mode(const MultiTemperatureMixture<P, U, N, M, S, L>& gas) noexcept
{ return gas.electron_mode; }

// Energy per unit mass e_i,m(T') of the components of every species i that
// mode m holds, at temperature temp (model units): the vibration, electronic
// or electron energy of the class comment, their sum when mode m holds
// both vibration and electronic energy of i, and 0 for a species with no
// component in mode m. One logarithm.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE std::array<SpecificEnergy<T, U>, N> species_mode_energy(
  const MultiTemperatureMixture<P, U, N, M, S, L>& gas, std::size_t m, Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    // NOLINTNEXTLINE(misc-static-assert): run-time mode index
    CANARY_ASSERT(m < M);
    const T t = temp.value();
    const T ln_t = math::log(t);
    std::array<SpecificEnergy<T, U>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        const ModeComponentSet c = gas.mode_components[m][i.value()];
        out[i.value()] = SpecificEnergy<T, U> {
          detail::holds_any(c) ? detail::component_energy(gas, i, c, t, ln_t).e : T(V(0)),
        };
    }
    return out;
}

// Detailed outcome of a requested charge-neutralization operation. None
// denotes success. NoNeutralisingSpecies means the electron convenience
// overload has no electron to select; a caller may select a negative ion.
enum class ChargeNeutralisationReason : unsigned char {
    None,
    NoNeutralisingSpecies,
    InvalidSpeciesCount,
    SpeciesIndexOutOfRange,
    NonFiniteFraction,
    NegativeFraction,
    NonFiniteMass,
    NonPositiveMass,
    NonFiniteCharge,
    NonNegativeNeutralisingCharge,
    NonFiniteChargePerMass,
    ChargePerMassUnderflow,
    NonFiniteWeightSum,
    NonFiniteChargeSum,
    NonFiniteMassPerCharge,
    MassPerChargeUnderflow,
    NonFiniteBalancingFraction,
    NegativeBalancingFraction,
    BalancingFractionUnderflow,
    NonFiniteNormalisationSum,
    NonPositiveNormalisationSum,
    NonFiniteNormalisedFraction,
};

// Fractions returned by charge neutralization and the calculation outcome.
// Invalid results retain the input fractions for diagnostics. This status
// describes the requested operation, not the admissibility of a charged fluid.
template <Scalar T, std::size_t N> struct ChargeNeutralFractions {
    std::array<MassFraction<T>, N> fractions;
    Status status;
    ChargeNeutralisationReason reason;
};

// Normalize a composition after adjusting the selected negative species b.
// Its weight becomes a = -(m_b/q_b) sum_{i != b}(q_i y_i/m_i); all weights
// are divided by a + sum_{i != b} y_i. Relative proportions of the other
// species are preserved. The selected species may be an electron or a
// negative ion. No species is selected implicitly by this overload.
//
// Requires fractions.extent(0) == N and non-throwing scalar operations and
// range access. Inputs need not sum to one and are not clamped. Invalid is
// returned for a mismatched species count, an out-of-range index, nonfinite
// or negative input weights, nonfinite charges, nonpositive/nonfinite masses,
// or a selected charge that is not negative. Nonfinite intermediate values,
// a negative required balancing weight, an underflowed nonzero charge ratio
// or balancing weight, and a nonpositive normalization sum also give Invalid.
// Validation uses scalar value parts; Dual derivatives are propagated but
// are not checked. Invalid retains the input fractions and reports the first
// failed check in reason. Checks follow species storage order.
// Valid returns normalized charge-neutral fractions, to floating-point
// accuracy. This is an optional composition operation, not a closure step.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P> &&
             ScalarWith<typename range_quantity_t<V>::value_type, P>
CANARY_HOST_DEVICE ChargeNeutralFractions<typename range_quantity_t<V>::value_type, N>
  charge_neutral_fractions(const MultiTemperatureMixture<P, U, N, M, S, L>& /*gas*/,
                           const SpeciesData<P, U, N>& species,
                           const V& fractions,
                           SpeciesIndex neutralising_species) noexcept
{
    using T = range_quantity_t<V>::value_type;
    using Value = ScalarTraits<T>::value_type;
    // NOLINTNEXTLINE(misc-static-assert): run-time extent
    CANARY_ASSERT(static_cast<std::size_t>(fractions.extent(0)) == N);
    ChargeNeutralFractions<T, N> out {.fractions = {},
                                      .status = Status::Invalid,
                                      .reason = ChargeNeutralisationReason::None};
    for (const SpeciesIndex i : species_indices(0, N)) {
        out.fractions[i.value()] = MassFraction<T> {fractions[i.value()].value()};
    }
    const auto b = neutralising_species.value();
    if (!std::cmp_equal(species.n, N)) {
        out.reason = ChargeNeutralisationReason::InvalidSpeciesCount;
        return out;
    }
    if (b >= N) {
        out.reason = ChargeNeutralisationReason::SpeciesIndexOutOfRange;
        return out;
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        const Value y = ScalarTraits<T>::value(out.fractions[i.value()].value());
        const P mass = species.mass[i.value()].value();
        if (!std::isfinite(y)) {
            out.reason = ChargeNeutralisationReason::NonFiniteFraction;
            return out;
        }
        if (y < Value(0)) {
            out.reason = ChargeNeutralisationReason::NegativeFraction;
            return out;
        }
        if (!std::isfinite(mass)) {
            out.reason = ChargeNeutralisationReason::NonFiniteMass;
            return out;
        }
        if (!(mass > P(0))) {
            out.reason = ChargeNeutralisationReason::NonPositiveMass;
            return out;
        }
        if (!std::isfinite(species.charge[i.value()].value())) {
            out.reason = ChargeNeutralisationReason::NonFiniteCharge;
            return out;
        }
    }
    if (!(species.charge[b].value() < P(0))) {
        out.reason = ChargeNeutralisationReason::NonNegativeNeutralisingCharge;
        return out;
    }
    T other_weight = T(Value(0));
    T charge = T(Value(0));
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (i.value() == b) { continue; }
        const T y = out.fractions[i.value()].value();
        const P q = species.charge[i.value()].value();
        const P charge_per_mass = q / species.mass[i.value()].value();
        if (!std::isfinite(charge_per_mass)) {
            out.reason = ChargeNeutralisationReason::NonFiniteChargePerMass;
            return out;
        }
        if (q != P(0) && charge_per_mass == P(0)) {
            out.reason = ChargeNeutralisationReason::ChargePerMassUnderflow;
            return out;
        }
        other_weight = other_weight + y;
        charge = charge + T(y * charge_per_mass);
    }
    if (!math::isfinite(other_weight)) {
        out.reason = ChargeNeutralisationReason::NonFiniteWeightSum;
        return out;
    }
    if (!math::isfinite(charge)) {
        out.reason = ChargeNeutralisationReason::NonFiniteChargeSum;
        return out;
    }
    const P mass_per_charge = species.mass[b].value() / species.charge[b].value();
    if (!std::isfinite(mass_per_charge)) {
        out.reason = ChargeNeutralisationReason::NonFiniteMassPerCharge;
        return out;
    }
    if (mass_per_charge == P(0)) {
        out.reason = ChargeNeutralisationReason::MassPerChargeUnderflow;
        return out;
    }
    const T balancing_weight = T(-(mass_per_charge * charge));
    const Value balancing_value = ScalarTraits<T>::value(balancing_weight);
    if (!std::isfinite(balancing_value)) {
        out.reason = ChargeNeutralisationReason::NonFiniteBalancingFraction;
        return out;
    }
    if (balancing_value < Value(0)) {
        out.reason = ChargeNeutralisationReason::NegativeBalancingFraction;
        return out;
    }
    if (ScalarTraits<T>::value(charge) != Value(0) && balancing_value == Value(0)) {
        out.reason = ChargeNeutralisationReason::BalancingFractionUnderflow;
        return out;
    }
    const T total = T(other_weight + balancing_weight);
    if (!math::isfinite(total)) {
        out.reason = ChargeNeutralisationReason::NonFiniteNormalisationSum;
        return out;
    }
    if (!(ScalarTraits<T>::value(total) > Value(0))) {
        out.reason = ChargeNeutralisationReason::NonPositiveNormalisationSum;
        return out;
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        const T weight = i.value() == b ? balancing_weight : out.fractions[i.value()].value();
        const T normalized = T(weight / total);
        if (!math::isfinite(normalized)) {
            out.reason = ChargeNeutralisationReason::NonFiniteNormalisedFraction;
            for (const SpeciesIndex j : species_indices(0, N)) {
                out.fractions[j.value()] = MassFraction<T> {fractions[j.value()].value()};
            }
            return out;
        }
        out.fractions[i.value()] = MassFraction<T> {normalized};
    }
    out.status = Status::Valid;
    return out;
}

// Charge neutralization using the model's electron species. The input extent
// must equal N, and scalar operations and range access must not throw.
// Without an electron, return the input unchanged with Invalid and
// NoNeutralisingSpecies; no negative ion is selected automatically, even if
// the input is already neutral. Otherwise the explicit-species overload
// supplies validation, normalization and diagnostics.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t S,
          std::size_t L,
          typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P> &&
             ScalarWith<typename range_quantity_t<V>::value_type, P>
CANARY_HOST_DEVICE ChargeNeutralFractions<typename range_quantity_t<V>::value_type, N>
  charge_neutral_fractions(const MultiTemperatureMixture<P, U, N, M, S, L>& gas,
                           const SpeciesData<P, U, N>& species,
                           const V& fractions) noexcept
{
    if (gas.electron) { return charge_neutral_fractions(gas, species, fractions, *gas.electron); }
    using T = range_quantity_t<V>::value_type;
    // NOLINTNEXTLINE(misc-static-assert): run-time extent
    CANARY_ASSERT(static_cast<std::size_t>(fractions.extent(0)) == N);
    ChargeNeutralFractions<T, N> out {.fractions = {},
                                      .status = Status::Invalid,
                                      .reason = ChargeNeutralisationReason::NoNeutralisingSpecies};
    for (const SpeciesIndex i : species_indices(0, N)) {
        out.fractions[i.value()] = MassFraction<T> {fractions[i.value()].value()};
    }
    return out;
}

namespace detail
{

template <std::size_t L> struct InputStorage<ElectronicLevelsInput<L>> {
    OwnedSpeciesName name {};
    int levels {};
    std::array<int, L> degeneracies {};
    std::array<double, L> temperatures {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const ElectronicLevelsInput<L>& input) :
        name(input.name), levels(input.levels), degeneracies(input.degeneracies),
        temperatures(input.temperatures)
    {
    }
    ElectronicLevelsInput<L> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .name = name.c_str(),
          .levels = levels,
          .degeneracies = degeneracies,
          .temperatures = temperatures,
        };
    }
};

}  // namespace detail

// Builds from owned and direct inputs, keeping borrowed views inside this call.
template <typename... Explicit, typename... Inputs>
    requires detail::HasOwnedInput<Inputs...>
auto make_multi_temperature_mixture(const Inputs&... inputs)
{
    return detail::with_inputs(
      [](const auto&... borrowed) {
          return make_multi_temperature_mixture<Explicit...>(borrowed...);
      },
      inputs...);
}

}  // namespace canary

#endif  // CANARY_MODELS_MULTI_TEMPERATURE_HPP
