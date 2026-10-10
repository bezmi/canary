#ifndef GASLIB_REACTIONS_HPP
#define GASLIB_REACTIONS_HPP

// Finite-rate chemistry of a mixture: a mechanism of reversible reactions
// with Arrhenius forward rates, backward rates from the equilibrium constant
// or from their own Arrhenius fits, and optional third bodies; per cell, the
// rates of progress of every reaction and the mass production rate of every
// species.
//
// Concentrations are number densities n_i = rho y_i / m_i (particles per
// unit volume), so that the same functions work in SI and in scaled units.
// For a reaction with reactant coefficients nu'_i and product coefficients
// nu''_i:
//   k_f = A T^n exp(-theta / T),
//   [M] = sum_i eps_i n_i with a third body, 1 without,
//   q_f = k_f [M] prod_i n_i^nu'_i,   q_b = k_b [M] prod_i n_i^nu''_i,
//   omega_i = m_i sum_r (nu''_ir - nu'_ir) (q_f,r - q_b,r).
// The backward rate constant is an Arrhenius fit of its own, or
// k_b = k_f / K_c with the equilibrium constant in concentration form
//   K_c = exp(-sum_i dnu_i g_i m_i / (k_B T)) (p0 / (k_B T))^(sum_i dnu_i),
// where dnu_i = nu''_i - nu'_i, g_i is the standard-state Gibbs energy per
// unit mass of species i and p0 the standard pressure. Every rate constant
// is evaluated at the translational temperature clamped to the mechanism's
// temperature limits.

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/config.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/species.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/units.hpp>
#include <limits>
#include <string_view>

namespace gaslib
{

// Host input of an Arrhenius rate constant k = A T^n exp(-theta / T), in SI
// per mole: A in (m^3/mol)^(order - 1) / s with T in K, where order is the
// sum of the coefficients of the reaction side the rate multiplies, plus 1
// with a third body; n dimensionless; activation_temperature theta in K.
struct ArrheniusInput {
    double A;
    double n;
    double activation_temperature;
};

// Host input of a species name and a number: a stoichiometric coefficient
// or a third-body efficiency.
struct SpeciesCoefficient {
    const char* name;
    double value;
};

// Host input of up to K species coefficients; only the first count entries
// are used.
template <std::size_t K>
    requires(K > 0)
struct SpeciesCoefficients {
    int count;
    std::array<SpeciesCoefficient, K> entries;
};

// How the backward rate constant of a reaction is found.
enum class BackwardRateModel : std::uint8_t {
    // k_b = k_f / K_c.
    equilibrium,
    // An Arrhenius rate constant of its own.
    arrhenius,
};

// Host input of one reversible reaction: its equation (for error messages),
// the reactant and product coefficients (positive integers), the forward
// rate, the backward rate model and, for BackwardRateModel::arrhenius, the
// backward rate, and whether a third body takes part, with the efficiencies
// of the listed species (species not listed have efficiency 1). With no
// third body, efficiencies must be empty.
template <std::size_t K>
    requires(K > 0)
struct ReactionInput {
    const char* equation;
    SpeciesCoefficients<K> reactants;
    SpeciesCoefficients<K> products;
    ArrheniusInput forward;
    BackwardRateModel backward_model;
    ArrheniusInput backward;
    bool third_body;
    SpeciesCoefficients<K> efficiencies;
};

// Host input of a mechanism of R reactions: the temperature limits of rate
// evaluation in K and the reactions, each with at most K species per list.
template <std::size_t R, std::size_t K>
    requires(R > 0 && K > 0)
struct MechanismInput {
    double temperature_lower;
    double temperature_upper;
    std::array<ReactionInput<K>, R> reactions;
};

// Arrhenius rate constant k = A T^n exp(-theta / T) in model units: T and
// theta in the temperature unit of U, and A per particle in
// (length^3 / count)^(order - 1) / time of U, times the temperature scale
// raised to n. A and n are plain values, since the dimension of A depends
// on the reaction order.
template <std::floating_point P, UnitSystem U> struct ArrheniusRate {
    P A;
    P n;
    Temperature<P, U> activation_temperature;
};

// One reaction of a mechanism over N species, indexed as the species of the
// thermo model: reactant and product coefficients, the sum of the
// coefficient changes sum_i (products[i] - reactants[i]), the forward rate,
// the backward rate model and rate (used with BackwardRateModel::arrhenius),
// whether a third body takes part and the efficiency of every species (used
// with a third body).
template <std::floating_point P, UnitSystem U, std::size_t N>
    requires(N > 0)
struct Reaction {
    std::array<int, N> reactants;
    std::array<int, N> products;
    int coefficient_change;
    ArrheniusRate<P, U> forward;
    BackwardRateModel backward_model;
    ArrheniusRate<P, U> backward;
    bool third_body;
    std::array<P, N> efficiencies;
};

// R reactions over the N species of a thermo model, in unit system U, as
// derived by make_reaction_mechanism: the reactions, the particle mass of
// every species, the temperature limits of rate evaluation, the standard
// pressure p0 of the equilibrium constants and the Boltzmann constant.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
    requires(N > 0 && R > 0)
struct ReactionMechanism {
    using parameter_type = P;
    using units = U;

    std::array<Reaction<P, U, N>, R> reactions;
    std::array<ParticleMass<P, U>, N> mass;
    Temperature<P, U> temperature_lower;
    Temperature<P, U> temperature_upper;
    Pressure<P, U> standard_pressure;
    BoltzmannConstant<P, U> boltzmann;
};

namespace detail
{

// The scales with double values, for converting rate constants whose scale
// factors may lie outside the range of P.
template <std::floating_point P> UnitScales<double, si> double_scales(const UnitScales<P, si>&)
{ return {}; }

template <std::floating_point P>
UnitScales<double, scaled> double_scales(const UnitScales<P, scaled>& s)
{
    return {
      .mass = static_cast<double>(s.mass),
      .density = static_cast<double>(s.density),
      .temperature = static_cast<double>(s.temperature),
      .length = static_cast<double>(s.length),
    };
}

// Index of name in names; names.size() when it is not there.
template <std::size_t N>
std::size_t species_index(const std::array<std::string_view, N>& names, std::string_view name)
{
    for (std::size_t i = 0; i < N; ++i) {
        if (names[i] == name) { return i; }
    }
    return N;
}

// Checks the entry count of a list and the names of its first count
// entries: each known and listed once. Returns the species indices.
template <std::size_t N, std::size_t K>
std::array<std::size_t, K> checked_indices(const SpeciesCoefficients<K>& list,
                                           const std::array<std::string_view, N>& names,
                                           std::string_view role,
                                           int min_count)
{
    if (list.count < min_count || list.count > static_cast<int>(K)) {
        throw ParameterError(std::format("number of {} must be between {} and {}, got {}",
                                         role,
                                         min_count,
                                         K,
                                         list.count));
    }
    std::array<std::size_t, K> index {};
    std::array<bool, N> seen {};
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const char* name = list.entries[e].name;
        if (name == nullptr) {
            throw ParameterError(std::format("{} entry {} has no name", role, e));
        }
        index[e] = species_index(names, name);
        if (index[e] == N) {
            throw ParameterError(std::format("unknown species \"{}\" in the {}", name, role));
        }
        if (seen[index[e]]) {
            throw ParameterError(
              std::format("species \"{}\" is listed twice in the {}", name, role));
        }
        seen[index[e]] = true;
    }
    return index;
}

// Stoichiometric coefficients per species from a list of at least one
// species. Each must be a positive integer.
template <std::size_t N, std::size_t K>
std::array<int, N> stoichiometric_coefficients(const SpeciesCoefficients<K>& list,
                                               const std::array<std::string_view, N>& names,
                                               std::string_view role)
{
    const std::array<std::size_t, K> index = checked_indices(list, names, role, 1);
    std::array<int, N> out {};
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const double nu = list.entries[e].value;
        if (!(nu >= 1.0) || !(nu <= static_cast<double>(std::numeric_limits<int>::max())) ||
            nu != std::floor(nu)) {
            throw ParameterError(
              std::format("coefficient of \"{}\" in the {} must be a positive integer, got {}",
                          list.entries[e].name,
                          role,
                          nu));
        }
        out[index[e]] = static_cast<int>(nu);
    }
    return out;
}

// Third-body efficiencies per species: 1 unless listed. Each listed value
// must be finite and non-negative.
template <std::floating_point P, std::size_t N, std::size_t K>
std::array<P, N> third_body_efficiencies(const SpeciesCoefficients<K>& list,
                                         const std::array<std::string_view, N>& names)
{
    const std::array<std::size_t, K> index =
      checked_indices(list, names, "third-body efficiencies", 0);
    std::array<P, N> out {};
    out.fill(P(1));
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const double eff = list.entries[e].value;
        if (!std::isfinite(eff) || eff < 0.0) {
            throw ParameterError(
              std::format("third-body efficiency of \"{}\" must be finite and non-negative, got {}",
                          list.entries[e].name,
                          eff));
        }
        out[index[e]] = static_cast<P>(eff);
    }
    return out;
}

// Rate constant of the given order in model units. A must be finite and
// positive, n finite and theta finite and non-negative. A per mole in SI
// becomes A / N_A^(order - 1) per particle; with T_SI = tau T, the model A
// is that times tau^n divided by the scale of (length^3 / count)^(order - 1)
// / time. The conversion is done in double.
template <std::floating_point P, UnitSystem U>
ArrheniusRate<P, U> arrhenius_rate(const ArrheniusInput& input,
                                   int order,
                                   const UnitScales<P, U>& scales,
                                   std::string_view role)
{
    if (!std::isfinite(input.A) || !(input.A > 0.0)) {
        throw ParameterError(
          std::format("{} rate: A must be finite and positive, got {}", role, input.A));
    }
    if (!std::isfinite(input.n)) {
        throw ParameterError(std::format("{} rate: n must be finite, got {}", role, input.n));
    }
    if (!std::isfinite(input.activation_temperature) || input.activation_temperature < 0.0) {
        throw ParameterError(
          std::format("{} rate: activation temperature must be finite and non-negative, got {}",
                      role,
                      input.activation_temperature));
    }
    const UnitScales<double, U> s = double_scales(scales);
    const Dimension dimension {
      .length = 3 * (order - 1),
      .time = -1,
      .count = -(order - 1),
    };
    const double tau = scale_factor<kind::temperature>(s);
    const double per_particle =
      input.A / std::pow(constants::avogadro<double>, static_cast<double>(order - 1));
    const auto a =
      static_cast<P>(per_particle * std::pow(tau, input.n) / scale_factor(dimension, s));
    if (!std::isfinite(a) || !(a > P(0))) {
        throw ParameterError(
          std::format("{} rate: A = {} is not representable in the model's units", role, input.A));
    }
    return {
      .A = a,
      .n = static_cast<P>(input.n),
      .activation_temperature =
        Temperature<P, U> {static_cast<P>(input.activation_temperature / tau)},
    };
}

template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t K>
Reaction<P, U, N> make_reaction(const ReactionInput<K>& input,
                                const std::array<std::string_view, N>& names,
                                const UnitScales<P, U>& scales)
{
    Reaction<P, U, N> r {};
    r.reactants = stoichiometric_coefficients(input.reactants, names, "reactants");
    r.products = stoichiometric_coefficients(input.products, names, "products");
    r.third_body = input.third_body;
    if (input.third_body) {
        r.efficiencies = third_body_efficiencies<P>(input.efficiencies, names);
    } else if (input.efficiencies.count != 0) {
        throw ParameterError("third-body efficiencies are given without a third body");
    }
    int reactant_order = 0;
    int product_order = 0;
    for (std::size_t i = 0; i < N; ++i) {
        reactant_order += r.reactants[i];
        product_order += r.products[i];
    }
    r.coefficient_change = product_order - reactant_order;
    const int third_body_order = input.third_body ? 1 : 0;
    r.forward = arrhenius_rate(input.forward, reactant_order + third_body_order, scales, "forward");
    r.backward_model = input.backward_model;
    if (input.backward_model == BackwardRateModel::arrhenius) {
        r.backward =
          arrhenius_rate(input.backward, product_order + third_body_order, scales, "backward");
    } else if (input.backward_model != BackwardRateModel::equilibrium) {
        throw ParameterError("unknown backward rate model");
    }
    return r;
}

}  // namespace detail

// Mechanism of the reactions of input over the N species of a thermo model,
// named in its species order (for example the name of each Nasa9Species),
// with their particle masses from species, in the unit system of the
// scales, with standard pressure p0. Rate constants are converted to
// per-particle model units, activation temperatures and the temperature
// limits to the temperature unit of U. Throws ParameterError when:
// - species.n is not N, or a name is listed twice;
// - the temperature limits are not finite with 0 < lower < upper, or p0 is
//   not finite and positive;
// and, with the reaction index and equation in the message, when:
// - the reactants or products are empty, name an unknown species or one
//   species twice, or a coefficient is not a positive integer;
// - a third-body efficiency names an unknown species or is negative or not
//   finite, or efficiencies are given without a third body;
// - an Arrhenius A is not finite and positive (also after conversion), n is
//   not finite, or an activation temperature is negative or not finite.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, std::size_t K>
    requires(N > 0 && R > 0 && K > 0)
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const std::array<std::string_view, N>& names,
                                                      const SpeciesData<P, U, N>& species,
                                                      const UnitScales<P, U>& scales,
                                                      Pressure<P, U> standard_pressure)
{
    if (species.n != static_cast<int>(N)) {
        throw ParameterError(
          std::format("species data must have {} species, got {}", N, species.n));
    }
    for (std::size_t i = 0; i < N; ++i) {
        if (detail::species_index(names, names[i]) != i) {
            throw ParameterError(std::format("species \"{}\" is named twice", names[i]));
        }
    }
    const double lower = input.temperature_lower;
    const double upper = input.temperature_upper;
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(lower > 0.0) || !(lower < upper)) {
        throw ParameterError(
          std::format("temperature limits must satisfy 0 < lower < upper, got {} and {}",
                      lower,
                      upper));
    }
    const P p0 = standard_pressure.value();
    if (!std::isfinite(p0) || !(p0 > P(0))) {
        throw ParameterError(
          std::format("standard pressure must be finite and positive, got {}", p0));
    }
    ReactionMechanism<P, U, N, R> mech {};
    for (std::size_t k = 0; k < R; ++k) {
        const ReactionInput<K>& reaction = input.reactions[k];
        try {
            mech.reactions[k] = detail::make_reaction(reaction, names, scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("reaction {} ({}): {}",
                                             k,
                                             reaction.equation != nullptr ? reaction.equation : "",
                                             error.what()));
        }
    }
    mech.mass = species.mass;
    mech.temperature_lower = convert<U>(Temperature<P, si> {static_cast<P>(lower)}, scales);
    mech.temperature_upper = convert<U>(Temperature<P, si> {static_cast<P>(upper)}, scales);
    mech.standard_pressure = standard_pressure;
    mech.boltzmann = boltzmann_constant(scales);
    return mech;
}

// As above, with the standard pressure 1 bar (1e5 Pa).
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, std::size_t K>
    requires(N > 0 && R > 0 && K > 0)
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const std::array<std::string_view, N>& names,
                                                      const SpeciesData<P, U, N>& species,
                                                      const UnitScales<P, U>& scales)
{
    return make_reaction_mechanism(input,
                                   names,
                                   species,
                                   scales,
                                   convert<U>(Pressure<P, si> {P(1.0e5)}, scales));
}

// A thermo model of N species whose species_properties(gas, temp) gives the
// SpeciesProperties of every species at a temperature with value type T.
template <typename G, typename T, std::size_t N>
concept HasSpeciesProperties =
  GasModelType<G> && Scalar<T> && requires(const G& gas, Temperature<T, typename G::units> temp) {
      {
          species_properties(gas, temp)
      } -> std::convertible_to<SpeciesProperties<T, typename G::units, N>>;
  };

// Forward and backward rates of progress of R reactions, in particles per
// unit volume and time.
template <Scalar T, UnitSystem U, std::size_t R>
    requires(R > 0)
struct RatesOfProgress {
    std::array<NumberRateDensity<T, U>, R> forward;
    std::array<NumberRateDensity<T, U>, R> backward;
};

// A thermo model and value type for which rates_of_progress and
// production_rates of a mechanism of N species with parameter type P can be
// evaluated.
template <typename G, typename T, typename P, std::size_t N>
concept ReactingModel = HasSpeciesProperties<G, T, N> && std::floating_point<P> &&
                        ScalarWith<T, P> && HasExp<T> && HasLog<T>;

// Rates of progress of every reaction of mech in a cell of density rho and
// translational temperature temp, with the normalised mass fractions of
// mix. The rate constants and the Gibbs energies of the equilibrium
// constants are evaluated at temp clamped to [temperature_lower,
// temperature_upper] of mech (on its value part; a clamped temperature has
// derivative zero). The Arrhenius expressions compute on the values of the
// quantities, since A has no fixed kind. The gas model is called once,
// through species_properties. The status of mix is not checked.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
GASLIB_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const V t_value = ScalarTraits<T>::value(temp.value());
    const auto lower = static_cast<V>(mech.temperature_lower.value());
    const auto upper = static_cast<V>(mech.temperature_upper.value());
    T t = temp.value();
    if (t_value < lower) {
        t = T(lower);
    } else if (upper < t_value) {
        t = T(upper);
    }

    const SpeciesProperties<T, U, N> props = species_properties(gas, Temperature<T, U> {t});
    const T ln_t = math::log(t);
    const T inv_t = T(V(1)) / t;
    const T inv_kt = T(V(1)) / T(mech.boltzmann.value() * t);
    // Logarithm of the standard-state number density p0 / (k_B T).
    const T ln_n0 = math::log(T(mech.standard_pressure.value() * inv_kt));

    // Number densities n_i and standard Gibbs energies per particle over
    // k_B T.
    std::array<T, N> n {};
    std::array<T, N> g {};
    for (std::size_t i = 0; i < N; ++i) {
        const P m = mech.mass[i].value();
        n[i] = T(rho.value() * mix.fractions[i].value() / m);
        g[i] = T(props.gibbs_energy[i].value() * m * inv_kt);
    }

    RatesOfProgress<T, U, R> out {};
    for (std::size_t k = 0; k < R; ++k) {
        const Reaction<P, U, N>& r = mech.reactions[k];
        // ln(k_f / A) = n ln T - theta / T.
        const T forward_exponent =
          T(r.forward.n * ln_t) - T(r.forward.activation_temperature.value() * inv_t);
        const T k_f = T(r.forward.A * math::exp(forward_exponent));
        T k_b = k_f;
        if (r.backward_model == BackwardRateModel::equilibrium) {
            // k_b = k_f / K_c = A exp(ln(k_f / A) + sum dnu_i g_i m_i / (k_B T)
            //                         - sum(dnu) ln(p0 / (k_B T))).
            T exponent = forward_exponent - T(V(r.coefficient_change) * ln_n0);
            for (std::size_t i = 0; i < N; ++i) {
                const int dnu = r.products[i] - r.reactants[i];
                if (dnu != 0) { exponent = exponent + T(V(dnu) * g[i]); }
            }
            k_b = T(r.forward.A * math::exp(exponent));
        } else {
            k_b = T(r.backward.A * math::exp(T(r.backward.n * ln_t) -
                                             T(r.backward.activation_temperature.value() * inv_t)));
        }
        T q_f = k_f;
        T q_b = k_b;
        if (r.third_body) {
            T m = T(r.efficiencies[0] * n[0]);
            for (std::size_t i = 1; i < N; ++i) { m = m + T(r.efficiencies[i] * n[i]); }
            q_f = q_f * m;
            q_b = q_b * m;
        }
        for (std::size_t i = 0; i < N; ++i) {
            for (int c = 0; c < r.reactants[i]; ++c) { q_f = q_f * n[i]; }
            for (int c = 0; c < r.products[i]; ++c) { q_b = q_b * n[i]; }
        }
        out.forward[k] = NumberRateDensity<T, U> {q_f};
        out.backward[k] = NumberRateDensity<T, U> {q_b};
    }
    return out;
}

// Mass production rate of every species, omega_i = m_i sum_r (nu''_ir -
// nu'_ir) (q_f,r - q_b,r), from the rates of progress of rates_of_progress
// with the same arguments.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
GASLIB_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const RatesOfProgress<T, U, R> q = rates_of_progress(mech, gas, mix, rho, temp);
    std::array<T, N> net {};
    net.fill(T(V(0)));
    for (std::size_t k = 0; k < R; ++k) {
        const Reaction<P, U, N>& r = mech.reactions[k];
        const T rate = q.forward[k].value() - q.backward[k].value();
        for (std::size_t i = 0; i < N; ++i) {
            const int dnu = r.products[i] - r.reactants[i];
            if (dnu != 0) { net[i] = net[i] + T(V(dnu) * rate); }
        }
    }
    std::array<MassProductionRate<T, U>, N> out {};
    for (std::size_t i = 0; i < N; ++i) {
        out[i] = MassProductionRate<T, U> {T(mech.mass[i].value() * net[i])};
    }
    return out;
}

}  // namespace gaslib

#endif  // GASLIB_REACTIONS_HPP
