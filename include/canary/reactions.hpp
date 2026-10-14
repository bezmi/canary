#ifndef CANARY_REACTIONS_HPP
#define CANARY_REACTIONS_HPP

// Finite-rate chemistry of a mixture: a mechanism of reversible reactions
// with Arrhenius forward rates, backward rates from the equilibrium constant
// or from their own Arrhenius fits, and optional third bodies; per cell, the
// rates of progress of every reaction and the mass production rate of every
// species, at one temperature, at two (the translational temperature T and
// the vibrational-electronic temperature T_ve of a two-temperature model) or
// at T and the M mode temperatures T_m of a multi-temperature model.
//
// Concentrations are number densities n_i = rho y_i / m_i (particles per
// unit volume), so that the same functions work in SI and in scaled units.
// T_l is T clamped to the mechanism's temperature limits. Each rate has a
// second temperature T_o: T_ve with two temperatures, the temperature of
// the mode the rate's other_temperature resolves to with M modes, and T_l
// with one temperature; T_o is not clamped. An Arrhenius rate constant with
// translational exponent s is evaluated at the rate-controlling temperature
// T_c = T_l^s T_o^(1 - s):
//   k(T_c) = A T_c^n exp(-theta / T_c).
// For a reaction with reactant coefficients nu'_i and product coefficients
// nu''_i:
//   k_f = k(T_c) with the forward A, n, theta and s,
//   [M] = sum_i eps_i n_i with a third body, 1 without,
//   q_f = k_f [M] prod_i n_i^nu'_i,   q_b = k_b [M] prod_i n_i^nu''_i,
//   omega_i = m_i sum_r (nu''_ir - nu'_ir) (q_f,r - q_b,r).
// The backward rate constant is either an Arrhenius rate of its own, at its
// own T_c, or k_b = k_f(T_eq) / K_c(T_eq) at the equilibrium temperature
// T_eq: T_l, or the backward rate's T_o when the backward input names an
// other_temperature. k_f(T_eq) is the forward Arrhenius expression at
// T_c = T_eq, and K_c is the equilibrium constant in concentration form
//   K_c = exp(-sum_i dnu_i g_i m_i / (k_B T_eq)) (p0 / (k_B T_eq))^(sum_i dnu_i),
// where dnu_i = nu''_i - nu'_i, g_i is the standard-state Gibbs energy per
// unit mass of species i at T_eq from the species_properties of the thermo
// model, and p0 the standard pressure. With T_eq = T_l, k_b does not depend
// on T_o.

#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <string_view>
#include <type_traits>

namespace canary
{

// Host input of an Arrhenius rate constant k = A T_c^n exp(-theta / T_c), in
// SI per mole: A in (m^3/mol)^(order - 1) / s with T_c in K, where order is
// the sum of the coefficients of the reaction side the rate multiplies, plus
// 1 with a third body; n dimensionless; activation_temperature theta in K;
// translational_exponent s in [0, 1], the exponent of the translational
// temperature in the rate-controlling temperature T_c = T^s T_o^(1 - s)
// (1: the rate depends on T only); other_temperature names the second
// temperature T_o for a thermo model with energy modes (see
// make_reaction_mechanism with a gas model): "electron" for the temperature
// of the mode holding the free electrons, or a species name for the
// temperature of the mode holding that species' vibration; null for the
// default. The backward input of a BackwardRateModel::equilibrium rate uses
// only other_temperature: when it is not null, the equilibrium backward rate
// is evaluated at that T_o.
struct ArrheniusInput {
    double A;
    double n;
    double activation_temperature;
    double translational_exponent = 1.0;
    const char* other_temperature = nullptr;
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

// Arrhenius rate constant k = A T_c^n exp(-theta / T_c) in model units: T_c
// and theta in the temperature unit of U, and A per particle in
// (length^3 / count)^(order - 1) / time of U, times the temperature scale
// raised to n. A and n are plain values, since the dimension of A depends
// on the reaction order. translational_exponent is s of
// T_c = T^s T_o^(1 - s); other_mode is the index of the energy mode whose
// temperature is T_o in an evaluation with M mode temperatures (0 for a
// mechanism built without a gas model).
template <std::floating_point P, UnitSystem U> struct ArrheniusRate {
    P A;
    P n;
    Temperature<P, U> activation_temperature;
    P translational_exponent;
    ModeIndex other_mode;
};

// One reaction of a mechanism over N species, indexed as the species of the
// thermo model: reactant and product coefficients, the sum of the
// coefficient changes sum_i (products[i] - reactants[i]), the forward rate,
// the backward rate model and rate (used with BackwardRateModel::arrhenius),
// whether an equilibrium backward rate is evaluated at the T_o of
// backward.other_mode instead of T_l, whether a third body takes part and
// the efficiency of every species (used with a third body).
template <std::floating_point P, UnitSystem U, std::size_t N>
    requires(N > 0)
struct Reaction {
    std::array<int, N> reactants;
    std::array<int, N> products;
    int coefficient_change;
    ArrheniusRate<P, U> forward;
    BackwardRateModel backward_model;
    ArrheniusRate<P, U> backward;
    bool equilibrium_at_other_temperature;
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

// Checks the entry count of a list and the names of its first count
// entries: each known and listed once. Returns the species indices.
template <std::size_t N, std::size_t K>
std::array<SpeciesIndex, K> checked_indices(const SpeciesCoefficients<K>& list,
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
    std::array<SpeciesIndex, K> index {};
    std::array<bool, N> seen {};
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const char* name = list.entries[e].name;
        if (name == nullptr) {
            throw ParameterError(std::format("{} entry {} has no name", role, e));
        }
        const auto found = species_index(names, name);
        if (!found) {
            throw ParameterError(std::format("unknown species \"{}\" in the {}", name, role));
        }
        index[e] = *found;
        if (seen[index[e].value()]) {
            throw ParameterError(
              std::format("species \"{}\" is listed twice in the {}", name, role));
        }
        seen[index[e].value()] = true;
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
    const std::array<SpeciesIndex, K> index = checked_indices(list, names, role, 1);
    std::array<int, N> out {};
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const double nu = list.entries[e].value;
        if (!(nu >= 1.0) || !(nu < std::ldexp(1.0, std::numeric_limits<int>::digits)) ||
            nu != std::floor(nu)) {
            throw ParameterError(
              std::format("coefficient of \"{}\" in the {} must be a positive integer, got {}",
                          list.entries[e].name,
                          role,
                          nu));
        }
        out[index[e].value()] = static_cast<int>(nu);
    }
    return out;
}

// Third-body efficiencies per species: 1 unless listed. Each listed value
// must be finite and non-negative.
template <std::floating_point P, std::size_t N, std::size_t K>
std::array<P, N> third_body_efficiencies(const SpeciesCoefficients<K>& list,
                                         const std::array<std::string_view, N>& names)
{
    const std::array<SpeciesIndex, K> index =
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
        out[index[e].value()] = static_cast<P>(eff);
    }
    return out;
}

// Host-side reaction orders include an optional third body. Accumulation
// is checked before addition; stored integer coefficients remain unchanged.
template <std::size_t N>
std::intmax_t
  reaction_order(const std::array<int, N>& coefficients, bool third_body, std::string_view role)
{
    std::intmax_t order = third_body ? 1 : 0;
    for (const int coefficient : coefficients) {
        if (coefficient < 0 || order > std::numeric_limits<std::intmax_t>::max() - coefficient) {
            throw ParameterError(std::format("{} reaction order is not representable", role));
        }
        order += coefficient;
    }
    return order;
}

// The integer exponents must be representable before constructing Dimension.
inline Dimension rate_dimension(std::intmax_t order, std::string_view role)
{
    constexpr std::intmax_t maximum_order = std::numeric_limits<int>::max() / 3 + 1;
    if (order < 1 || order > maximum_order) {
        throw ParameterError(std::format(
          "{} rate: reaction order {} is outside the representable dimension range [1, {}]",
          role,
          order,
          maximum_order));
    }
    const auto power = static_cast<int>(order - 1);
    return {.length = 3 * power, .time = -1, .count = -power};
}

// Rate constant of the given order in model units. A must be finite and
// positive, n finite, theta finite and non-negative and s finite in [0, 1].
// A per mole in SI
// becomes A / N_A^(order - 1) per particle; with T_SI = tau T, the model A
// is that times tau^n divided by the scale of (length^3 / count)^(order - 1)
// / time. Conversion uses at least double precision and retains long double.
template <std::floating_point P, UnitSystem U>
ArrheniusRate<P, U> arrhenius_rate(const ArrheniusInput& input,
                                   std::intmax_t order,
                                   const UnitScales<P, U>& scales,
                                   std::string_view role)
{
    using W = std::common_type_t<double, P>;
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
    const double s_t = input.translational_exponent;
    if (!std::isfinite(s_t) || s_t < 0.0 || s_t > 1.0) {
        throw ParameterError(
          std::format("{} rate: translational exponent must be finite and between 0 and 1, got {}",
                      role,
                      s_t));
    }
    const auto& s = scales;
    const Dimension dimension = rate_dimension(order, role);
    const auto tau = scale_factor<kind::temperature>(s);
    const auto per_particle =
      input.A / std::pow(constants::avogadro<double>, static_cast<W>(order - 1));
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
      .translational_exponent = static_cast<P>(s_t),
      .other_mode = ModeIndex {0},
    };
}

// The other modes of a mechanism built without a gas model: every T_o is
// the single second temperature of a two-temperature evaluation (mode 0).
struct SingleOtherTemperature {
    template <std::size_t K>
    ModeIndex operator()(const ArrheniusInput& /*rate*/,
                         const SpeciesCoefficients<K>& /*reactants*/,
                         std::string_view /*role*/) const
    { return ModeIndex {0}; }
};

// Reaction of input with the other mode of each rate from
// other_mode(rate, reactants, role), which throws ParameterError for a rate
// whose other temperature does not resolve.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t K, typename OtherMode>
Reaction<P, U, N> make_reaction(const ReactionInput<K>& input,
                                const std::array<std::string_view, N>& names,
                                const UnitScales<P, U>& scales,
                                const OtherMode& other_mode)
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
    const std::intmax_t reactant_order = reaction_order(r.reactants, input.third_body, "reactant");
    const std::intmax_t product_order = reaction_order(r.products, input.third_body, "product");
    const std::intmax_t change = product_order - reactant_order;
    if (change < std::numeric_limits<int>::min() || change > std::numeric_limits<int>::max()) {
        throw ParameterError(
          std::format("net stoichiometric coefficient {} is not representable", change));
    }
    r.coefficient_change = static_cast<int>(change);
    r.forward = arrhenius_rate(input.forward, reactant_order, scales, "forward");
    r.forward.other_mode = other_mode(input.forward, input.reactants, "forward");
    r.backward_model = input.backward_model;
    r.equilibrium_at_other_temperature = false;
    if (input.backward_model == BackwardRateModel::arrhenius) {
        r.backward = arrhenius_rate(input.backward, product_order, scales, "backward");
        r.backward.other_mode = other_mode(input.backward, input.reactants, "backward");
    } else if (input.backward_model == BackwardRateModel::equilibrium) {
        if (input.backward.other_temperature != nullptr) {
            r.equilibrium_at_other_temperature = true;
            r.backward.other_mode = other_mode(input.backward, input.reactants, "backward");
        }
    } else {
        throw ParameterError("unknown backward rate model");
    }
    return r;
}

// Builds the reactions of input into mech; see make_reaction_mechanism.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t K,
          typename OtherMode>
ReactionMechanism<P, U, N, R> build_mechanism(const MechanismInput<R, K>& input,
                                              const std::array<std::string_view, N>& names,
                                              const SpeciesData<P, U, N>& species,
                                              const UnitScales<P, U>& scales,
                                              Pressure<P, U> standard_pressure,
                                              const OtherMode& other_mode)
{
    if (species.n != static_cast<int>(N)) {
        throw ParameterError(
          std::format("species data must have {} species, got {}", N, species.n));
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (species_index(names, names[i.value()]) != OptionalSpeciesIndex {i}) {
            throw ParameterError(std::format("species \"{}\" is named twice", names[i.value()]));
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
            mech.reactions[k] = make_reaction(reaction, names, scales, other_mode);
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
//   not finite, an activation temperature is negative or not finite, or a
//   translational exponent is not finite or outside [0, 1] (the backward
//   rate is checked only with BackwardRateModel::arrhenius).
// names, species and the thermodynamic model used during rate evaluation
// must use the same species order. Construction checks names and counts,
// but SpeciesData carries no names to verify this ordering.
// Every rate's other mode is 0, whatever its other_temperature: each T_o is
// the single second temperature T_ve of a two-temperature evaluation, and a
// backward input of an equilibrium rate with an other_temperature gives an
// equilibrium backward rate at T_ve.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, std::size_t K>
    requires(N > 0 && R > 0 && K > 0)
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const std::array<std::string_view, N>& names,
                                                      const SpeciesData<P, U, N>& species,
                                                      const UnitScales<P, U>& scales,
                                                      Pressure<P, U> standard_pressure)
{
    return detail::build_mechanism(input,
                                   names,
                                   species,
                                   scales,
                                   standard_pressure,
                                   detail::SingleOtherTemperature {});
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

// A thermo model with energy modes: G::mode_count modes, mode m with its
// own temperature T_m. species_vibration_mode(gas, i) is the mode that holds
// the vibration of species i (in the model's species order) and
// electron_energy_mode(gas) the mode that holds the energy of the free
// electrons; each returns absence when there is none.
template <typename G>
concept HasEnergyModes = GasModelType<G> && requires(const G& gas, SpeciesIndex i) {
    { G::mode_count } -> std::convertible_to<std::size_t>;
    { species_vibration_mode(gas, i) } -> std::same_as<OptionalModeIndex>;
    { electron_energy_mode(gas) } -> std::same_as<OptionalModeIndex>;
};

namespace detail
{

// The other modes of a mechanism built with the gas model *gas, over the
// species *names (in the species order of *gas); see
// make_reaction_mechanism. Both must outlive the builder call.
template <HasEnergyModes G, std::size_t N> struct ModelOtherTemperature {
    const G* gas;
    const std::array<std::string_view, N>* names;

    template <std::size_t K>
    ModeIndex operator()(const ArrheniusInput& rate,
                         const SpeciesCoefficients<K>& reactants,
                         std::string_view role) const
    {
        const char* name = rate.other_temperature;
        if (name == nullptr) {
            if (G::mode_count == 1) { return ModeIndex {0}; }
            for (std::size_t e = 0; e < static_cast<std::size_t>(reactants.count); ++e) {
                const auto index = species_index(*names, reactants.entries[e].name);
                if (!index) { throw ParameterError("unknown reactant species for mode selection"); }
                const OptionalModeIndex mode = species_vibration_mode(*gas, *index);
                if (mode.has_value()) { return *mode; }
            }
            if (rate.translational_exponent != 1.0) {
                throw ParameterError(std::format(
                  "{} rate: no reactant has a vibration component to give the other temperature",
                  role));
            }
            return ModeIndex {0};
        }
        if (std::string_view {name} == "electron") {
            const OptionalModeIndex mode = electron_energy_mode(*gas);
            if (!mode) {
                throw ParameterError(std::format(
                  "{} rate: other temperature \"electron\": the model has no free electrons",
                  role));
            }
            return *mode;
        }
        const auto i = species_index(*names, name);
        if (!i) {
            throw ParameterError(
              std::format(R"({} rate: other temperature "{}" is not "electron" or a species)",
                          role,
                          name));
        }
        const OptionalModeIndex mode = species_vibration_mode(*gas, *i);
        if (!mode) {
            throw ParameterError(std::format(
              "{} rate: other temperature \"{}\": the species has no vibration component",
              role,
              name));
        }
        return *mode;
    }
};

}  // namespace detail

// As make_reaction_mechanism with a standard pressure, with each rate's
// other temperature resolved against the energy modes of gas, whose
// species are those of names in the same order. The other mode of a rate
// is:
// - for other_temperature "electron": electron_energy_mode(gas);
// - for a species name: species_vibration_mode(gas) of that species;
// - without other_temperature: 0 when gas has one mode; otherwise the
//   vibration mode of the first species in the reactant list (in input
//   order) that has one, or 0 when none has one and s = 1.
// The backward rate of BackwardRateModel::arrhenius resolves the same way,
// with the reactant list; the backward input of an equilibrium rate is
// resolved only when it names an other_temperature. In addition to the
// errors above, throws ParameterError, with the reaction index and
// equation, when other_temperature is "electron" and gas has no electron
// mode, names neither "electron" nor a species, or names a species without
// a vibration mode, or when it is null, gas has more than one mode, s is not
// 1 and no reactant has a vibration mode.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t K,
          HasEnergyModes G>
    requires(N > 0 && R > 0 && K > 0) && std::same_as<typename G::units, U>
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const std::array<std::string_view, N>& names,
                                                      const SpeciesData<P, U, N>& species,
                                                      const UnitScales<P, U>& scales,
                                                      Pressure<P, U> standard_pressure,
                                                      const G& gas)
{
    return detail::build_mechanism(
      input,
      names,
      species,
      scales,
      standard_pressure,
      detail::ModelOtherTemperature<G, N> {.gas = &gas, .names = &names});
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

namespace detail
{

// temp clamped to [temperature_lower, temperature_upper] of mech, on its
// value part; a clamped temperature has derivative zero.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE T clamped_temperature(const ReactionMechanism<P, U, N, R>& mech,
                                         Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const V t_value = ScalarTraits<T>::value(temp.value());
    const auto lower = static_cast<V>(mech.temperature_lower.value());
    const auto upper = static_cast<V>(mech.temperature_upper.value());
    if (t_value < lower) { return T(lower); }
    if (upper < t_value) { return T(upper); }
    return temp.value();
}

// ln(k / A) of an Arrhenius rate at T_l, from ln T_l and 1 / T_l:
// n ln T_l - theta / T_l.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T arrhenius_exponent(const ArrheniusRate<P, U>& rate, T ln_t, T inv_t)
{ return T(rate.n * ln_t) - T(rate.activation_temperature.value() * inv_t); }

// ln(k / A) of an Arrhenius rate at T_c = T_l^s T_ve^(1 - s), from ln T_l,
// 1 / T_l and ln_ratio = ln(T_ve / T_l): with w = (1 - s) ln_ratio,
// ln T_c = ln T_l + w and 1 / T_c = exp(-w) / T_l. w = 0 gives the value of
// the expression at T_l.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T
  arrhenius_exponent(const ArrheniusRate<P, U>& rate, T ln_t, T inv_t, T ln_ratio)
{
    const T w = T((P(1) - rate.translational_exponent) * ln_ratio);
    return T(rate.n * T(ln_t + w)) -
           T(rate.activation_temperature.value() * T(inv_t * math::exp(T(-w))));
}

// Index into the M other temperatures of rates_of_progress_at of a rate
// with the given other mode: 0 with one other temperature (T_ve of a
// two-temperature evaluation, or the only mode), the mode otherwise.
template <std::size_t M> CANARY_HOST_DEVICE constexpr ModeIndex other_index(ModeIndex mode)
{
    if constexpr (M == 1) {
        return ModeIndex {0};
    } else {
        // NOLINTNEXTLINE(misc-static-assert): run-time mode index
        CANARY_ASSERT(mode.value() < M);
        return mode;
    }
}

// ln(k / A) of an Arrhenius rate in rates_of_progress_at: at T_l with no
// other temperature (M = 0); otherwise at T_c for s != 1, from ln_ratio =
// ln(T_o / T_l) of the rate's other temperature, and at T_l for s = 1, so
// that a rate that depends on T only computes no exponential of w.
template <std::size_t M, Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T
  rate_exponent(const ArrheniusRate<P, U>& rate, T ln_t, T inv_t, const std::array<T, M>& ln_ratio)
{
    if constexpr (M > 0) {
        if (rate.translational_exponent != P(1)) {
            return arrhenius_exponent(rate,
                                      ln_t,
                                      inv_t,
                                      ln_ratio[other_index<M>(rate.other_mode).value()]);
        }
    }
    return arrhenius_exponent(rate, ln_t, inv_t);
}

// ln(p0 / (k_B T)) and the standard Gibbs energies per particle over k_B T,
// g_i m_i / (k_B T), at a temperature T, for the equilibrium constants.
template <Scalar T, std::size_t N> struct EquilibriumTerms {
    T ln_n0;
    std::array<T, N> g;
};

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE EquilibriumTerms<T, N>
  equilibrium_terms(const ReactionMechanism<P, U, N, R>& mech, const G& gas, T t)
{
    using V = ScalarTraits<T>::value_type;
    const SpeciesProperties<T, U, N> props = species_properties(gas, Temperature<T, U> {t});
    const T inv_kt = T(V(1)) / T(mech.boltzmann.value() * t);
    EquilibriumTerms<T, N> out {};
    out.ln_n0 = math::log(T(mech.standard_pressure.value() * inv_kt));
    for (const SpeciesIndex i : species_indices(0, N)) {
        out.g[i.value()] =
          T(props.gibbs_energy[i.value()].value() * mech.mass[i.value()].value() * inv_kt);
    }
    return out;
}

// ln(k_b / A_f) of an equilibrium backward rate, from ln(k_f / A_f) at the
// equilibrium temperature and the terms of equilibrium_terms there:
// ln(k_f / A_f) + sum dnu_i g_i m_i / (k_B T) - sum(dnu) ln(p0 / (k_B T)).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T equilibrium_exponent(const Reaction<P, U, N>& r,
                                          T forward_exponent,
                                          const EquilibriumTerms<T, N>& terms)
{
    using V = ScalarTraits<T>::value_type;
    T exponent = forward_exponent - T(V(r.coefficient_change) * terms.ln_n0);
    for (const SpeciesIndex i : species_indices(0, N)) {
        const int dnu = r.products[i.value()] - r.reactants[i.value()];
        if (dnu != 0) { exponent = exponent + T(V(dnu) * terms.g[i.value()]); }
    }
    return exponent;
}

// Rates of progress at the clamped translational temperature t (T_l) and
// the M other temperatures other (unclamped): T_ve of a two-temperature
// evaluation (M = 1) or the mode temperatures; with M = 0 every rate is
// evaluated at T_l. The gas model is called once at t and once at each
// other temperature at which some reaction has its equilibrium backward
// rate.
template <std::size_t M,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress_at(const ReactionMechanism<P, U, N, R>& mech,
                       const G& gas,
                       const NormalisedComposition<T, U, N>& mix,
                       Density<T, U> rho,
                       T t,
                       const std::array<T, M>& other)
{
    using V = ScalarTraits<T>::value_type;
    const T ln_t = math::log(t);
    const T inv_t = T(V(1)) / t;
    const EquilibriumTerms<T, N> at_t = equilibrium_terms(mech, gas, t);

    // Number densities n_i.
    std::array<T, N> n {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        n[i.value()] =
          T(rho.value() * mix.fractions[i.value()].value() / mech.mass[i.value()].value());
    }

    // ln(T_o / T_l) of every other temperature and, at each other
    // temperature that has equilibrium backward rates, ln T_o, 1 / T_o and
    // the equilibrium terms.
    std::array<T, M> ln_ratio {};
    std::array<T, M> ln_other {};
    std::array<T, M> inv_other {};
    std::array<EquilibriumTerms<T, N>, M> at_other {};
    if constexpr (M > 0) {
        std::array<bool, M> used {};
        for (const Reaction<P, U, N>& r : mech.reactions) {
            if (r.backward_model == BackwardRateModel::equilibrium &&
                r.equilibrium_at_other_temperature) {
                used[other_index<M>(r.backward.other_mode).value()] = true;
            }
        }
        for (const ModeIndex m : mode_indices(0, M)) {
            ln_ratio[m.value()] = math::log(T(other[m.value()] / t));
            if (!used[m.value()]) { continue; }
            ln_other[m.value()] = math::log(other[m.value()]);
            inv_other[m.value()] = T(V(1)) / other[m.value()];
            at_other[m.value()] = equilibrium_terms(mech, gas, other[m.value()]);
        }
    }

    RatesOfProgress<T, U, R> out {};
    for (std::size_t k = 0; k < R; ++k) {
        const Reaction<P, U, N>& r = mech.reactions[k];
        // ln(k_f / A) at T_l, and at the rate-controlling temperature.
        const T forward_exponent = arrhenius_exponent(r.forward, ln_t, inv_t);
        T controlled_exponent = forward_exponent;
        if constexpr (M > 0) {
            if (r.forward.translational_exponent != P(1)) {
                controlled_exponent =
                  arrhenius_exponent(r.forward,
                                     ln_t,
                                     inv_t,
                                     ln_ratio[other_index<M>(r.forward.other_mode).value()]);
            }
        }
        const T k_f = T(r.forward.A * math::exp(controlled_exponent));
        T k_b = k_f;
        if (r.backward_model == BackwardRateModel::equilibrium) {
            bool at_other_temperature = false;
            if constexpr (M > 0) { at_other_temperature = r.equilibrium_at_other_temperature; }
            T exponent = forward_exponent;
            if (at_other_temperature) {
                const ModeIndex m = other_index<M>(r.backward.other_mode);
                exponent = equilibrium_exponent(
                  r,
                  arrhenius_exponent(r.forward, ln_other[m.value()], inv_other[m.value()]),
                  at_other[m.value()]);
            } else {
                exponent = equilibrium_exponent(r, forward_exponent, at_t);
            }
            k_b = T(r.forward.A * math::exp(exponent));
        } else {
            k_b = T(r.backward.A * math::exp(rate_exponent<M>(r.backward, ln_t, inv_t, ln_ratio)));
        }
        T q_f = k_f;
        T q_b = k_b;
        if (r.third_body) {
            T m = T(r.efficiencies[0] * n[0]);
            for (const SpeciesIndex i : species_indices(1, N)) {
                m = m + T(r.efficiencies[i.value()] * n[i.value()]);
            }
            q_f = q_f * m;
            q_b = q_b * m;
        }
        for (const SpeciesIndex i : species_indices(0, N)) {
            for (int c = 0; c < r.reactants[i.value()]; ++c) { q_f = q_f * n[i.value()]; }
            for (int c = 0; c < r.products[i.value()]; ++c) { q_b = q_b * n[i.value()]; }
        }
        out.forward[k] = NumberRateDensity<T, U> {q_f};
        out.backward[k] = NumberRateDensity<T, U> {q_b};
    }
    return out;
}

// omega_i = m_i sum_r (nu''_ir - nu'_ir) (q_f,r - q_b,r) from the rates of
// progress q of mech.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  mass_production_rates(const ReactionMechanism<P, U, N, R>& mech,
                        const RatesOfProgress<T, U, R>& q)
{
    using V = ScalarTraits<T>::value_type;
    std::array<T, N> net {};
    net.fill(T(V(0)));
    for (std::size_t k = 0; k < R; ++k) {
        const Reaction<P, U, N>& r = mech.reactions[k];
        const T rate = q.forward[k].value() - q.backward[k].value();
        for (const SpeciesIndex i : species_indices(0, N)) {
            const int dnu = r.products[i.value()] - r.reactants[i.value()];
            if (dnu != 0) { net[i.value()] = net[i.value()] + T(V(dnu) * rate); }
        }
    }
    std::array<MassProductionRate<T, U>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        out[i.value()] =
          MassProductionRate<T, U> {T(mech.mass[i.value()].value() * net[i.value()])};
    }
    return out;
}

}  // namespace detail

// Rates of progress of every reaction of mech in a cell of density rho,
// translational temperature temp and vibrational-electronic temperature
// temp_ve, with the normalised mass fractions of mix. Every rate constant
// and the Gibbs energies of the equilibrium constants use T_l, temp clamped
// to [temperature_lower, temperature_upper] of mech (on its value part; a
// clamped temperature has derivative zero), except that temp_ve is T_o of
// every rate: it is not clamped and enters through the rate-controlling
// temperatures of the Arrhenius rates and the equilibrium backward rates
// at T_o (see the top of this file). The Arrhenius expressions compute on
// the values of the quantities, since A has no fixed kind. The gas model is
// called through species_properties at T_l, and at T_ve when a reaction
// has its equilibrium backward rate at T_o. The status of mix is not
// checked.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp,
                    Temperature<T, U> temp_ve)
{
    return detail::rates_of_progress_at<1>(mech,
                                           gas,
                                           mix,
                                           rho,
                                           detail::clamped_temperature(mech, temp),
                                           std::array<T, 1> {temp_ve.value()});
}

// Rates of progress with the M mode temperatures mode_temp of the cell
// (for example MultiTemperatureState::mode_temp): as above, with T_o of
// each rate the temperature of its other mode (ArrheniusRate::other_mode;
// with M = 1, the only mode). The gas model is called through
// species_properties at T_l and at each mode temperature at which a
// reaction has its equilibrium backward rate.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G,
          std::size_t M>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> && (M > 0)
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp,
                    const std::array<ModeTemperature<T, U>, M>& mode_temp)
{
    std::array<T, M> other {};
    for (const ModeIndex m : mode_indices(0, M)) {
        other[m.value()] = mode_temp[m.value()].value();
    }
    return detail::rates_of_progress_at<M>(mech,
                                           gas,
                                           mix,
                                           rho,
                                           detail::clamped_temperature(mech, temp),
                                           other);
}

// Rates of progress at one temperature: every T_o is T_l, so that every
// rate-controlling temperature and every equilibrium temperature is T_l.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp)
{
    return detail::rates_of_progress_at<0>(mech,
                                           gas,
                                           mix,
                                           rho,
                                           detail::clamped_temperature(mech, temp),
                                           std::array<T, 0> {});
}

// Mass production rate of every species, omega_i = m_i sum_r (nu''_ir -
// nu'_ir) (q_f,r - q_b,r), from the rates of progress of rates_of_progress
// with the same arguments.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp,
                   Temperature<T, U> temp_ve)
{
    return detail::mass_production_rates(mech,
                                         rates_of_progress(mech, gas, mix, rho, temp, temp_ve));
}

// As above with the M mode temperatures of the cell.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G,
          std::size_t M>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> && (M > 0)
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp,
                   const std::array<ModeTemperature<T, U>, M>& mode_temp)
{
    return detail::mass_production_rates(mech,
                                         rates_of_progress(mech, gas, mix, rho, temp, mode_temp));
}

// As above at one temperature.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp)
{ return detail::mass_production_rates(mech, rates_of_progress(mech, gas, mix, rho, temp)); }

namespace detail
{

template <> struct InputStorage<ArrheniusInput> {
    double A {};
    double n {};
    double activation_temperature {};
    double translational_exponent {};
    OwnedTemperatureLabel other_temperature {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const ArrheniusInput& input) :
        A(input.A), n(input.n), activation_temperature(input.activation_temperature),
        translational_exponent(input.translational_exponent),
        other_temperature(input.other_temperature)
    {
    }
    ArrheniusInput borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .A = A,
          .n = n,
          .activation_temperature = activation_temperature,
          .translational_exponent = translational_exponent,
          .other_temperature = other_temperature.c_str(),
        };
    }
};

template <> struct InputStorage<SpeciesCoefficient> {
    OwnedSpeciesName name {};
    double value {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const SpeciesCoefficient& input) : name(input.name), value(input.value) {}
    SpeciesCoefficient borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .name = name.c_str(),
          .value = value,
        };
    }
};

template <std::size_t K> struct InputStorage<SpeciesCoefficients<K>> {
    int count {};
    InputStorage<std::array<SpeciesCoefficient, K>> entries {};
    using Scratch = std::tuple<typename InputStorage<std::array<SpeciesCoefficient, K>>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const SpeciesCoefficients<K>& input) :
        count(input.count), entries(input.entries)
    {
    }
    SpeciesCoefficients<K> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .count = count,
          .entries = entries.borrow(std::get<0>(scratch)),
        };
    }
};

template <std::size_t K> struct InputStorage<ReactionInput<K>> {
    OwnedReactionEquation equation {};
    InputStorage<SpeciesCoefficients<K>> reactants {};
    InputStorage<SpeciesCoefficients<K>> products {};
    InputStorage<ArrheniusInput> forward {};
    BackwardRateModel backward_model {};
    InputStorage<ArrheniusInput> backward {};
    bool third_body {};
    InputStorage<SpeciesCoefficients<K>> efficiencies {};
    using Scratch = std::tuple<typename InputStorage<SpeciesCoefficients<K>>::Scratch,
                               typename InputStorage<SpeciesCoefficients<K>>::Scratch,
                               typename InputStorage<ArrheniusInput>::Scratch,
                               typename InputStorage<ArrheniusInput>::Scratch,
                               typename InputStorage<SpeciesCoefficients<K>>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const ReactionInput<K>& input) :
        equation(input.equation), reactants(input.reactants), products(input.products),
        forward(input.forward), backward_model(input.backward_model), backward(input.backward),
        third_body(input.third_body), efficiencies(input.efficiencies)
    {
    }
    ReactionInput<K> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .equation = equation.c_str(),
          .reactants = reactants.borrow(std::get<0>(scratch)),
          .products = products.borrow(std::get<1>(scratch)),
          .forward = forward.borrow(std::get<2>(scratch)),
          .backward_model = backward_model,
          .backward = backward.borrow(std::get<3>(scratch)),
          .third_body = third_body,
          .efficiencies = efficiencies.borrow(std::get<4>(scratch)),
        };
    }
};

template <std::size_t R, std::size_t K> struct InputStorage<MechanismInput<R, K>> {
    double temperature_lower {};
    double temperature_upper {};
    InputStorage<std::array<ReactionInput<K>, R>> reactions {};
    using Scratch = std::tuple<typename InputStorage<std::array<ReactionInput<K>, R>>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const MechanismInput<R, K>& input) :
        temperature_lower(input.temperature_lower), temperature_upper(input.temperature_upper),
        reactions(input.reactions)
    {
    }
    MechanismInput<R, K> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .temperature_lower = temperature_lower,
          .temperature_upper = temperature_upper,
          .reactions = reactions.borrow(std::get<0>(scratch)),
        };
    }
};

}  // namespace detail

// Builds from owned and direct inputs, keeping borrowed views inside this call.
template <typename... Explicit, typename... Inputs>
    requires detail::HasOwnedInput<Inputs...>
auto make_reaction_mechanism(const Inputs&... inputs)
{
    return detail::with_inputs(
      [](const auto&... borrowed) { return make_reaction_mechanism<Explicit...>(borrowed...); },
      inputs...);
}

}  // namespace canary

#endif  // CANARY_REACTIONS_HPP
