#ifndef CANARY_ENERGY_EXCHANGE_HPP
#define CANARY_ENERGY_EXCHANGE_HPP

// Exchange of vibrational-electronic energy in a two-temperature mixture:
// per cell, the source Q of vibrational-electronic energy per unit volume and
// time, Q = Q_VT + Q_CV; and the energy exchange between the modes of a
// multi-temperature mixture, with free electrons (see
// MultiTemperatureExchange).
//
// Species i has particle mass m_i, number density n_i = rho y_i / m_i and
// vibrational-electronic energy per unit mass e_ve,i (from the thermo
// model). T and T_ve are the translational and vibrational-electronic
// temperatures; neither is clamped here.
//
// Vibration-translation exchange (Landau-Teller), over the pairs (p, q) of a
// relaxing species p and a collision partner q:
//   Q_VT = sum_pairs rho y_p (e_ve,p(T) - e_ve,p(T_ve)) / tau_pq,
//   tau_pq = tau_MW,pq + tau_corr,pq.
// tau_MW is the correlation of Millikan and White (J. Chem. Phys. 39, 3209,
// 1963) with the bath pressure x_q p of the partner, x_q = n_q / sum_j n_j
// and p = sum_j n_j k_B T:
//   tau_MW = (p_atm / (x_q p)) exp(a (T^(-1/3) - b) - 18.42)
// in seconds, with T in K, p in Pa and p_atm = 101325 Pa. tau_corr is a
// high-temperature correction (RelaxationCorrection), zero for a pair
// without one:
//   tau_corr = 1 / (sigma c n_q),
//   c = sqrt(8 k_B T / (pi mu_pq)),  mu_pq = m_p m_q / (m_p + m_q).
// Both terms are proportional to 1 / n_q, so
//   1 / tau_pq = n_q / (exp(a' T^(-1/3) + c_pq) / (k_B T) + T^beta / K_pq)
// with a', c_pq, beta and K_pq constants of the pair (see
// VibrationTranslationPair); a partner with n_q = 0 contributes exactly 0.
//
// Chemistry-vibration coupling (impartial dissociation: Treanor and Marrone,
// Phys. Fluids 5, 1962; Knab, Fruhauf and Messerschmid, J. Thermophys. Heat
// Transf. 9(2), 1995, eqs. 36-39 with U -> infinity), over couplings of a
// reaction r and a molecule s with net coefficient change dnu = nu''_s -
// nu'_s and rates of progress q_f, q_b:
//   production = max(dnu, 0) q_f + max(-dnu, 0) q_b,
//   loss = max(-dnu, 0) q_f + max(dnu, 0) q_b,
//   Q_CV = sum_couplings production G(0) - loss G(1 / T_ve - 1 / T),
// where G is the mean vibrational energy per molecule of a harmonic
// oscillator truncated at the dissociation energy, with dissociation equally
// likely from every level:
//   G(x) = k_B (theta_v / (exp(theta_v x) - 1) - theta_D / (exp(theta_D x) - 1)),
//   G(0) = k_B (theta_D - theta_v) / 2,
// with the characteristic vibrational temperature theta_v and theta_D = D /
// R_u for the dissociation energy D per mole.

#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <numbers>
#include <string>
#include <string_view>
#include <utility>

namespace canary
{

// Host input of the Millikan-White constants of a pair: a in K^(1/3) (any
// finite value), b in K^(-1/3).
struct MillikanWhiteInput {
    double a;
    double b;
};

// High-temperature corrections of the vibration-translation relaxation
// time, tau_corr = 1 / (sigma c n_q), with the cross-section sigma:
// - none: no correction (tau_corr = 0);
// - park_1985: sigma = sigma_0 (C. Park, Nonequilibrium Hypersonic
//   Aerothermodynamics, 1990);
// - park_1993: sigma = sigma_0 (T_ref / T)^2 (C. Park, J. Thermophys. Heat
//   Transf. 7(3), 1993, p. 387);
// - kim_2021: sigma = sigma_0 T^omega, T in K (J. G. Kim and S. M. Jo, Int.
//   J. Heat Mass Transf. 169, 2021).
enum class RelaxationCorrection : std::uint8_t { none, park_1985, park_1993, kim_2021 };

// Host input of a high-temperature correction: the model, sigma_0 in m^2,
// T_ref in K (park_1993 only) and omega (kim_2021 only). Values a model
// does not use are ignored.
struct RelaxationCorrectionInput {
    RelaxationCorrection model;
    double cross_section = 0.0;
    double reference_temperature = 0.0;
    double exponent = 0.0;
};

// Host input of one vibration-translation pair: the relaxing species p, the
// collision partner q, the Millikan-White constants and the high-temperature
// correction.
struct VibrationTranslationInput {
    const char* species;
    const char* partner;
    MillikanWhiteInput millikan_white;
    RelaxationCorrectionInput correction;
};

// Host input of one chemistry-vibration coupling: the index of a reaction in
// the mechanism, the molecule whose vibrational energy the reaction changes,
// its dissociation energy D in J/mol and its characteristic vibrational
// temperature theta_v in K.
struct ChemistryVibrationInput {
    int reaction;
    const char* species;
    double dissociation_energy;
    double vibrational_temperature;
};

// Host input of the energy exchange of a mixture: V vibration-translation
// pairs and C chemistry-vibration couplings, either of which may be 0.
template <std::size_t V, std::size_t C> struct EnergyExchangeInput {
    std::array<VibrationTranslationInput, V> vibration_translation;
    std::array<ChemistryVibrationInput, C> chemistry_vibration;
};

// Cross-section of an elastic collision of an electron with a heavy
// partner:
// - gnoffo_neutral: a neutral partner, sigma = a + b T' + c T'^2 in m^2
//   with T' = min(T_e, maximum_temperature) in K (P. A. Gnoffo, R. N. Gupta
//   and J. L. Shinn, NASA TP-2867, 1989);
// - coulomb: a charged partner, sigma = (2 pi / 9) L^2 ln(Lambda) with
//   L = e^2 / (4 pi eps_0 k_B T_e), Lambda = 9 / (4 pi n_e' L^3) and
//   n_e' = max(n_e, 1 m^-3).
enum class ElectronCrossSection : std::uint8_t { gnoffo_neutral, coulomb };

// Host input of one elastic electron-heavy pair: the partner, the
// cross-section model and, for gnoffo_neutral, a in m^2, b in m^2/K, c in
// m^2/K^2 and the maximum temperature in K (ignored for coulomb).
struct ElectronTranslationInput {
    const char* partner;
    ElectronCrossSection cross_section;
    double a;
    double b;
    double c;
    double maximum_temperature;
};

// Host input of Candler's vibration-electron relaxation time:
// log10(p_e tau / (1 atm s)) = f(x) with x = log10(T_e / 1 K),
// f_low = a_low x^2 + b_low x + c_low for T_e <= blend_lower,
// f_high = a_high x^2 + b_high x + c_high for T_e >= blend_upper, and
// w f_low + (1 - w) f_high between, with w = cos^2(pi / 2 (T_e -
// blend_lower) / (blend_upper - blend_lower)); blend bounds in K.
struct CandlerFitInput {
    double a_low;
    double b_low;
    double c_low;
    double a_high;
    double b_high;
    double c_high;
    double blend_lower;
    double blend_upper;
};

// Host input of one vibration-electron exchange: the species whose
// vibration relaxes with the free electrons, and its relaxation time.
struct VibrationElectronInput {
    const char* species;
    CandlerFitInput relaxation_time;
};

// Host input of one electron-impact ionisation: the index of the reaction
// in the mechanism, the neutral species it ionises (to the species named
// with "+" appended) and its ionisation energy in J/mol.
struct ElectronImpactIonisationInput {
    int reaction;
    const char* species;
    double ionisation_energy;
};

// Host input of the exchanges with the free electrons: E elastic
// electron-heavy pairs, VE vibration-electron exchanges and I
// electron-impact ionisations, any of which may be 0.
template <std::size_t E, std::size_t VE, std::size_t I> struct ElectronExchangeInput {
    std::array<ElectronTranslationInput, E> electron_translation;
    std::array<VibrationElectronInput, VE> vibration_electron;
    std::array<ElectronImpactIonisationInput, I> electron_impact_ionisation;
};

// One vibration-translation pair in model units: the indices of the
// relaxing species and of the partner in the species order of the thermo
// model, and the constants of
//   1 / tau = n_q / (exp(a' T^(-1/3) + c) / (k_B T) + T^beta / K)
// with T, k_B and n_q in model units and tau in the time unit of U, from
// tau_MW and tau_corr n_q = T^beta / K. With tau_T the temperature scale,
// t_ref the time scale and n_ref the number density scale of U (all 1 in
// SI):
// - millikan_white_a: a' = a tau_T^(-1/3);
// - millikan_white_offset: c = ln(p_atm (1 s) in model units) - a b -
//   18.42;
// - correction_exponent: beta = -1/2 (park_1985), 3/2 (park_1993),
//   -1/2 - omega (kim_2021), 0 without a correction;
// - correction_coefficient: 1 / K, where K = K_SI n_ref t_ref / tau_T^beta
//   and K_SI = sigma_0 sqrt(8 k_B / (pi mu_pq)) (times T_ref^2 for
//   park_1993) in SI; 0 without a correction.
template <std::floating_point P, UnitSystem U> struct VibrationTranslationPair {
    std::size_t species;
    std::size_t partner;
    P millikan_white_a;
    P millikan_white_offset;
    P correction_exponent;
    P correction_coefficient;
};

// One chemistry-vibration coupling in model units: the reaction index, the
// species index, theta_D = D / R_u and theta_v.
template <std::floating_point P, UnitSystem U> struct ChemistryVibrationCoupling {
    std::size_t reaction;
    std::size_t species;
    Temperature<P, U> dissociation_temperature;
    Temperature<P, U> vibrational_temperature;
};

// Energy exchange over the N species of a thermo model, in unit system U,
// as derived by make_energy_exchange: V vibration-translation pairs, C
// chemistry-vibration couplings, the particle mass of every species and the
// Boltzmann constant.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t V, std::size_t C>
    requires(N > 0)
struct EnergyExchange {
    using parameter_type = P;
    using units = U;

    std::array<VibrationTranslationPair<P, U>, V> vibration_translation;
    std::array<ChemistryVibrationCoupling<P, U>, C> chemistry_vibration;
    std::array<ParticleMass<P, U>, N> mass;
    BoltzmannConstant<P, U> boltzmann;
};

namespace detail
{

// Index of the named species, checked: throws ParameterError when name is
// null or unknown.
template <std::size_t N>
std::size_t known_species(const std::array<std::string_view, N>& names,
                          const char* name,
                          std::string_view role)
{
    if (name == nullptr) { throw ParameterError(std::format("{} has no name", role)); }
    const std::size_t index = species_index(names, name);
    if (index == N) { throw ParameterError(std::format("unknown {} \"{}\"", role, name)); }
    return index;
}

// A finite, positive value converted to P; throws ParameterError otherwise,
// naming what.
template <std::floating_point P> P positive_parameter(double value, std::string_view what)
{
    const auto converted = static_cast<P>(value);
    if (!std::isfinite(value) || !(value > 0.0) || !std::isfinite(converted) ||
        !(converted > P(0))) {
        throw ParameterError(std::format("{} must be finite and positive, got {}", what, value));
    }
    return converted;
}

template <std::floating_point P, UnitSystem U, std::size_t N>
VibrationTranslationPair<P, U>
  make_vibration_translation_pair(const VibrationTranslationInput& input,
                                  const std::array<std::string_view, N>& names,
                                  const SpeciesData<P, U, N>& species,
                                  const UnitScales<P, U>& scales)
{
    VibrationTranslationPair<P, U> pair {};
    pair.species = known_species(names, input.species, "species");
    pair.partner = known_species(names, input.partner, "partner");
    const double a = input.millikan_white.a;
    const double b = input.millikan_white.b;
    if (!std::isfinite(a)) {
        throw ParameterError(std::format("Millikan-White a must be finite, got {}", a));
    }
    if (!std::isfinite(b)) {
        throw ParameterError(std::format("Millikan-White b must be finite, got {}", b));
    }
    const UnitScales<double, U> s = double_scales(scales);
    const double tau = scale_factor<kind::temperature>(s);
    // p_atm times one second, a pressure times a time.
    const Dimension pressure_time {.mass = 1, .length = -1, .time = -1};
    const double atmosphere_second = 101325.0 / scale_factor(pressure_time, s);
    pair.millikan_white_a = static_cast<P>(a / std::cbrt(tau));
    pair.millikan_white_offset = static_cast<P>(std::log(atmosphere_second) - (a * b) - 18.42);
    if (!std::isfinite(pair.millikan_white_a) || !std::isfinite(pair.millikan_white_offset)) {
        throw ParameterError(
          std::format("Millikan-White a = {}, b = {} are not representable in the model's units",
                      a,
                      b));
    }
    const RelaxationCorrectionInput& correction = input.correction;
    // K_SI / sqrt(8 k_B / (pi mu)) and beta.
    double k_factor = 0.0;
    double beta = 0.0;
    switch (correction.model) {
    case RelaxationCorrection::none: break;
    case RelaxationCorrection::park_1985:
        k_factor = positive_parameter<double>(correction.cross_section, "correction cross-section");
        beta = -0.5;
        break;
    case RelaxationCorrection::park_1993: {
        const auto sigma =
          positive_parameter<double>(correction.cross_section, "correction cross-section");
        const auto t_ref = positive_parameter<double>(correction.reference_temperature,
                                                      "correction reference temperature");
        k_factor = sigma * t_ref * t_ref;
        beta = 1.5;
        break;
    }
    case RelaxationCorrection::kim_2021:
        k_factor = positive_parameter<double>(correction.cross_section, "correction cross-section");
        if (!std::isfinite(correction.exponent)) {
            throw ParameterError(
              std::format("correction exponent must be finite, got {}", correction.exponent));
        }
        beta = -0.5 - correction.exponent;
        break;
    default: throw ParameterError("unknown relaxation-time correction model");
    }
    pair.correction_exponent = static_cast<P>(beta);
    pair.correction_coefficient = P(0);
    if (correction.model != RelaxationCorrection::none) {
        const double mass_scale = scale_factor<kind::particle_mass>(s);
        const double m_p = static_cast<double>(species.mass[pair.species].value()) * mass_scale;
        const double m_q = static_cast<double>(species.mass[pair.partner].value()) * mass_scale;
        const double reduced_mass = m_p * m_q / (m_p + m_q);
        const double k_si = k_factor * std::sqrt(8.0 * constants::boltzmann<double> /
                                                 (std::numbers::pi * reduced_mass));
        // Number density times time.
        const Dimension density_time {.length = -3, .time = 1, .count = 1};
        const double k_model = k_si * scale_factor(density_time, s) / std::pow(tau, beta);
        pair.correction_coefficient =
          positive_parameter<P>(1.0 / k_model, "correction coefficient 1 / K");
    }
    return pair;
}

template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
ChemistryVibrationCoupling<P, U>
  make_chemistry_vibration_coupling(const ChemistryVibrationInput& input,
                                    const std::array<std::string_view, N>& names,
                                    const ReactionMechanism<P, U, N, R>& mechanism,
                                    const UnitScales<P, U>& scales)
{
    ChemistryVibrationCoupling<P, U> coupling {};
    if (input.reaction < 0 || static_cast<std::size_t>(input.reaction) >= R) {
        throw ParameterError(
          std::format("reaction index must be between 0 and {}, got {}", R - 1, input.reaction));
    }
    coupling.reaction = static_cast<std::size_t>(input.reaction);
    coupling.species = known_species(names, input.species, "species");
    const Reaction<P, U, N>& r = mechanism.reactions[coupling.reaction];
    if (r.products[coupling.species] == r.reactants[coupling.species]) {
        throw ParameterError(std::format("species \"{}\" has no net change in reaction {}",
                                         input.species,
                                         input.reaction));
    }
    const double d = input.dissociation_energy;
    const double theta_v = input.vibrational_temperature;
    if (!std::isfinite(d) || !(d > 0.0)) {
        throw ParameterError(
          std::format("dissociation energy must be finite and positive, got {}", d));
    }
    if (!std::isfinite(theta_v) || !(theta_v > 0.0)) {
        throw ParameterError(
          std::format("vibrational temperature must be finite and positive, got {}", theta_v));
    }
    const double tau = scale_factor<kind::temperature>(double_scales(scales));
    const double theta_d = d / constants::universal_gas_constant<double>;
    if (!(theta_d > theta_v)) {
        throw ParameterError(
          std::format("dissociation temperature must exceed vibrational temperature, got {} and {}",
                      theta_d,
                      theta_v));
    }
    coupling.dissociation_temperature =
      Temperature<P, U> {positive_parameter<P>(theta_d / tau, "dissociation temperature")};
    coupling.vibrational_temperature =
      Temperature<P, U> {positive_parameter<P>(theta_v / tau, "vibrational temperature")};
    if (!(coupling.dissociation_temperature.value() > coupling.vibrational_temperature.value())) {
        throw ParameterError(std::format("stored dissociation temperature must exceed stored "
                                         "vibrational temperature, got {} and {}",
                                         coupling.dissociation_temperature.value(),
                                         coupling.vibrational_temperature.value()));
    }
    return coupling;
}

// Throws ParameterError when species.n is not N or a name is listed twice.
template <std::floating_point P, UnitSystem U, std::size_t N>
void check_species_names(const std::array<std::string_view, N>& names,
                         const SpeciesData<P, U, N>& species)
{
    if (species.n != static_cast<int>(N)) {
        throw ParameterError(
          std::format("species data must have {} species, got {}", N, species.n));
    }
    for (std::size_t i = 0; i < N; ++i) {
        if (species_index(names, names[i]) != i) {
            throw ParameterError(std::format("species \"{}\" is named twice", names[i]));
        }
    }
}

inline const char* name_or_empty(const char* name) { return name != nullptr ? name : ""; }

// The vibration-translation pairs and chemistry-vibration couplings of input
// into vibration_translation and chemistry_vibration. check_species(i) is
// called with the index of each relaxing species and of each coupling's
// species, and may throw ParameterError. Errors name the entry.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t V,
          std::size_t C,
          typename Check>
void make_vibrational_exchanges(
  const EnergyExchangeInput<V, C>& input,
  const std::array<std::string_view, N>& names,
  const SpeciesData<P, U, N>& species,
  const ReactionMechanism<P, U, N, R>& mechanism,
  const UnitScales<P, U>& scales,
  const Check& check_species,
  std::array<VibrationTranslationPair<P, U>, V>& vibration_translation,
  std::array<ChemistryVibrationCoupling<P, U>, C>& chemistry_vibration)
{
    for (std::size_t e = 0; e < V; ++e) {
        const VibrationTranslationInput& pair = input.vibration_translation[e];
        try {
            vibration_translation[e] =
              make_vibration_translation_pair(pair, names, species, scales);
            check_species(vibration_translation[e].species);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("vibration-translation pair {} ({}, {}): {}",
                                             e,
                                             name_or_empty(pair.species),
                                             name_or_empty(pair.partner),
                                             error.what()));
        }
    }
    for (std::size_t c = 0; c < C; ++c) {
        const ChemistryVibrationInput& coupling = input.chemistry_vibration[c];
        try {
            chemistry_vibration[c] =
              make_chemistry_vibration_coupling(coupling, names, mechanism, scales);
            check_species(chemistry_vibration[c].species);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("chemistry-vibration coupling {} ({}): {}",
                                             c,
                                             name_or_empty(coupling.species),
                                             error.what()));
        }
    }
}

}  // namespace detail

// Energy exchange of input over the N species of a thermo model, named in
// its species order, with their particle masses from species, and the
// reactions of mechanism (over the same species), in the unit system of the
// scales. Throws ParameterError when species.n is not N or a name is listed
// twice, and, with "vibration-translation pair <index> (<species>,
// <partner>)" or "chemistry-vibration coupling <index> (<species>)" in the
// message, when:
// - a species or partner name is null or unknown;
// - a Millikan-White a or b is not finite;
// - a correction is used with a cross-section, a reference temperature
//   (park_1993) that is not finite and positive, or an exponent (kim_2021)
//   that is not finite, or its model is unknown;
// - a reaction index is outside [0, R), or the species has no net change of
//   coefficient in the reaction;
// - a dissociation energy or vibrational temperature is not finite and
//   positive, or dissociation temperature does not exceed vibrational
//   temperature before and after conversion to the stored scalar and units;
// - a derived constant is not finite (and, for the correction coefficient
//   and the coupling temperatures, positive) in P.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t V,
          std::size_t C>
    requires(N > 0 && R > 0)
EnergyExchange<P, U, N, V, C> make_energy_exchange(const EnergyExchangeInput<V, C>& input,
                                                   const std::array<std::string_view, N>& names,
                                                   const SpeciesData<P, U, N>& species,
                                                   const ReactionMechanism<P, U, N, R>& mechanism,
                                                   const UnitScales<P, U>& scales)
{
    detail::check_species_names(names, species);
    EnergyExchange<P, U, N, V, C> exchange {};
    detail::make_vibrational_exchanges(
      input,
      names,
      species,
      mechanism,
      scales,
      [](std::size_t /*i*/) {},
      exchange.vibration_translation,
      exchange.chemistry_vibration);
    exchange.mass = species.mass;
    exchange.boltzmann = boltzmann_constant(scales);
    return exchange;
}

// A thermo model of N species whose
// species_vibrational_electronic_energy(gas, temp) gives e_ve,i of every
// species at a temperature with value type T.
template <typename G, typename T, std::size_t N>
concept HasVibrationalElectronicEnergy =
  GasModelType<G> && Scalar<T> && requires(const G& gas, Temperature<T, typename G::units> temp) {
      {
          species_vibrational_electronic_energy(gas, temp)
      } -> std::convertible_to<std::array<SpecificEnergy<T, typename G::units>, N>>;
  };

// A thermo model and value type for which vibration_translation_source of
// an exchange over N species with parameter type P can be evaluated.
template <typename G, typename T, typename P, std::size_t N>
concept RelaxingModel = HasVibrationalElectronicEnergy<G, T, N> && std::floating_point<P> &&
                        ScalarWith<T, P> && HasExp<T> && HasLog<T> && HasSqrt<T>;

namespace detail
{

// A value and its derivative.
template <std::floating_point V> struct ValueAndDerivative {
    V value;
    V derivative;
};

// h(y) = 1 / (exp(y) - 1) - 1 / y + 1 / 2 and dh/dy, from the series
// h(y) = sum_k B_2k y^(2k - 1) / (2k)! with the Bernoulli numbers B_2k,
// k = 1..9. For |y| <= 1/2 the first omitted term is below 1e-20 of the
// first term.
template <std::floating_point V>
CANARY_HOST_DEVICE constexpr ValueAndDerivative<V> bernoulli_remainder(V y)
{
    // B_2k / (2k)!.
    constexpr std::array<V, 9> c {
      V(1.0 / 12.0),
      V(-1.0 / 720.0),
      V(1.0 / 30240.0),
      V(-1.0 / 1209600.0),
      V(1.0 / 47900160.0),
      V(-691.0 / 1307674368000.0),
      V(1.0 / 74724249600.0),
      V(-3617.0 / (510.0 * 20922789888000.0)),
      V(43867.0 / (798.0 * 6402373705728000.0)),
    };
    const V z = y * y;
    V h = c[8];
    V dh = V(17) * c[8];
    for (std::size_t k = 8; k-- > 0;) {
        h = c[k] + (z * h);
        dh = (V((2 * k) + 1) * c[k]) + (z * dh);
    }
    return {.value = y * h, .derivative = dh};
}

// g(x) = G(x) / k_B = theta_v / (exp(theta_v x) - 1) - theta_d /
// (exp(theta_d x) - 1) and dg/dx, for every real x, without overflow:
// - |x| max(theta_v, theta_d) <= 1/2: g = (theta_d - theta_v) / 2 +
//   theta_v h(theta_v x) - theta_d h(theta_d x), with h of
//   bernoulli_remainder, in which the divergent terms 1 / x cancel exactly
//   (g(0) = (theta_d - theta_v) / 2);
// - otherwise, at |x|, theta / (exp(theta |x|) - 1) = theta exp(-theta |x|)
//   / (-expm1(-theta |x|)), and for x < 0, g(x) = (theta_d - theta_v) -
//   g(|x|), so that dg/dx is even.
// The relative error is a few ulp when theta_d is not close to theta_v.
template <std::floating_point V>
CANARY_HOST_DEVICE ValueAndDerivative<V> truncated_oscillator_energy(V theta_v, V theta_d, V x)
{
    const V half_difference = (theta_d - theta_v) / V(2);
    const V abs_x = x < V(0) ? -x : x;
    const V theta_max = theta_v < theta_d ? theta_d : theta_v;
    if (theta_max * abs_x <= V(0.5)) {
        const ValueAndDerivative<V> h_v = bernoulli_remainder(theta_v * x);
        const ValueAndDerivative<V> h_d = bernoulli_remainder(theta_d * x);
        return {
          .value = half_difference + (theta_v * h_v.value) - (theta_d * h_d.value),
          .derivative = (theta_v * theta_v * h_v.derivative) - (theta_d * theta_d * h_d.derivative),
        };
    }
    const V e_v = std::exp(-theta_v * abs_x);
    const V d_v = -std::expm1(-theta_v * abs_x);
    const V e_d = std::exp(-theta_d * abs_x);
    const V d_d = -std::expm1(-theta_d * abs_x);
    const V g = (theta_v * e_v / d_v) - (theta_d * e_d / d_d);
    const V dg = (theta_d * theta_d * e_d / (d_d * d_d)) - (theta_v * theta_v * e_v / (d_v * d_v));
    return {.value = x < V(0) ? (theta_d - theta_v) - g : g, .derivative = dg};
}

// sum over the partners q of each species i of 1 / tau_iq at translational
// temperature t, with the number densities n of every species (see
// VibrationTranslationPair). The bath pressure of partner q in tau_MW is
// n_q k_B t_bath, with t_bath = p / (k_B sum_j n_j) for the mixture
// pressure p (t_bath = t without free electrons).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t V>
    requires ScalarWith<T, P> && HasExp<T> && HasLog<T> && HasSqrt<T>
CANARY_HOST_DEVICE std::array<T, N>
  inverse_relaxation_times(const std::array<VibrationTranslationPair<P, U>, V>& pairs,
                           BoltzmannConstant<P, U> boltzmann,
                           const std::array<T, N>& n,
                           T t,
                           T t_bath)
{
    using Real = ScalarTraits<T>::value_type;
    const T ln_t = math::log(t);
    const T t_inverse_cbrt = math::exp(T(-ln_t / T(Real(3))));
    const T sqrt_t = math::sqrt(t);
    const T kt = T(boltzmann.value() * t_bath);
    std::array<T, N> inverse_tau {};
    inverse_tau.fill(T(Real(0)));
    for (const VibrationTranslationPair<P, U>& pair : pairs) {
        const T millikan_white =
          T(math::exp(T(T(pair.millikan_white_a * t_inverse_cbrt) + pair.millikan_white_offset)) /
            kt);
        T correction = T(Real(0));
        if (pair.correction_coefficient != P(0)) {
            // T^beta, with the exponents of park_1993 and park_1985 from sqrt(T).
            T power = T(Real(0));
            if (pair.correction_exponent == P(1.5)) {
                power = T(t * sqrt_t);
            } else if (pair.correction_exponent == P(-0.5)) {
                power = T(T(Real(1)) / sqrt_t);
            } else {
                power = math::exp(T(pair.correction_exponent * ln_t));
            }
            correction = T(pair.correction_coefficient * power);
        }
        inverse_tau[pair.species] =
          inverse_tau[pair.species] + T(n[pair.partner] / T(millikan_white + correction));
    }
    return inverse_tau;
}

// production G(0) - loss G(1 / t_m - 1 / t), over k_B, of one coupling
// with the rates of progress of mechanism, at translational temperature t
// and the temperature t_m of the molecule's vibration. G(1 / t_m - 1 / t)
// is evaluated on value parts by truncated_oscillator_energy; a dual number
// receives the first derivative dG/dx times the derivative of x. With t_m
// equal to t, G(1 / t_m - 1 / t) equals G(0) exactly.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T coupling_energy_rate(const ChemistryVibrationCoupling<P, U>& coupling,
                                          const ReactionMechanism<P, U, N, R>& mechanism,
                                          const RatesOfProgress<T, U, R>& rates,
                                          T t,
                                          T t_m)
{
    using Real = ScalarTraits<T>::value_type;
    // x = 1 / t_m - 1 / t, its value part, and x minus its value part: zero
    // for a floating-point type, the derivative part for a dual number.
    const T x = T(t - t_m) / T(t * t_m);
    const Real x_value = ScalarTraits<T>::value(x);
    const T x_derivative = x - T(x_value);
    // NOLINTNEXTLINE(misc-static-assert): run-time reaction index
    CANARY_ASSERT(coupling.reaction < R);
    const Reaction<P, U, N>& r = mechanism.reactions[coupling.reaction];
    const int dnu = r.products[coupling.species] - r.reactants[coupling.species];
    const auto made = static_cast<Real>(dnu > 0 ? dnu : 0);
    const auto destroyed = static_cast<Real>(dnu < 0 ? -dnu : 0);
    const T q_f = rates.forward[coupling.reaction].value();
    const T q_b = rates.backward[coupling.reaction].value();
    const T production = T(made * q_f) + T(destroyed * q_b);
    const T loss = T(destroyed * q_f) + T(made * q_b);
    const auto theta_v = static_cast<Real>(coupling.vibrational_temperature.value());
    const auto theta_d = static_cast<Real>(coupling.dissociation_temperature.value());
    const ValueAndDerivative<Real> g = truncated_oscillator_energy(theta_v, theta_d, x_value);
    const T g_vanish = T(g.value) + (x_derivative * T(g.derivative));
    const T g_appear = T((theta_d - theta_v) / Real(2));
    return T(production * g_appear) - T(loss * g_vanish);
}

// rho y_i and n_i = rho y_i / m_i of every species.
template <Scalar T, std::size_t N> struct PartialDensities {
    std::array<T, N> mass;
    std::array<T, N> number;
};

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE PartialDensities<T, N>
  partial_densities(const std::array<ParticleMass<P, U>, N>& mass,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho)
{
    PartialDensities<T, N> out {};
    for (std::size_t i = 0; i < N; ++i) {
        out.mass[i] = T(rho.value() * mix.fractions[i].value());
        out.number[i] = T(out.mass[i] / mass[i].value());
    }
    return out;
}

}  // namespace detail

// Q_VT in a cell of density rho, translational temperature temp and
// vibrational-electronic temperature temp_ve, with the normalised mass
// fractions of mix (see the top of this file). The gas model is called
// twice, through species_vibrational_electronic_energy at temp and at
// temp_ve. The status of mix is not checked.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t V,
          std::size_t C,
          typename G>
    requires RelaxingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE PowerDensity<T, U>
  vibration_translation_source(const EnergyExchange<P, U, N, V, C>& exchange,
                               const G& gas,
                               const NormalisedComposition<T, U, N>& mix,
                               Density<T, U> rho,
                               Temperature<T, U> temp,
                               Temperature<T, U> temp_ve)
{
    using Real = ScalarTraits<T>::value_type;
    const detail::PartialDensities<T, N> d = detail::partial_densities(exchange.mass, mix, rho);
    const std::array<T, N> inverse_tau =
      detail::inverse_relaxation_times(exchange.vibration_translation,
                                       exchange.boltzmann,
                                       d.number,
                                       temp.value(),
                                       temp.value());
    const std::array<SpecificEnergy<T, U>, N> e_t =
      species_vibrational_electronic_energy(gas, temp);
    const std::array<SpecificEnergy<T, U>, N> e_ve =
      species_vibrational_electronic_energy(gas, temp_ve);
    T q = T(Real(0));
    for (std::size_t i = 0; i < N; ++i) {
        q = q + T(d.mass[i] * T(e_t[i].value() - e_ve[i].value()) * inverse_tau[i]);
    }
    return PowerDensity<T, U> {q};
}

// Q_CV of the couplings of exchange, with rates the rates of progress of
// mechanism in the cell, at translational temperature temp and
// vibrational-electronic temperature temp_ve (see the top of this file and
// detail::coupling_energy_rate).
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t V,
          std::size_t C>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE PowerDensity<T, U>
  chemistry_vibration_source(const EnergyExchange<P, U, N, V, C>& exchange,
                             const ReactionMechanism<P, U, N, R>& mechanism,
                             const RatesOfProgress<T, U, R>& rates,
                             Temperature<T, U> temp,
                             Temperature<T, U> temp_ve)
{
    using Real = ScalarTraits<T>::value_type;
    T q = T(Real(0));
    for (const ChemistryVibrationCoupling<P, U>& coupling : exchange.chemistry_vibration) {
        q = q +
            detail::coupling_energy_rate(coupling, mechanism, rates, temp.value(), temp_ve.value());
    }
    return PowerDensity<T, U> {T(exchange.boltzmann.value() * q)};
}

// Species mass production rates and source of vibrational-electronic energy
// of a cell of a two-temperature mixture.
template <Scalar T, UnitSystem U, std::size_t N>
    requires(N > 0)
struct TwoTemperatureSources {
    std::array<MassProductionRate<T, U>, N> species;
    PowerDensity<T, U> vibrational_electronic;
};

// production_rates of mechanism and Q_VT + Q_CV of exchange in a cell of
// density rho, translational temperature temp and vibrational-electronic
// temperature temp_ve, with the normalised mass fractions of mix. The rates
// of progress are evaluated once (rates_of_progress with temp_ve) and used
// for both.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t V,
          std::size_t C,
          typename G>
    requires ReactingModel<G, T, P, N> && RelaxingModel<G, T, P, N> &&
             std::same_as<typename G::units, U>
CANARY_HOST_DEVICE TwoTemperatureSources<T, U, N>
  two_temperature_sources(const ReactionMechanism<P, U, N, R>& mechanism,
                          const EnergyExchange<P, U, N, V, C>& exchange,
                          const G& gas,
                          const NormalisedComposition<T, U, N>& mix,
                          Density<T, U> rho,
                          Temperature<T, U> temp,
                          Temperature<T, U> temp_ve)
{
    const RatesOfProgress<T, U, R> rates =
      rates_of_progress(mechanism, gas, mix, rho, temp, temp_ve);
    const PowerDensity<T, U> q_vt =
      vibration_translation_source(exchange, gas, mix, rho, temp, temp_ve);
    const PowerDensity<T, U> q_cv =
      chemistry_vibration_source(exchange, mechanism, rates, temp, temp_ve);
    return {
      .species = detail::mass_production_rates(mechanism, rates),
      .vibrational_electronic = PowerDensity<T, U> {T(q_vt.value() + q_cv.value())},
    };
}

// Energy exchange of a multi-temperature mixture with M energy modes.
//
// The source of each mode is in energy per unit volume and time. Species i
// has the vibration mode m_i (the mode holding its vibration) and
// e_i,m(T') is the energy per unit mass of its components in mode m at the
// temperature T' (species_mode_energy). T is the translational temperature,
// T_m the temperature of mode m, T_e that of the electron mode, n_e and m_e
// the number density and particle mass of e-, and p_e = n_e k_B T_e. The
// exchanges are:
// - vibration-translation, into mode m_p:
//   Q = rho y_p (e_p,m_p(T) - e_p,m_p(T_m_p)) sum_q 1 / tau_pq,
//   with tau_pq as for EnergyExchange, except that the bath pressure x_q p
//   of tau_MW has the mixture pressure p = sum over heavy j of n_j k_B T +
//   p_e and x_q = n_q / sum_j n_j over all species, e- included;
// - chemistry-vibration, into mode m_s: Q_CV of EnergyExchange with T_ve
//   replaced by T_m_s;
// - elastic electron-heavy (electron-translation; Gnoffo, Gupta and Shinn,
//   NASA TP-2867, 1989; Glass and Liu, J. Fluid Mech. 85, 1978), into the
//   electron mode:
//   Q_ET = 3 n_e k_B (T - T_e) sum_q (m_e / m_q) n_q sigma_eq c_e,
//   c_e = sqrt(8 k_B T_e / (pi m_e)), sigma_eq of ElectronCrossSection;
// - vibration-electron (Landau-Teller with Candler's relaxation time), of
//   species p: Q = rho y_p (e_p,m_p(T_e) - e_p,m_p(T_m_p)) / tau_pe, added
//   to mode m_p and subtracted from the electron mode, with
//   1 / tau_pe = (p_e / 1 atm) 10^(-f) / (1 s) (see CandlerFitInput), so
//   that p_e = 0 gives no exchange;
// - electron-impact ionisation (Gnoffo, Gupta and Shinn 1989), into the
//   electron mode: Q_CE = -sum (q_f,r - q_b,r) I_s / N_A over the
//   ionisations, with the rates of progress of reaction r and the
//   ionisation energy I_s per mole of the ionised species s.
// No temperature is clamped.

// One elastic electron-heavy pair in model units: the index of the partner,
// the cross-section model and, for gnoffo_neutral, the coefficients of
// sigma = a + b T' + c T'^2 with T' = min(T_e, maximum_temperature), sigma
// an area per particle (length^2 / count of U) and T' in the temperature
// unit of U.
template <std::floating_point P, UnitSystem U> struct ElectronTranslationPair {
    std::size_t partner;
    ElectronCrossSection cross_section;
    P a;
    P b;
    P c;
    Temperature<P, U> maximum_temperature;
};

// Candler's fit in model units: 1 / tau = p_e 10^(-f(x)) with
// x = log10(T_e), p_e, T_e and tau in model units, f_low = a_low x^2 +
// b_low x + c_low and f_high likewise, blended between the blend bounds as
// in CandlerFitInput. With s = log10(tau_T) for the temperature scale tau_T
// and K = t_ref p_ref / (1 atm (1 s)) for the time scale t_ref and the
// pressure scale p_ref (in SI, s = 0 and K = 1 / 101325), each fit
// a x^2 + b x + c of the input becomes
// a x^2 + (2 a s + b) x + (a s^2 + b s + c - log10 K).
template <std::floating_point P, UnitSystem U> struct CandlerFit {
    P a_low;
    P b_low;
    P c_low;
    P a_high;
    P b_high;
    P c_high;
    Temperature<P, U> blend_lower;
    Temperature<P, U> blend_upper;
};

// One vibration-electron exchange in model units: the index of the species
// and its relaxation time.
template <std::floating_point P, UnitSystem U> struct VibrationElectronCoupling {
    std::size_t species;
    CandlerFit<P, U> relaxation_time;
};

// One electron-impact ionisation in model units: the reaction index, the
// index of the ionised species and its ionisation energy per particle,
// I / N_A, in the energy unit of U.
template <std::floating_point P, UnitSystem U> struct ElectronImpactIonisationCoupling {
    std::size_t reaction;
    std::size_t species;
    P ionisation_energy;
};

// Energy exchange over the N species and M energy modes of a thermo model,
// in unit system U, as derived by make_multi_temperature_exchange: V
// vibration-translation pairs, C chemistry-vibration couplings, E elastic
// electron-heavy pairs, VE vibration-electron exchanges and I
// electron-impact ionisations; the vibration mode of every species (-1
// without one); the particle mass of every species; the Boltzmann constant;
// the index of e- and the mode that holds its energy (both -1 without
// free electrons); the constants of the Coulomb cross-section
// sigma = (coulomb_cross_section / T_e^2) ln(coulomb_logarithm_factor T_e^3
// / n_e'), that is (2 pi / 9) C^2 and 9 / (4 pi C^3) with
// C = L T_e = e^2 / (4 pi eps_0 k_B), in model units (a cross-section is an
// area per particle, length^2 / count, and Lambda counts particles); and
// minimum_electron_density, 1 m^-3 in the number density unit of U.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
    requires(N > 0 && M > 0)
struct MultiTemperatureExchange {
    using parameter_type = P;
    using units = U;
    static constexpr std::size_t mode_count = M;

    std::array<VibrationTranslationPair<P, U>, V> vibration_translation;
    std::array<ChemistryVibrationCoupling<P, U>, C> chemistry_vibration;
    std::array<ElectronTranslationPair<P, U>, E> electron_translation;
    std::array<VibrationElectronCoupling<P, U>, VE> vibration_electron;
    std::array<ElectronImpactIonisationCoupling<P, U>, I> electron_impact_ionisation;
    std::array<int, N> vibration_mode;
    std::array<ParticleMass<P, U>, N> mass;
    BoltzmannConstant<P, U> boltzmann;
    int electron;
    int electron_mode;
    P coulomb_cross_section;
    P coulomb_logarithm_factor;
    P minimum_electron_density;
};

namespace detail
{

// A finite value converted to P; throws ParameterError otherwise, naming
// what.
template <std::floating_point P> P finite_parameter(double value, std::string_view what)
{
    const auto converted = static_cast<P>(value);
    if (!std::isfinite(value) || !std::isfinite(converted)) {
        throw ParameterError(std::format("{} must be finite, got {}", what, value));
    }
    return converted;
}

// A quadratic on [0, maximum_temperature] has its minimum at an endpoint
// or at its vertex when convex. Widen stored values for the host-side check.
template <std::floating_point P>
void check_neutral_cross_section(P a, P b, P c, P maximum_temperature)
{
    const long double wa = a;
    const long double wb = b;
    const long double wc = c;
    const long double upper = maximum_temperature;
    const auto check = [&](long double temperature) {
        const long double sigma = std::fma(std::fma(wc, temperature, wb), temperature, wa);
        if (!std::isfinite(sigma) || sigma < 0.0L) {
            throw ParameterError(
              std::format("neutral cross-section must be finite and non-negative over its "
                          "temperature interval, got {} at {}",
                          sigma,
                          temperature));
        }
    };
    check(0.0L);
    check(upper);
    if (wc > 0.0L && wb < 0.0L) {
        const long double vertex = (-0.5L * wb) / wc;
        if (vertex > 0.0L && vertex < upper) { check(vertex); }
    }
}

template <std::floating_point P, UnitSystem U, std::size_t N>
ElectronTranslationPair<P, U>
  make_electron_translation_pair(const ElectronTranslationInput& input,
                                 const std::array<std::string_view, N>& names,
                                 const SpeciesData<P, U, N>& species,
                                 const UnitScales<P, U>& scales)
{
    ElectronTranslationPair<P, U> pair {};
    pair.partner = known_species(names, input.partner, "partner");
    if (names[pair.partner] == "e-") {
        throw ParameterError("the partner must be a heavy species, not e-");
    }
    const bool charged = species.charge[pair.partner].value() != P(0);
    pair.cross_section = input.cross_section;
    switch (input.cross_section) {
    case ElectronCrossSection::coulomb:
        if (!charged) { throw ParameterError("a coulomb partner must be charged"); }
        return pair;
    case ElectronCrossSection::gnoffo_neutral: break;
    default: throw ParameterError("unknown electron cross-section model");
    }
    if (charged) { throw ParameterError("a gnoffo_neutral partner must be neutral"); }
    const UnitScales<double, U> s = double_scales(scales);
    const double tau = scale_factor<kind::temperature>(s);
    // An area per particle.
    const double area = scale_factor(Dimension {.length = 2, .count = -1}, s);
    const auto a = finite_parameter<double>(input.a, "cross-section coefficient a");
    const auto b = finite_parameter<double>(input.b, "cross-section coefficient b");
    const auto c = finite_parameter<double>(input.c, "cross-section coefficient c");
    const auto t_max =
      positive_parameter<double>(input.maximum_temperature, "cross-section maximum temperature");
    check_neutral_cross_section(a, b, c, t_max);
    pair.a = finite_parameter<P>(a / area, "cross-section coefficient a in the model's units");
    pair.b =
      finite_parameter<P>(b * tau / area, "cross-section coefficient b in the model's units");
    pair.c =
      finite_parameter<P>(c * tau * tau / area, "cross-section coefficient c in the model's units");
    pair.maximum_temperature = Temperature<P, U> {
      positive_parameter<P>(t_max / tau, "cross-section maximum temperature in the model's units"),
    };
    check_neutral_cross_section(pair.a, pair.b, pair.c, pair.maximum_temperature.value());
    return pair;
}

template <std::floating_point P, UnitSystem U, std::size_t N>
VibrationElectronCoupling<P, U>
  make_vibration_electron_coupling(const VibrationElectronInput& input,
                                   const std::array<std::string_view, N>& names,
                                   const UnitScales<P, U>& scales)
{
    VibrationElectronCoupling<P, U> coupling {};
    coupling.species = known_species(names, input.species, "species");
    const CandlerFitInput& fit = input.relaxation_time;
    const std::array<double, 8> values {
      fit.a_low,
      fit.b_low,
      fit.c_low,
      fit.a_high,
      fit.b_high,
      fit.c_high,
      fit.blend_lower,
      fit.blend_upper,
    };
    for (const double v : values) {
        if (!std::isfinite(v)) {
            throw ParameterError(
              std::format("relaxation-time parameters must be finite, got {}", v));
        }
    }
    if (!(fit.blend_upper > fit.blend_lower)) {
        throw ParameterError(
          std::format("relaxation-time blend bounds must satisfy lower < upper, got {} and {}",
                      fit.blend_lower,
                      fit.blend_upper));
    }
    const UnitScales<double, U> s = double_scales(scales);
    const double tau = scale_factor<kind::temperature>(s);
    const double shift = std::log10(tau);
    const double log_k = std::log10(scale_factor(Dimension {.time = 1}, s) *
                                    scale_factor<kind::pressure>(s) / 101325.0);
    const auto convert_fit = [&](double a, double b, double c) {
        return std::array<P, 3> {
          finite_parameter<P>(a, "relaxation-time coefficient in the model's units"),
          finite_parameter<P>((2.0 * a * shift) + b,
                              "relaxation-time coefficient in the model's units"),
          finite_parameter<P>((((a * shift) + b) * shift) + c - log_k,
                              "relaxation-time coefficient in the model's units"),
        };
    };
    const std::array<P, 3> low = convert_fit(fit.a_low, fit.b_low, fit.c_low);
    const std::array<P, 3> high = convert_fit(fit.a_high, fit.b_high, fit.c_high);
    coupling.relaxation_time = {
      .a_low = low[0],
      .b_low = low[1],
      .c_low = low[2],
      .a_high = high[0],
      .b_high = high[1],
      .c_high = high[2],
      .blend_lower = Temperature<P, U> {static_cast<P>(fit.blend_lower / tau)},
      .blend_upper = Temperature<P, U> {static_cast<P>(fit.blend_upper / tau)},
    };
    if (!(coupling.relaxation_time.blend_upper.value() >
          coupling.relaxation_time.blend_lower.value())) {
        throw ParameterError("relaxation-time blend bounds are not distinct in P");
    }
    return coupling;
}

template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
ElectronImpactIonisationCoupling<P, U>
  make_ionisation_coupling(const ElectronImpactIonisationInput& input,
                           const std::array<std::string_view, N>& names,
                           const ReactionMechanism<P, U, N, R>& mechanism,
                           const UnitScales<P, U>& scales)
{
    ElectronImpactIonisationCoupling<P, U> coupling {};
    if (input.reaction < 0 || static_cast<std::size_t>(input.reaction) >= R) {
        throw ParameterError(
          std::format("reaction index must be between 0 and {}, got {}", R - 1, input.reaction));
    }
    coupling.reaction = static_cast<std::size_t>(input.reaction);
    coupling.species = known_species(names, input.species, "species");
    const std::string ion = std::format("{}+", input.species);
    const std::size_t ion_index = species_index(names, ion);
    if (ion_index == N) {
        throw ParameterError(
          std::format(R"(species "{}" has no singly charged ion "{}")", input.species, ion));
    }
    const Reaction<P, U, N>& r = mechanism.reactions[coupling.reaction];
    const int lost = r.reactants[coupling.species] - r.products[coupling.species];
    const int gained = r.products[ion_index] - r.reactants[ion_index];
    if (lost != 1 || gained != 1) {
        throw ParameterError(std::format(R"(reaction {} does not turn one "{}" into one "{}")",
                                         input.reaction,
                                         input.species,
                                         ion));
    }
    const auto energy = positive_parameter<double>(input.ionisation_energy, "ionisation energy");
    // Energy per particle.
    const Dimension particle_energy {.mass = 1, .length = 2, .time = -2, .count = -1};
    coupling.ionisation_energy = positive_parameter<P>(
      energy / constants::avogadro<double> / scale_factor(particle_energy, double_scales(scales)),
      "ionisation energy in the model's units");
    return coupling;
}

}  // namespace detail

// Energy exchange of input and electrons over the N species of the thermo
// model gas with energy modes, named in its species order, with their
// particle masses and charges from species and the reactions of mechanism
// (over the same species), in the unit system of the scales. The vibration
// mode of each species is species_vibration_mode(gas); e- is the species
// named "e-", and the model has free electrons when it has e- and
// electron_energy_mode(gas) is not -1. Throws ParameterError for the
// errors of make_energy_exchange (with the same entry names) and, with
// "electron-translation pair <index> (<partner>)", "vibration-electron
// exchange <index> (<species>)" or "electron-impact ionisation <index>
// (<species>)" in the message, when:
// - a vibration-translation, chemistry-vibration or vibration-electron
//   species has no vibration mode;
// - an electron exchange is given and the model has no free electrons;
// - a name is null or unknown; an electron-translation partner is e-, is
//   not charged for coulomb or is charged for gnoffo_neutral, or the
//   cross-section model is unknown;
// - a gnoffo_neutral coefficient is not finite, its maximum temperature
//   is not finite and positive, or its quadratic is negative/nonfinite on
//   [0, maximum_temperature] before or after parameter conversion;
// - a relaxation-time parameter is not finite or blend_upper <=
//   blend_lower;
// - an ionisation's reaction index is outside [0, R), the species named
//   with "+" appended does not exist, the reaction does not take one of the
//   species and give one of its ion (net change -1 and +1), or the
//   ionisation energy is not finite and positive;
// - a derived constant is not finite (and positive where it must be) in P.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I,
          HasEnergyModes G>
    requires(N > 0 && R > 0) && std::same_as<typename G::units, U>
MultiTemperatureExchange<P, U, N, G::mode_count, V, C, E, VE, I>
  make_multi_temperature_exchange(const EnergyExchangeInput<V, C>& input,
                                  const ElectronExchangeInput<E, VE, I>& electrons,
                                  const std::array<std::string_view, N>& names,
                                  const SpeciesData<P, U, N>& species,
                                  const ReactionMechanism<P, U, N, R>& mechanism,
                                  const G& gas,
                                  const UnitScales<P, U>& scales)
{
    constexpr std::size_t modes = G::mode_count;
    detail::check_species_names(names, species);
    MultiTemperatureExchange<P, U, N, modes, V, C, E, VE, I> exchange {};
    for (std::size_t i = 0; i < N; ++i) {
        const int mode = species_vibration_mode(gas, i);
        exchange.vibration_mode[i] =
          std::cmp_greater_equal(mode, 0) && std::cmp_less(mode, modes) ? mode : -1;
    }
    const auto check_vibration = [&](std::size_t i) {
        if (exchange.vibration_mode[i] < 0) {
            throw ParameterError(
              std::format("species \"{}\" has no vibration component in a mode", names[i]));
        }
    };
    detail::make_vibrational_exchanges(input,
                                       names,
                                       species,
                                       mechanism,
                                       scales,
                                       check_vibration,
                                       exchange.vibration_translation,
                                       exchange.chemistry_vibration);

    const std::size_t e = detail::species_index(names, "e-");
    const int electron_mode = electron_energy_mode(gas);
    const bool has_electrons =
      e < N && std::cmp_greater_equal(electron_mode, 0) && std::cmp_less(electron_mode, modes);
    exchange.electron = has_electrons ? static_cast<int>(e) : -1;
    exchange.electron_mode = has_electrons ? electron_mode : -1;
    const auto check_electrons = [&] {
        if (!has_electrons) { throw ParameterError("the model has no free electrons"); }
    };
    for (std::size_t k = 0; k < E; ++k) {
        const ElectronTranslationInput& entry = electrons.electron_translation[k];
        try {
            check_electrons();
            exchange.electron_translation[k] =
              detail::make_electron_translation_pair(entry, names, species, scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("electron-translation pair {} ({}): {}",
                                             k,
                                             detail::name_or_empty(entry.partner),
                                             error.what()));
        }
    }
    for (std::size_t k = 0; k < VE; ++k) {
        const VibrationElectronInput& entry = electrons.vibration_electron[k];
        try {
            check_electrons();
            exchange.vibration_electron[k] =
              detail::make_vibration_electron_coupling(entry, names, scales);
            check_vibration(exchange.vibration_electron[k].species);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("vibration-electron exchange {} ({}): {}",
                                             k,
                                             detail::name_or_empty(entry.species),
                                             error.what()));
        }
    }
    for (std::size_t k = 0; k < I; ++k) {
        const ElectronImpactIonisationInput& entry = electrons.electron_impact_ionisation[k];
        try {
            check_electrons();
            exchange.electron_impact_ionisation[k] =
              detail::make_ionisation_coupling(entry, names, mechanism, scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("electron-impact ionisation {} ({}): {}",
                                             k,
                                             detail::name_or_empty(entry.species),
                                             error.what()));
        }
    }

    const UnitScales<double, U> s = detail::double_scales(scales);
    const double e_charge = constants::elementary_charge<double>;
    const double coulomb_si = e_charge * e_charge /
                              (4.0 * std::numbers::pi * constants::vacuum_permittivity<double> *
                               constants::boltzmann<double>);
    const double pi = std::numbers::pi;
    exchange.coulomb_cross_section = detail::positive_parameter<P>(
      (2.0 * pi / 9.0) * coulomb_si * coulomb_si /
        scale_factor(Dimension {.length = 2, .temperature = 2, .count = -1}, s),
      "Coulomb cross-section constant in the model's units");
    exchange.coulomb_logarithm_factor = detail::positive_parameter<P>(
      9.0 / (4.0 * pi * coulomb_si * coulomb_si * coulomb_si) /
        scale_factor(Dimension {.length = -3, .temperature = -3, .count = 1}, s),
      "Coulomb logarithm constant in the model's units");
    exchange.minimum_electron_density =
      detail::positive_parameter<P>(1.0 / scale_factor(Dimension {.length = -3, .count = 1}, s),
                                    "1 m^-3 in the model's number density unit");
    exchange.mass = species.mass;
    exchange.boltzmann = boltzmann_constant(scales);
    return exchange;
}

// A thermo model with energy modes (HasEnergyModes) of N species whose
// species_mode_energy(gas, m, temp) gives e_i,m of every species at a
// temperature with value type T.
template <typename G, typename T, std::size_t N>
concept HasModeEnergies =
  HasEnergyModes<G> && Scalar<T> &&
  requires(const G& gas, std::size_t m, Temperature<T, typename G::units> temp) {
      {
          species_mode_energy(gas, m, temp)
      } -> std::convertible_to<std::array<SpecificEnergy<T, typename G::units>, N>>;
  };

// A thermo model and value type for which the mode sources of a
// MultiTemperatureExchange over N species with parameter type P can be
// evaluated.
template <typename G, typename T, typename P, std::size_t N>
concept ModeRelaxingModel = HasModeEnergies<G, T, N> && std::floating_point<P> &&
                            ScalarWith<T, P> && HasExp<T> && HasLog<T> && HasSqrt<T>;

namespace detail
{

// e_i,m_i(T) - e_i,m_i(T_m_i) (translation) and e_i,m_i(T_e) - e_i,m_i(T_m_i)
// (electron) of every species i with a vibration mode m_i, computed for the
// modes of the relaxing species of the vibration-translation pairs
// (translation, with Translation) and of the vibration-electron exchanges
// (electron, with Electron); 0 for every other species.
template <Scalar T, std::size_t N> struct RelaxationEnergies {
    std::array<T, N> translation;
    std::array<T, N> electron;
};

template <bool Translation,
          bool Electron,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I,
          typename G>
    requires ModeRelaxingModel<G, T, P, N>
CANARY_HOST_DEVICE RelaxationEnergies<T, N>
  relaxation_energies(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                      const G& gas,
                      const MultiTemperatureState<T, U, M>& state)
{
    using Real = ScalarTraits<T>::value_type;
    std::array<bool, M> translation_mode {};
    std::array<bool, M> electron_mode {};
    if constexpr (Translation) {
        for (const VibrationTranslationPair<P, U>& pair : exchange.vibration_translation) {
            translation_mode[static_cast<std::size_t>(exchange.vibration_mode[pair.species])] =
              true;
        }
    }
    if constexpr (Electron) {
        for (const VibrationElectronCoupling<P, U>& coupling : exchange.vibration_electron) {
            electron_mode[static_cast<std::size_t>(exchange.vibration_mode[coupling.species])] =
              true;
        }
    }
    RelaxationEnergies<T, N> out {};
    out.translation.fill(T(Real(0)));
    out.electron.fill(T(Real(0)));
    for (std::size_t m = 0; m < M; ++m) {
        if (!translation_mode[m] && !electron_mode[m]) { continue; }
        const std::array<SpecificEnergy<T, U>, N> at_mode =
          species_mode_energy(gas, m, Temperature<T, U> {state.mode_temp[m].value()});
        if (translation_mode[m]) {
            const std::array<SpecificEnergy<T, U>, N> at_t =
              species_mode_energy(gas, m, state.temp);
            for (std::size_t i = 0; i < N; ++i) {
                if (exchange.vibration_mode[i] == static_cast<int>(m)) {
                    out.translation[i] = T(at_t[i].value() - at_mode[i].value());
                }
            }
        }
        if (electron_mode[m]) {
            // NOLINTNEXTLINE(misc-static-assert): run-time electron mode
            CANARY_ASSERT(exchange.electron_mode >= 0);
            const auto e = static_cast<std::size_t>(exchange.electron_mode);
            const std::array<SpecificEnergy<T, U>, N> at_e =
              species_mode_energy(gas, m, Temperature<T, U> {state.mode_temp[e].value()});
            for (std::size_t i = 0; i < N; ++i) {
                if (exchange.vibration_mode[i] == static_cast<int>(m)) {
                    out.electron[i] = T(at_e[i].value() - at_mode[i].value());
                }
            }
        }
    }
    return out;
}

template <Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M> zero_mode_sources()
{
    using Real = ScalarTraits<T>::value_type;
    std::array<PowerDensity<T, U>, M> out {};
    out.fill(PowerDensity<T, U> {T(Real(0))});
    return out;
}

template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
    requires ScalarWith<T, P> && HasExp<T> && HasLog<T> && HasSqrt<T>
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M>
  vibration_translation_terms(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                              const PartialDensities<T, N>& d,
                              const std::array<T, N>& energy_difference,
                              const MultiTemperatureState<T, U, M>& state)
{
    using Real = ScalarTraits<T>::value_type;
    // t_bath = (T sum over heavy j of n_j + T_e n_e) / sum_j n_j.
    const T t = state.temp.value();
    T t_bath = t;
    if (exchange.electron >= 0) {
        const auto e = static_cast<std::size_t>(exchange.electron);
        const T t_e = state.mode_temp[static_cast<std::size_t>(exchange.electron_mode)].value();
        T heavy = T(Real(0));
        for (std::size_t i = 0; i < N; ++i) {
            if (i != e) { heavy = heavy + d.number[i]; }
        }
        const T total = T(heavy + d.number[e]);
        t_bath = T(T(T(heavy * t) + T(d.number[e] * t_e)) / total);
    }
    const std::array<T, N> inverse_tau = inverse_relaxation_times(exchange.vibration_translation,
                                                                  exchange.boltzmann,
                                                                  d.number,
                                                                  t,
                                                                  t_bath);
    std::array<T, M> q {};
    q.fill(T(Real(0)));
    for (std::size_t i = 0; i < N; ++i) {
        const int m = exchange.vibration_mode[i];
        if (m < 0) { continue; }
        const auto k = static_cast<std::size_t>(m);
        q[k] = q[k] + T(d.mass[i] * energy_difference[i] * inverse_tau[i]);
    }
    std::array<PowerDensity<T, U>, M> out {};
    for (std::size_t m = 0; m < M; ++m) { out[m] = PowerDensity<T, U> {q[m]}; }
    return out;
}

// log10(p_e tau_pe) of Candler's fit at x = log10(T_e) and temperature
// t_e; between the blend bounds the weight w is evaluated on value parts,
// and a dual number receives dw/dT_e times the derivative of T_e.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T candler_exponent(const CandlerFit<P, U>& fit, T x, T t_e)
{
    using Real = ScalarTraits<T>::value_type;
    const T f_low = T(T(T(fit.a_low * x) + fit.b_low) * x) + fit.c_low;
    const T f_high = T(T(T(fit.a_high * x) + fit.b_high) * x) + fit.c_high;
    const Real t_value = ScalarTraits<T>::value(t_e);
    const auto lower = static_cast<Real>(fit.blend_lower.value());
    const auto upper = static_cast<Real>(fit.blend_upper.value());
    if (t_value <= lower) { return f_low; }
    if (t_value >= upper) { return f_high; }
    const Real scale = std::numbers::pi_v<Real> / (Real(2) * (upper - lower));
    const Real phi = scale * (t_value - lower);
    const Real cos_phi = std::cos(phi);
    const Real w_value = cos_phi * cos_phi;
    const Real dw = -scale * std::sin(Real(2) * phi);
    const T w = T(w_value) + (T(t_e - T(t_value)) * T(dw));
    return T(w * f_low) + T(T(T(Real(1)) - w) * f_high);
}

// Mode sources of the vibration-electron exchanges, with energy_difference
// the electron differences of relaxation_energies. Requires free electrons.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
    requires ScalarWith<T, P> && HasExp<T> && HasLog<T>
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M>
  vibration_electron_terms(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                           const PartialDensities<T, N>& d,
                           const std::array<T, N>& energy_difference,
                           const MultiTemperatureState<T, U, M>& state)
{
    using Real = ScalarTraits<T>::value_type;
    // NOLINTNEXTLINE(misc-static-assert): run-time electron index
    CANARY_ASSERT(exchange.electron >= 0 && exchange.electron_mode >= 0);
    const auto e = static_cast<std::size_t>(exchange.electron);
    const auto e_mode = static_cast<std::size_t>(exchange.electron_mode);
    const T t_e = state.mode_temp[e_mode].value();
    const T p_e = T(T(d.number[e] * exchange.boltzmann.value()) * t_e);
    const Real ln10 = std::numbers::ln10_v<Real>;
    const T x = T(math::log(t_e) / T(ln10));
    std::array<T, M> q {};
    q.fill(T(Real(0)));
    for (const VibrationElectronCoupling<P, U>& coupling : exchange.vibration_electron) {
        const std::size_t p = coupling.species;
        const T f = candler_exponent(coupling.relaxation_time, x, t_e);
        const T inverse_tau = T(p_e * math::exp(T(-(f * T(ln10)))));
        const T rate = T(d.mass[p] * energy_difference[p] * inverse_tau);
        const auto m = static_cast<std::size_t>(exchange.vibration_mode[p]);
        q[m] = q[m] + rate;
        q[e_mode] = q[e_mode] - rate;
    }
    std::array<PowerDensity<T, U>, M> out {};
    for (std::size_t m = 0; m < M; ++m) { out[m] = PowerDensity<T, U> {q[m]}; }
    return out;
}

}  // namespace detail

// Mode sources of the vibration-translation pairs of exchange in a cell of
// the state, with the normalised mass fractions of mix (see
// MultiTemperatureExchange). The gas model is called through
// species_mode_energy at T and at T_m for every mode of a relaxing species.
// The status of mix is not checked.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I,
          typename G>
    requires ModeRelaxingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             (G::mode_count == M)
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M>
  vibration_translation_source(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                               const G& gas,
                               const NormalisedComposition<T, U, N>& mix,
                               const MultiTemperatureState<T, U, M>& state)
{
    return detail::vibration_translation_terms(
      exchange,
      detail::partial_densities(exchange.mass, mix, state.rho),
      detail::relaxation_energies<true, false>(exchange, gas, state).translation,
      state);
}

// Mode sources of the chemistry-vibration couplings of exchange, with rates
// the rates of progress of mechanism in a cell of the state: Q_CV of each
// coupling, at T and the temperature of the molecule's vibration mode, into
// that mode.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t R,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M>
  chemistry_vibration_source(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                             const ReactionMechanism<P, U, N, R>& mechanism,
                             const RatesOfProgress<T, U, R>& rates,
                             const MultiTemperatureState<T, U, M>& state)
{
    using Real = ScalarTraits<T>::value_type;
    std::array<T, M> q {};
    q.fill(T(Real(0)));
    for (const ChemistryVibrationCoupling<P, U>& coupling : exchange.chemistry_vibration) {
        const auto m = static_cast<std::size_t>(exchange.vibration_mode[coupling.species]);
        q[m] = q[m] + detail::coupling_energy_rate(coupling,
                                                   mechanism,
                                                   rates,
                                                   state.temp.value(),
                                                   state.mode_temp[m].value());
    }
    std::array<PowerDensity<T, U>, M> out {};
    for (std::size_t m = 0; m < M; ++m) {
        out[m] = PowerDensity<T, U> {T(exchange.boltzmann.value() * q[m])};
    }
    return out;
}

// Mode sources of the elastic electron-heavy pairs of exchange in a cell of
// the state, with the normalised mass fractions of mix: Q_ET into the
// electron mode, 0 in every other mode and without free electrons. T' of
// gnoffo_neutral and n_e' of coulomb are limited on their value parts; a
// limited value has derivative zero.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
    requires ScalarWith<T, P> && HasLog<T> && HasSqrt<T>
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M>
  electron_translation_source(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                              const NormalisedComposition<T, U, N>& mix,
                              const MultiTemperatureState<T, U, M>& state)
{
    using Real = ScalarTraits<T>::value_type;
    std::array<PowerDensity<T, U>, M> out = detail::zero_mode_sources<T, U, M>();
    if (E == 0 || exchange.electron < 0) { return out; }
    const auto e = static_cast<std::size_t>(exchange.electron);
    const auto e_mode = static_cast<std::size_t>(exchange.electron_mode);
    const P m_e = exchange.mass[e].value();
    const P k = exchange.boltzmann.value();
    const T t_e = state.mode_temp[e_mode].value();
    const Real t_e_value = ScalarTraits<T>::value(t_e);
    const detail::PartialDensities<T, N> d =
      detail::partial_densities(exchange.mass, mix, state.rho);
    const T n_e = d.number[e];
    // Coulomb cross-section (2 pi / 9) L^2 ln(9 / (4 pi n_e' L^3)), L = C / T_e.
    const T n_e_limited = ScalarTraits<T>::value(n_e) < exchange.minimum_electron_density
                            ? T(exchange.minimum_electron_density)
                            : n_e;
    const T inverse_t_e = T(T(Real(1)) / t_e);
    const T coulomb =
      T(T(exchange.coulomb_cross_section * T(inverse_t_e * inverse_t_e)) *
        math::log(T(T(exchange.coulomb_logarithm_factor * T(T(t_e * t_e) * t_e)) / n_e_limited)));
    T sum = T(Real(0));
    for (const ElectronTranslationPair<P, U>& pair : exchange.electron_translation) {
        T sigma = coulomb;
        if (pair.cross_section == ElectronCrossSection::gnoffo_neutral) {
            const P t_max = pair.maximum_temperature.value();
            const T t_limited = t_e_value > t_max ? T(t_max) : t_e;
            sigma = T(T(T(pair.c * t_limited) + pair.b) * t_limited) + pair.a;
        }
        sum =
          sum + T(T((m_e / exchange.mass[pair.partner].value()) * d.number[pair.partner]) * sigma);
    }
    const T c_e = math::sqrt(T(T((Real(8) * k) * t_e) / (std::numbers::pi_v<Real> * m_e)));
    out[e_mode] = PowerDensity<T, U> {
      T(T(T(Real(3) * n_e) * T(k * T(state.temp.value() - t_e))) * T(c_e * sum)),
    };
    return out;
}

// Mode sources of the vibration-electron exchanges of exchange in a cell of
// the state, with the normalised mass fractions of mix: each exchange adds
// Q to the vibration mode of its species and subtracts it from the
// electron mode. The gas model is called through species_mode_energy at
// T_m and T_e for every mode of a relaxing species. Between the blend
// bounds the blend weight is evaluated on value parts; a dual number
// receives dw/dT_e times the derivative of T_e.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I,
          typename G>
    requires ModeRelaxingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             (G::mode_count == M)
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M>
  vibration_electron_source(const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                            const G& gas,
                            const NormalisedComposition<T, U, N>& mix,
                            const MultiTemperatureState<T, U, M>& state)
{
    std::array<PowerDensity<T, U>, M> out = detail::zero_mode_sources<T, U, M>();
    if (VE == 0 || exchange.electron < 0) { return out; }
    const detail::RelaxationEnergies<T, N> energies =
      detail::relaxation_energies<false, true>(exchange, gas, state);
    return detail::vibration_electron_terms(
      exchange,
      detail::partial_densities(exchange.mass, mix, state.rho),
      energies.electron,
      state);
}

// Mode sources of the electron-impact ionisations of exchange, with rates
// the rates of progress of the mechanism in the cell: Q_CE into the
// electron mode, 0 in every other mode and without free electrons.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t R,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE std::array<PowerDensity<T, U>, M> electron_impact_ionisation_source(
  const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
  const RatesOfProgress<T, U, R>& rates)
{
    using Real = ScalarTraits<T>::value_type;
    std::array<PowerDensity<T, U>, M> out = detail::zero_mode_sources<T, U, M>();
    if (I == 0 || exchange.electron < 0) { return out; }
    T q = T(Real(0));
    for (const ElectronImpactIonisationCoupling<P, U>& coupling :
         exchange.electron_impact_ionisation) {
        // NOLINTNEXTLINE(misc-static-assert): run-time reaction index
        CANARY_ASSERT(coupling.reaction < R);
        const T net =
          T(rates.forward[coupling.reaction].value() - rates.backward[coupling.reaction].value());
        q = q - T(net * coupling.ionisation_energy);
    }
    out[static_cast<std::size_t>(exchange.electron_mode)] = PowerDensity<T, U> {q};
    return out;
}

// Species mass production rates and the source of every energy mode of a
// cell of a multi-temperature mixture.
template <Scalar T, UnitSystem U, std::size_t N, std::size_t M>
    requires(N > 0 && M > 0)
struct MultiTemperatureSources {
    std::array<MassProductionRate<T, U>, N> species;
    std::array<PowerDensity<T, U>, M> modes;
};

// production_rates of mechanism with the mode temperatures of the state,
// and the sum of the mode sources of exchange (vibration-translation,
// chemistry-vibration, electron-translation, vibration-electron and
// electron-impact ionisation, added in this order) in a cell of the state,
// with the normalised mass fractions of mix. The rates of progress are
// evaluated once and the mode energies once per mode and temperature.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t R,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I,
          typename G>
    requires ReactingModel<G, T, P, N> && ModeRelaxingModel<G, T, P, N> &&
             std::same_as<typename G::units, U> && (G::mode_count == M)
CANARY_HOST_DEVICE MultiTemperatureSources<T, U, N, M>
  multi_temperature_sources(const ReactionMechanism<P, U, N, R>& mechanism,
                            const MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>& exchange,
                            const G& gas,
                            const NormalisedComposition<T, U, N>& mix,
                            const MultiTemperatureState<T, U, M>& state)
{
    const RatesOfProgress<T, U, R> rates =
      rates_of_progress(mechanism, gas, mix, state.rho, state.temp, state.mode_temp);
    const detail::PartialDensities<T, N> d =
      detail::partial_densities(exchange.mass, mix, state.rho);
    const detail::RelaxationEnergies<T, N> energies =
      detail::relaxation_energies<true, true>(exchange, gas, state);
    const auto vt = detail::vibration_translation_terms(exchange, d, energies.translation, state);
    const auto cv = chemistry_vibration_source(exchange, mechanism, rates, state);
    const auto et = electron_translation_source(exchange, mix, state);
    std::array<PowerDensity<T, U>, M> ve = detail::zero_mode_sources<T, U, M>();
    if (exchange.electron >= 0) {
        ve = detail::vibration_electron_terms(exchange, d, energies.electron, state);
    }
    const auto ce = electron_impact_ionisation_source(exchange, rates);
    MultiTemperatureSources<T, U, N, M> out {};
    out.species = detail::mass_production_rates(mechanism, rates);
    for (std::size_t m = 0; m < M; ++m) {
        out.modes[m] = PowerDensity<T, U> {
          T(T(T(T(vt[m].value() + cv[m].value()) + et[m].value()) + ve[m].value()) + ce[m].value()),
        };
    }
    return out;
}

}  // namespace canary

#endif  // CANARY_ENERGY_EXCHANGE_HPP
