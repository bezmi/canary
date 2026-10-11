#ifndef GASLIB_ENERGY_EXCHANGE_HPP
#define GASLIB_ENERGY_EXCHANGE_HPP

// Exchange of vibrational-electronic energy in a two-temperature mixture:
// per cell, the source Q of vibrational-electronic energy per unit volume and
// time, Q = Q_VT + Q_CV.
//
// Species i has particle mass m_i, number density n_i = rho y_i / m_i and
// vibrational-electronic energy per unit mass e_ve,i (from the thermo
// model). T and T_ve are the translational and vibrational-electronic
// temperatures; neither is clamped here.
//
// Vibration-translation exchange (Landau-Teller), over the pairs (p, q) of a
// relaxing species p and a collision partner q:
//   Q_VT = sum_pairs rho y_p (e_ve,p(T) - e_ve,p(T_ve)) / tau_pq,
//   tau_pq = tau_MW,pq + tau_P,pq.
// tau_MW is the correlation of Millikan and White (J. Chem. Phys. 39, 3209,
// 1963) with the bath pressure x_q p of the partner, x_q = n_q / sum_j n_j
// and p = sum_j n_j k_B T:
//   tau_MW = (p_atm / (x_q p)) exp(a (T^(-1/3) - b) - 18.42)
// in seconds, with T in K, p in Pa and p_atm = 101325 Pa. tau_P is Park's
// high-temperature correction (C. Park, J. Thermophys. Heat Transf. 7(3),
// 1993, p. 387), zero for a pair without it:
//   tau_P = 1 / (sigma c n_q),  sigma = sigma_ref (T_ref / T)^2,
//   c = sqrt(8 k_B T / (pi mu_pq)),  mu_pq = m_p m_q / (m_p + m_q).
// Both terms are proportional to 1 / n_q, so
//   1 / tau_pq = n_q / (exp(a' T^(-1/3) + c_pq) / (k_B T) + T^(3/2) / K_pq)
// with a', c_pq and K_pq constants of the pair (see
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
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/config.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/models/two_temperature.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/reactions.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/species.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/units.hpp>
#include <numbers>
#include <string_view>

namespace gaslib
{

// Host input of the Millikan-White constants of a pair: a in K^(1/3), b in
// K^(-1/3).
struct MillikanWhiteInput {
    double a;
    double b;
};

// Host input of Park's high-temperature correction: the cross-section
// sigma_ref in m^2 at the reference temperature T_ref in K.
struct ParkCorrectionInput {
    double cross_section;
    double reference_temperature;
};

// Host input of one vibration-translation pair: the relaxing species p, the
// collision partner q, the Millikan-White constants, and whether Park's
// correction is used, with its constants (ignored without it).
struct VibrationTranslationInput {
    const char* species;
    const char* partner;
    MillikanWhiteInput millikan_white;
    bool park_correction;
    ParkCorrectionInput park;
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

// One vibration-translation pair in model units: the indices of the
// relaxing species and of the partner in the species order of the thermo
// model, and the constants of
//   1 / tau = n_q / (exp(a' T^(-1/3) + c) / (k_B T) + T^(3/2) / K)
// with T, k_B and n_q in model units and tau in the time unit of U. With
// tau_T the temperature scale, t_ref the time scale and n_ref the number
// density scale of U (all 1 in SI):
// - millikan_white_a: a' = a tau_T^(-1/3);
// - millikan_white_offset: c = ln(p_atm (1 s) in model units) - a b -
//   18.42;
// - park_coefficient: 1 / K, where K = K_SI n_ref t_ref / tau_T^(3/2) and
//   K_SI = sigma_ref T_ref^2 sqrt(8 k_B / (pi mu_pq)) in SI; 0 without
//   Park's correction.
template <std::floating_point P, UnitSystem U> struct VibrationTranslationPair {
    std::size_t species;
    std::size_t partner;
    P millikan_white_a;
    P millikan_white_offset;
    P park_coefficient;
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
    if (!std::isfinite(a) || !(a > 0.0)) {
        throw ParameterError(
          std::format("Millikan-White a must be finite and positive, got {}", a));
    }
    if (!std::isfinite(b)) {
        throw ParameterError(std::format("Millikan-White b must be finite, got {}", b));
    }
    const UnitScales<double, U> s = double_scales(scales);
    const double tau = scale_factor<kind::temperature>(s);
    // p_atm times one second, a pressure times a time.
    const Dimension pressure_time {.mass = 1, .length = -1, .time = -1};
    const double atmosphere_second = 101325.0 / scale_factor(pressure_time, s);
    pair.millikan_white_a = positive_parameter<P>(a / std::cbrt(tau), "Millikan-White a");
    pair.millikan_white_offset = static_cast<P>(std::log(atmosphere_second) - (a * b) - 18.42);
    if (!std::isfinite(pair.millikan_white_offset)) {
        throw ParameterError(
          std::format("Millikan-White a = {}, b = {} are not representable in the model's units",
                      a,
                      b));
    }
    if (input.park_correction) {
        const auto sigma =
          positive_parameter<double>(input.park.cross_section, "Park cross-section");
        const auto t_ref = positive_parameter<double>(input.park.reference_temperature,
                                                      "Park reference temperature");
        const double mass_scale = scale_factor<kind::particle_mass>(s);
        const double m_p = static_cast<double>(species.mass[pair.species].value()) * mass_scale;
        const double m_q = static_cast<double>(species.mass[pair.partner].value()) * mass_scale;
        const double reduced_mass = m_p * m_q / (m_p + m_q);
        const double k_si =
          sigma * t_ref * t_ref *
          std::sqrt(8.0 * constants::boltzmann<double> / (std::numbers::pi * reduced_mass));
        // Number density times time.
        const Dimension density_time {.length = -3, .time = 1, .count = 1};
        const double k_model = k_si * scale_factor(density_time, s) / (tau * std::sqrt(tau));
        pair.park_coefficient = positive_parameter<P>(1.0 / k_model, "Park coefficient 1 / K");
    } else {
        pair.park_coefficient = P(0);
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
    coupling.dissociation_temperature =
      Temperature<P, U> {positive_parameter<P>(theta_d / tau, "dissociation temperature")};
    coupling.vibrational_temperature =
      Temperature<P, U> {positive_parameter<P>(theta_v / tau, "vibrational temperature")};
    return coupling;
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
// - a Millikan-White a is not finite and positive or b is not finite;
// - Park's correction is used with a cross-section or reference temperature
//   that is not finite and positive;
// - a reaction index is outside [0, R), or the species has no net change of
//   coefficient in the reaction;
// - a dissociation energy or vibrational temperature is not finite and
//   positive;
// - a derived constant is not finite and positive in P.
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
    if (species.n != static_cast<int>(N)) {
        throw ParameterError(
          std::format("species data must have {} species, got {}", N, species.n));
    }
    for (std::size_t i = 0; i < N; ++i) {
        if (detail::species_index(names, names[i]) != i) {
            throw ParameterError(std::format("species \"{}\" is named twice", names[i]));
        }
    }
    const auto name_or_empty = [](const char* name) { return name != nullptr ? name : ""; };
    EnergyExchange<P, U, N, V, C> exchange {};
    for (std::size_t e = 0; e < V; ++e) {
        const VibrationTranslationInput& pair = input.vibration_translation[e];
        try {
            exchange.vibration_translation[e] =
              detail::make_vibration_translation_pair(pair, names, species, scales);
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
            exchange.chemistry_vibration[c] =
              detail::make_chemistry_vibration_coupling(coupling, names, mechanism, scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("chemistry-vibration coupling {} ({}): {}",
                                             c,
                                             name_or_empty(coupling.species),
                                             error.what()));
        }
    }
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
GASLIB_HOST_DEVICE constexpr ValueAndDerivative<V> bernoulli_remainder(V y)
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
GASLIB_HOST_DEVICE ValueAndDerivative<V> truncated_oscillator_energy(V theta_v, V theta_d, V x)
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
GASLIB_HOST_DEVICE PowerDensity<T, U>
  vibration_translation_source(const EnergyExchange<P, U, N, V, C>& exchange,
                               const G& gas,
                               const NormalisedComposition<T, U, N>& mix,
                               Density<T, U> rho,
                               Temperature<T, U> temp,
                               Temperature<T, U> temp_ve)
{
    using Real = ScalarTraits<T>::value_type;
    const T t = temp.value();
    const T t_inverse_cbrt = math::exp(T(-math::log(t) / T(Real(3))));
    const T t_three_halves = T(t * math::sqrt(t));
    const T kt = T(exchange.boltzmann.value() * t);

    // rho y_i and n_i.
    std::array<T, N> partial {};
    std::array<T, N> n {};
    for (std::size_t i = 0; i < N; ++i) {
        partial[i] = T(rho.value() * mix.fractions[i].value());
        n[i] = T(partial[i] / exchange.mass[i].value());
    }
    // sum over the partners q of species i of 1 / tau_iq.
    std::array<T, N> inverse_tau {};
    inverse_tau.fill(T(Real(0)));
    for (const VibrationTranslationPair<P, U>& pair : exchange.vibration_translation) {
        const T millikan_white =
          T(math::exp(T(T(pair.millikan_white_a * t_inverse_cbrt) + pair.millikan_white_offset)) /
            kt);
        const T park = T(pair.park_coefficient * t_three_halves);
        inverse_tau[pair.species] =
          inverse_tau[pair.species] + T(n[pair.partner] / T(millikan_white + park));
    }
    const std::array<SpecificEnergy<T, U>, N> e_t =
      species_vibrational_electronic_energy(gas, temp);
    const std::array<SpecificEnergy<T, U>, N> e_ve =
      species_vibrational_electronic_energy(gas, temp_ve);
    T q = T(Real(0));
    for (std::size_t i = 0; i < N; ++i) {
        q = q + T(partial[i] * T(e_t[i].value() - e_ve[i].value()) * inverse_tau[i]);
    }
    return PowerDensity<T, U> {q};
}

// Q_CV of the couplings of exchange, with rates the rates of progress of
// mechanism in the cell, at translational temperature temp and
// vibrational-electronic temperature temp_ve (see the top of this file).
// G(1 / T_ve - 1 / T) is evaluated on value parts by
// detail::truncated_oscillator_energy; a dual number receives the first
// derivative dG/dx times the derivative of x. With temp_ve equal to temp,
// G(1 / T_ve - 1 / T) equals G(0) exactly.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t V,
          std::size_t C>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE PowerDensity<T, U>
  chemistry_vibration_source(const EnergyExchange<P, U, N, V, C>& exchange,
                             const ReactionMechanism<P, U, N, R>& mechanism,
                             const RatesOfProgress<T, U, R>& rates,
                             Temperature<T, U> temp,
                             Temperature<T, U> temp_ve)
{
    using Real = ScalarTraits<T>::value_type;
    const T t = temp.value();
    const T t_ve = temp_ve.value();
    // x = 1 / T_ve - 1 / T, its value part, and x minus its value part: zero
    // for a floating-point type, the derivative part for a dual number.
    const T x = T(t - t_ve) / T(t * t_ve);
    const Real x_value = ScalarTraits<T>::value(x);
    const T x_derivative = x - T(x_value);
    T q = T(Real(0));
    for (const ChemistryVibrationCoupling<P, U>& coupling : exchange.chemistry_vibration) {
        // NOLINTNEXTLINE(misc-static-assert): run-time reaction index
        GASLIB_ASSERT(coupling.reaction < R);
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
        const detail::ValueAndDerivative<Real> g =
          detail::truncated_oscillator_energy(theta_v, theta_d, x_value);
        const T g_vanish = T(g.value) + (x_derivative * T(g.derivative));
        const T g_appear = T((theta_d - theta_v) / Real(2));
        q = q + T(production * g_appear) - T(loss * g_vanish);
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
GASLIB_HOST_DEVICE TwoTemperatureSources<T, U, N>
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

}  // namespace gaslib

#endif  // GASLIB_ENERGY_EXCHANGE_HPP
