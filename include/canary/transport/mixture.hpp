// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TRANSPORT_MIXTURE_HPP
#define CANARY_TRANSPORT_MIXTURE_HPP

// Transport coefficients of a mixture of N species: viscosity, thermal
// conductivity and diffusion coefficients, by two models.
//
// Species i has particle mass m_i and mole fraction
// x_i = (y_i / m_i) / sum_j (y_j / m_j), from the normalised mass fractions
// y_i of the mixture value. p is the pressure and T the translational
// temperature of the state, k_B the Boltzmann constant and
// mu_ij = m_i m_j / (m_i + m_j) the reduced mass of a pair. No mole fraction
// is thresholded: a species with x_i = 0 contributes zero to every sum.
//
// WilkeTransport. Species viscosities mu_i(T) and conductivities k_i(T) are
// NASA CEA transport fits, mixed by Wilke's rule (J. Chem. Phys. 18, 517,
// 1950):
//   mu = sum_i x_i mu_i / sum_j x_j phi_ij,
//   phi_ij = [1 + (mu_i / mu_j)^(1/2) (m_j / m_i)^(1/4)]^2 / [8 (1 + m_i / m_j)]^(1/2),
//   k = sum_i x_i k_i / sum_j x_j psi_ij,
// with psi_ij chosen by ConductivityMixing; phi_ii = psi_ii = 1. The binary
// diffusion coefficients are the first Chapman-Enskog approximation for
// Lennard-Jones molecules (Hirschfelder, Curtiss and Bird, 1954):
//   D_ij = (3/16) (2 pi (k_B T)^3 / mu_ij)^(1/2) / (p pi sigma_ij^2 Omega_D(T*)),
//   sigma_ij = (sigma_i + sigma_j) / 2, eps_ij = (eps_i eps_j)^(1/2), T* = T / eps_ij,
// with eps the well depth over k_B and the collision integral fit of
// Neufeld, Janzen and Aziz (J. Chem. Phys. 57, 1100, 1972):
//   Omega_D = A T*^(-B) + C exp(-D T*) + E exp(-F T*) + G exp(-H T*).
// The model has no conductivity of internal energy modes.
//
// GuptaYosTransport. The mixing rules of Gupta, Yos, Thompson and Lee (NASA
// RP-1232, 1990) and Gnoffo, Gupta and Shinn (NASA TP-2867, 1989), with the
// collision integrals pi Omega(l,l)_ij(T) of every pair:
//   Delta(1)_ij = (8/3) [2 mu_ij / (pi k_B T)]^(1/2) pi Omega(1,1)_ij,
//   Delta(2)_ij = (16/5) [2 mu_ij / (pi k_B T)]^(1/2) pi Omega(2,2)_ij,
//   mu = sum_i m_i x_i / sum_j x_j Delta(2)_ij,
//   k_tr = (15/4) k_B sum_i x_i / sum_j alpha_ij x_j Delta(2)_ij,
//   alpha_ij = 1 + (1 - m_i / m_j) (0.45 - 2.54 m_i / m_j) / (1 + m_i / m_j)^2,
//   k_rot = k_B sum_{i molecule} x_i / sum_j x_j Delta(1)_ij,
//   k_ve = k_B sum_{i molecule} c_i x_i / sum_j x_j Delta(1)_ij,
//   c_i = max(min(cp_i / R_i, 9/2) - 7/2, 0),
//   D_ij = k_B T / (p Delta(1)_ij).
// Molecules are the species with rotational degrees of freedom, R_i =
// k_B / m_i, and cp_i is the isobaric specific heat of species i from the
// thermo model at the vibrational-electronic temperature T_ve (T for a
// single-temperature model). k_ve is the conductivity of the
// vibrational-electronic energy; a single-temperature model adds it to
// k_tr + k_rot. A model with energy modes gives k_ve as the conductivity of
// its mode, and must have exactly one mode that holds the vibrational and
// electronic energy of every species (HasVibrationalElectronicMode).
//
// Mixture-averaged diffusion coefficients, for both models:
//   D_i = (1 - x_i) / sum_{j != i} x_j / D_ij,
// with 1 - x_i evaluated as sum_{j != i} x_j, and D_i = 0 when the sum in
// the denominator is zero (no other species present).

#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/detail/parameters.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <canary/transport/transport.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <numbers>
#include <string_view>
#include <type_traits>

namespace canary
{

// A CEA transport fit of up to S segments, in SI: for breaks[k] <= T <=
// breaks[k+1] (T in K), ln(x / SI unit) = A ln T + B / T + C / T^2 + D with
// coefficients[k] = {A, B, C, D}; x is a viscosity in Pa s or a thermal
// conductivity in W/(m K). Only the first `segments` entries are used.
template <std::size_t S> struct TransportCurveInput {
    int segments;
    std::array<double, S + 1> breaks;
    std::array<std::array<double, 4>, S> coefficients;
};

// Lennard-Jones collision diameter sigma in m and well depth epsilon / k_B in K.
struct LennardJonesInput {
    double sigma;
    double epsilon;
};

// Host input of the transport data of one species: its name, its CEA
// viscosity and conductivity fits and its Lennard-Jones parameters.
template <std::size_t S> struct SpeciesTransportInput {
    const char* name;
    TransportCurveInput<S> viscosity;
    TransportCurveInput<S> conductivity;
    LennardJonesInput lennard_jones;
};

// ln(pi Omega / m^2) = A (ln T)^3 + B (ln T)^2 + C ln T + D, T in K.
struct CollisionFitInput {
    double A;
    double B;
    double C;
    double D;
};

// The collision integrals of a pair of species, named in either order.
struct CollisionPairInput {
    const char* first;
    const char* second;
    CollisionFitInput pi_omega_11;
    CollisionFitInput pi_omega_22;
};

// How WilkeTransport mixes the species conductivities.
enum class ConductivityMixing : std::uint8_t {
    // psi_ij as phi_ij with k_i / k_j in place of mu_i / mu_j.
    conductivity_ratio,
    // psi_ij = phi_ij (Mason and Saxena, Phys. Fluids 1, 361, 1958, without
    // their factor 1.065).
    viscosity_ratio,
};

// A CEA transport fit in model units, as derived by make_transport_curve:
// for breaks[k] <= T <= breaks[k + 1], ln(x / model unit) = A ln(T / model
// temperature unit) + B / T + C / T^2 + D with coefficients[k] = {A, B, C,
// D}. The first segment k (in order) whose closed interval contains T is
// used, so a break belongs to the segment below it. low and high are the
// values at breaks[0] and breaks[segments], which the curve keeps below and
// above its ends.
template <std::floating_point P, std::size_t S>
    requires(S > 0)
struct TransportCurve {
    int segments;
    std::array<P, S + 1> breaks;
    std::array<std::array<P, 4>, S> coefficients;
    P low;
    P high;
};

namespace detail
{

// Value of a curve at temperature t (model units), with ln_t = ln t and
// inv_t = 1 / t computed by the caller once for all species. The segment
// is chosen on the value part of t; beyond the ends the value is the
// constant low or high, with derivative zero.
template <Scalar T, std::floating_point P, std::size_t S>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T transport_curve_value(const TransportCurve<P, S>& curve, T t, T ln_t, T inv_t)
{
    using V = ScalarTraits<T>::value_type;
    const V tv = ScalarTraits<T>::value(t);
    const auto last = static_cast<std::size_t>(curve.segments);
    if (tv < curve.breaks[0]) { return T(static_cast<V>(curve.low)); }
    if (tv > curve.breaks[last]) { return T(static_cast<V>(curve.high)); }
    std::size_t k = 0;
    while (k + 1 < last && !(tv <= curve.breaks[k + 1])) { ++k; }
    const std::array<P, 4>& c = curve.coefficients[k];
    const T reciprocal = T(T(c[1] + T(c[2] * inv_t)) * inv_t);
    return math::exp(T(T(T(c[0] * ln_t) + reciprocal) + c[3]));
}

// Checks of make_transport_curve; throws ParameterError.
template <std::size_t S> void check_transport_curve(const TransportCurveInput<S>& input)
{
    if (input.segments < 1 || static_cast<std::size_t>(input.segments) > S) {
        throw ParameterError(
          std::format("segment count must be between 1 and {}, got {}", S, input.segments));
    }
    const auto n = static_cast<std::size_t>(input.segments);
    for (std::size_t i = 0; i <= n; ++i) {
        if (!std::isfinite(input.breaks[i]) || !(input.breaks[i] > 0.0) ||
            (i > 0 && !(input.breaks[i] > input.breaks[i - 1]))) {
            throw ParameterError(std::format(
              "segment boundaries must be finite, positive and increasing, got {} at index {}",
              input.breaks[i],
              i));
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t k = 0; k < 4; ++k) {
            if (!std::isfinite(input.coefficients[j][k])) {
                throw ParameterError(
                  std::format("coefficient {} of segment {} is not finite", k, j));
            }
        }
    }
}

// The temperatures at which f(T) = A ln T + B / T + C / T^2 + D, with
// coefficients = {A, B, C, D}, can take its extremes on the closed segment
// [low, high]: the two ends and the roots of A T^2 - B T - 2 C = 0 strictly
// inside the segment. An entry without an interior root repeats the lower
// end.
template <std::floating_point W, std::floating_point P>
std::array<W, 4> stationary_points(const std::array<P, 4>& coefficients, P low, P high)
{
    const W a = coefficients[0];
    const W b = coefficients[1];
    const W c = coefficients[2];
    std::array<W, 4> points {W(low), W(high), W(low), W(low)};
    const auto keep_interior = [&points, low, high](std::size_t entry, W root) {
        if (root > low && root < high) { points[entry] = root; }
    };
    if (a == W(0)) {
        if (b != W(0)) { keep_interior(2, W(-2) * c / b); }
        return points;
    }
    const W discriminant = (b * b) + (W(8) * a * c);
    if (!(discriminant >= W(0))) { return points; }
    // q has the larger magnitude of the two numerators b +- sqrt(discriminant),
    // so no cancellation occurs; the second root follows from the product of
    // the roots, -2 C / A.
    const W q = (b + std::copysign(std::sqrt(discriminant), b)) / W(2);
    keep_interior(2, q / a);
    if (q != W(0)) { keep_interior(3, W(-2) * c / q); }
    return points;
}

}  // namespace detail

// Curve of a fit in a model's units, with temperature_scale tau the SI value
// of the model's temperature unit in K and value_scale s the SI value of the
// model's unit of the fitted quantity (both 1 in SI). The coefficients are
// folded in W, at least double precision, and then rounded to P:
//   A' = A, B' = B / tau, C' = C / tau^2, D' = D + A ln tau - ln s,
// and the breaks are divided by tau. The range of the stored fit is then
// certified: on every segment the fit is evaluated in W, from the stored
// coefficients and breaks, at the two ends and at the stationary points
// inside the segment, where its extremes lie. Throws ParameterError when:
// - a scale is not finite and positive, the segment count is not in [1, S],
//   the used breaks are not finite, positive and increasing, or a used
//   coefficient is not finite;
// - a stored break is not finite and positive in P, two neighbouring breaks
//   are equal in P, or a stored coefficient is not finite in P;
// - the fit is not finite and positive in P at one of the certified
//   temperatures.
template <std::size_t S, std::floating_point P>
    requires(S > 0)
TransportCurve<P, S>
  make_transport_curve(const TransportCurveInput<S>& input, P temperature_scale, P value_scale)
{
    if (!std::isfinite(temperature_scale) || !(temperature_scale > P(0))) {
        throw ParameterError(
          std::format("temperature scale must be finite and positive, got {}", temperature_scale));
    }
    if (!std::isfinite(value_scale) || !(value_scale > P(0))) {
        throw ParameterError(
          std::format("value scale must be finite and positive, got {}", value_scale));
    }
    detail::check_transport_curve(input);
    using W = std::common_type_t<double, P>;
    const auto tau = static_cast<W>(temperature_scale);
    const auto ln_tau = std::log(tau);
    const auto ln_scale = std::log(static_cast<W>(value_scale));
    const auto n = static_cast<std::size_t>(input.segments);
    TransportCurve<P, S> curve {};
    curve.segments = input.segments;
    for (std::size_t i = 0; i <= n; ++i) {
        curve.breaks[i] =
          detail::positive_parameter<P>(input.breaks[i] / tau,
                                        std::format("segment boundary {} in the model's units", i));
        if (i > 0 && !(curve.breaks[i] > curve.breaks[i - 1])) {
            throw ParameterError(
              std::format("segment boundaries {} and {} collapse after conversion", i - 1, i));
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        const std::array<double, 4>& c = input.coefficients[j];
        const std::array<W, 4> folded {
          c[0],
          c[1] / tau,
          c[2] / (tau * tau),
          c[3] + (c[0] * ln_tau) - ln_scale,
        };
        for (std::size_t k = 0; k < 4; ++k) {
            curve.coefficients[j][k] = detail::finite_parameter<P>(
              folded[k],
              std::format("coefficient {} of segment {} in the model's units", k, j));
        }
    }
    for (std::size_t j = 0; j < n; ++j) {
        const std::array<P, 4>& c = curve.coefficients[j];
        for (const W t : detail::stationary_points<W>(c, curve.breaks[j], curve.breaks[j + 1])) {
            const W inverse = W(1) / t;
            const W value = std::exp((W(c[0]) * std::log(t)) +
                                     ((W(c[1]) + (W(c[2]) * inverse)) * inverse) + W(c[3]));
            (void)detail::positive_parameter<P>(
              value,
              std::format("the fit at temperature {} in segment {}", t, j));
        }
    }
    // The ends lie inside [breaks[0], breaks[segments]], so these calls
    // evaluate the fit and do not read low or high.
    const auto value_at = [&curve](P t) {
        return detail::transport_curve_value(curve, t, std::log(t), P(1) / t);
    };
    curve.low =
      detail::positive_parameter<P>(value_at(curve.breaks[0]), "the fit at the lower end");
    curve.high =
      detail::positive_parameter<P>(value_at(curve.breaks[n]), "the fit at the upper end");
    return curve;
}

// Constants of the binary diffusion coefficient of a pair in model units:
// D_ij = factor T^(3/2) / (p Omega_D(T / well_depth)), with the well depth
// eps_ij / k_B in the temperature unit of U and log_well_depth its natural
// logarithm.
template <std::floating_point P, UnitSystem U> struct BinaryDiffusionPair {
    P factor;
    Temperature<P, U> well_depth;
    P log_well_depth;
};

// Wilke mixture transport of N species with CEA fits of up to S segments,
// in unit system U, as derived by make_wilke_transport:
// - viscosity, conductivity: the fits of every species;
// - mass: particle masses m_i;
// - mass_ratio[i][j] = (m_j / m_i)^(1/4) and
//   mixing_factor[i][j] = [8 (1 + m_i / m_j)]^(-1/2), the constant factors of
//   phi_ij;
// - mixing: the conductivity mixing rule;
// - diffusion[i][j]: the constants of D_ij (symmetric).
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires(N > 0 && S > 0)
struct WilkeTransport {
    using parameter_type = P;
    using units = U;

    std::array<TransportCurve<P, S>, N> viscosity;
    std::array<TransportCurve<P, S>, N> conductivity;
    std::array<ParticleMass<P, U>, N> mass;
    std::array<std::array<P, N>, N> mass_ratio;
    std::array<std::array<P, N>, N> mixing_factor;
    ConductivityMixing mixing;
    std::array<std::array<BinaryDiffusionPair<P, U>, N>, N> diffusion;
};

// ln(Delta / model unit) of a pair as a cubic in L = ln(T / model
// temperature unit): ((a L + b) L + c) L + d. Delta has the dimension of a
// particle mass over a dynamic viscosity.
template <std::floating_point P> struct CollisionTermFit {
    P a;
    P b;
    P c;
    P d;
};

// Delta(1) and Delta(2) of a pair.
template <std::floating_point P> struct CollisionTerms {
    CollisionTermFit<P> delta_1;
    CollisionTermFit<P> delta_2;
};

// Gupta-Yos mixture transport of N species in unit system U, as derived by
// make_gupta_yos_transport:
// - collision[i][j]: Delta(1)_ij and Delta(2)_ij (symmetric);
// - mass: particle masses m_i;
// - gas_constant: R_i = k_B / m_i;
// - alpha[i][j]: alpha_ij;
// - rotational_degrees: d_i, with molecules the species with d_i > 0;
// - boltzmann: k_B.
template <std::floating_point P, UnitSystem U, std::size_t N>
    requires(N > 0)
struct GuptaYosTransport {
    using parameter_type = P;
    using units = U;

    std::array<std::array<CollisionTerms<P>, N>, N> collision;
    std::array<ParticleMass<P, U>, N> mass;
    std::array<SpecificGasConstant<P, U>, N> gas_constant;
    std::array<std::array<P, N>, N> alpha;
    std::array<int, N> rotational_degrees;
    BoltzmannConstant<P, U> boltzmann;
};

namespace detail
{

// Checks the species data and the names of a mixture transport builder: the
// species count, that no name is listed twice and that every particle mass
// is finite and positive.
template <std::floating_point P, UnitSystem U, std::size_t N>
void check_transport_species(const std::array<std::string_view, N>& names,
                             const SpeciesData<P, U, N>& species)
{
    if (species.n != N) {
        throw ParameterError(
          std::format("species data must have {} species, got {}", N, species.n));
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (species_index(names, names[i.value()]) != OptionalSpeciesIndex {i}) {
            throw ParameterError(std::format("species \"{}\" is named twice", names[i.value()]));
        }
        const P mass = species.mass[i.value()].value();
        if (!std::isfinite(mass) || !(mass > P(0))) {
            throw ParameterError(
              std::format("species {}: particle mass must be finite and positive, got {}",
                          names[i.value()],
                          mass));
        }
    }
}

// Index into input of the entry of every species of names. Each entry must
// name a known species, and each species must have exactly one entry.
template <std::size_t N, std::size_t S, std::size_t K>
std::array<std::size_t, N> transport_entries(const std::array<SpeciesTransportInput<S>, K>& input,
                                             const std::array<std::string_view, N>& names)
{
    std::array<std::size_t, N> entry {};
    entry.fill(K);
    for (std::size_t e = 0; e < K; ++e) {
        const char* name = input[e].name;
        if (name == nullptr) {
            throw ParameterError(std::format("transport data entry {} has no name", e));
        }
        const auto i = species_index(names, name);
        if (!i) {
            throw ParameterError(std::format("unknown species \"{}\" in the transport data", name));
        }
        if (entry[(*i).value()] != K) {
            throw ParameterError(
              std::format("species \"{}\" is listed twice in the transport data", name));
        }
        entry[(*i).value()] = e;
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (entry[i.value()] == K) {
            throw ParameterError(
              std::format("no transport data for species \"{}\"", names[i.value()]));
        }
    }
    return entry;
}

// A curve of a species, with role (viscosity or conductivity) in the
// message of a ParameterError.
template <std::size_t S, std::floating_point P>
TransportCurve<P, S> species_transport_curve(const TransportCurveInput<S>& input,
                                             P temperature_scale,
                                             P value_scale,
                                             std::string_view role)
{
    try {
        return make_transport_curve(input, temperature_scale, value_scale);
    } catch (const ParameterError& error) {
        throw ParameterError(std::format("{}: {}", role, error.what()));
    }
}

// Reduced mass m_i m_j / (m_i + m_j) of a pair in SI, from the particle
// masses in the unit system of the scales.
template <std::floating_point W, std::floating_point P, UnitSystem U>
ParticleMass<W, si>
  reduced_mass(ParticleMass<W, U> first, ParticleMass<W, U> second, const UnitScales<P, U>& scales)
{
    const W a = convert<si>(first, scales).value();
    const W b = convert<si>(second, scales).value();
    return ParticleMass<W, si> {a * b / (a + b)};
}

}  // namespace detail

// Wilke transport of the N species of a thermo model, named in its species
// order, with their particle masses from species, from the transport data
// input (looked up by name, in any order), in the unit system of the scales,
// with the conductivity mixing rule mixing. The fits are folded into the
// model's units by make_transport_curve, and the constant factor of D_ij,
//   (3/16) (2 pi k_B^3 / mu_ij)^(1/2) / (pi sigma_ij^2)
// in SI, is converted to model units of T, p and D. Every stored number is
// computed in at least double precision and checked after its conversion to
// P. Throws ParameterError when:
// - species.n is not N, a name is listed twice, or mixing is not a
//   ConductivityMixing value;
// - an entry of input has no name, names an unknown species or a species
//   listed before, or a species has no entry;
// - the temperature, viscosity or conductivity scale is not finite and
//   positive in P;
// - with the species name in the message, its particle mass is not finite
//   and positive, make_transport_curve rejects its viscosity or conductivity
//   fit, or its Lennard-Jones sigma or epsilon is not finite and positive;
// - with "species <i> and <j>" in the message, a constant of the pair (the
//   factors of phi_ij, the factor of D_ij, the well depth or its logarithm)
//   is not finite, or not positive where it must be, in P.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, std::size_t K>
    requires(N > 0 && S > 0 && K > 0)
WilkeTransport<P, U, N, S>
  make_wilke_transport(const std::array<SpeciesTransportInput<S>, K>& input,
                       const std::array<std::string_view, N>& names,
                       const SpeciesData<P, U, N>& species,
                       const UnitScales<P, U>& scales,
                       ConductivityMixing mixing)
{
    using W = std::common_type_t<double, P>;
    detail::check_transport_species(names, species);
    if (mixing != ConductivityMixing::conductivity_ratio &&
        mixing != ConductivityMixing::viscosity_ratio) {
        throw ParameterError("unknown conductivity mixing rule");
    }
    const std::array<std::size_t, N> entry = detail::transport_entries(input, names);
    const auto& s = scales;
    const auto tau = scale_factor<kind::temperature>(s);
    const P temperature_scale = detail::positive_parameter<P>(tau, "temperature scale");
    const P mu_scale = detail::positive_parameter<P>(scale_factor<kind::dynamic_viscosity>(s),
                                                     "dynamic viscosity scale");
    const P k_scale = detail::positive_parameter<P>(scale_factor<kind::thermal_conductivity>(s),
                                                    "thermal conductivity scale");
    // T^(3/2) / (p D) in model units is tau^(3/2) / (p_ref D_ref) in SI.
    const auto diffusion_scale =
      tau * std::sqrt(tau) /
      (scale_factor<kind::pressure>(s) * scale_factor<kind::diffusion_coefficient>(s));

    WilkeTransport<P, U, N, S> model {};
    model.mixing = mixing;
    model.mass = species.mass;
    std::array<LennardJonesInput, N> lennard_jones {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        const SpeciesTransportInput<S>& data = input[entry[i.value()]];
        try {
            model.viscosity[i.value()] = detail::species_transport_curve(data.viscosity,
                                                                         temperature_scale,
                                                                         mu_scale,
                                                                         "viscosity");
            model.conductivity[i.value()] = detail::species_transport_curve(data.conductivity,
                                                                            temperature_scale,
                                                                            k_scale,
                                                                            "conductivity");
            const LennardJonesInput& lj = data.lennard_jones;
            if (!std::isfinite(lj.sigma) || !(lj.sigma > 0.0)) {
                throw ParameterError(
                  std::format("Lennard-Jones sigma must be finite and positive, got {}", lj.sigma));
            }
            if (!std::isfinite(lj.epsilon) || !(lj.epsilon > 0.0)) {
                throw ParameterError(
                  std::format("Lennard-Jones epsilon must be finite and positive, got {}",
                              lj.epsilon));
            }
            lennard_jones[i.value()] = lj;
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", names[i.value()], error.what()));
        }
    }
    const double boltzmann = constants::boltzmann<double>;
    for (const SpeciesIndex i : species_indices(0, N)) {
        const ParticleMass<W, U> m_i {species.mass[i.value()]};
        for (const SpeciesIndex j : species_indices(0, N)) {
            const ParticleMass<W, U> m_j {species.mass[j.value()]};
            try {
                model.mass_ratio[i.value()][j.value()] =
                  detail::positive_parameter<P>(std::sqrt(std::sqrt(m_j.value() / m_i.value())),
                                                "mass ratio factor");
                model.mixing_factor[i.value()][j.value()] = detail::positive_parameter<P>(
                  1.0 / std::sqrt(8.0 * (1.0 + (m_i.value() / m_j.value()))),
                  "mixing factor");
                const ParticleMass<W, si> mu_ij = detail::reduced_mass(m_i, m_j, scales);
                const double sigma =
                  (lennard_jones[i.value()].sigma + lennard_jones[j.value()].sigma) / 2.0;
                const auto well_depth =
                  convert<U>(Temperature<W, si> {std::sqrt(lennard_jones[i.value()].epsilon *
                                                           lennard_jones[j.value()].epsilon)},
                             scales);
                const auto factor = (3.0 / 16.0) *
                                    std::sqrt(2.0 * std::numbers::pi * boltzmann * boltzmann *
                                              boltzmann / mu_ij.value()) /
                                    (std::numbers::pi * sigma * sigma);
                model.diffusion[i.value()][j.value()] = {
                  .factor = detail::positive_parameter<P>(factor * diffusion_scale,
                                                          "binary diffusion factor"),
                  .well_depth =
                    detail::positive_thermodynamic_parameter<P>(well_depth,
                                                                "Lennard-Jones well depth"),
                  .log_well_depth =
                    detail::finite_parameter<P>(std::log(well_depth.value()),
                                                "logarithm of the Lennard-Jones well depth"),
                };
            } catch (const ParameterError& error) {
                throw ParameterError(std::format("species {} and {}: {}",
                                                 names[i.value()],
                                                 names[j.value()],
                                                 error.what()));
            }
        }
    }
    return model;
}

// As above, with ConductivityMixing::conductivity_ratio.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, std::size_t K>
    requires(N > 0 && S > 0 && K > 0)
WilkeTransport<P, U, N, S>
  make_wilke_transport(const std::array<SpeciesTransportInput<S>, K>& input,
                       const std::array<std::string_view, N>& names,
                       const SpeciesData<P, U, N>& species,
                       const UnitScales<P, U>& scales)
{
    return make_wilke_transport(input,
                                names,
                                species,
                                scales,
                                ConductivityMixing::conductivity_ratio);
}

namespace detail
{

// Index of the named species of a collision pair; throws ParameterError
// when the name is null or unknown.
template <std::size_t N>
SpeciesIndex pair_species(const std::array<std::string_view, N>& names, const char* name)
{
    if (name == nullptr) { throw ParameterError("a species has no name"); }
    const auto i = species_index(names, name);
    if (!i) { throw ParameterError(std::format("unknown species \"{}\"", name)); }
    return *i;
}

// Checks that the coefficients of a fit are finite; throws ParameterError
// naming the fit.
inline void check_collision_fit(const CollisionFitInput& fit, std::string_view role)
{
    const std::array<double, 4> c {fit.A, fit.B, fit.C, fit.D};
    const std::array<char, 4> name {'A', 'B', 'C', 'D'};
    for (std::size_t k = 0; k < 4; ++k) {
        if (!std::isfinite(c[k])) {
            throw ParameterError(
              std::format("coefficient {} of {} must be finite, got {}", name[k], role, c[k]));
        }
    }
}

// ln(Delta) in model units of a fit of pi Omega in SI, with factor 8/3 for
// Delta(1) or 16/5 for Delta(2): with ln T_SI = L + a (a = ln tau), the
// cubic of the fit in L + a, the term -(L + a) / 2 of T^(-1/2) and the
// constant ln(factor (2 mu_ij / (pi k_B))^(1/2)) - ln(Delta scale), whose
// SI terms are in log_constant, are expanded into a cubic in L. Throws
// ParameterError, naming the coefficient and the fit (role), when a
// coefficient of the cubic is not finite in P.
template <std::floating_point P>
CollisionTermFit<P> collision_term_fit(const CollisionFitInput& fit,
                                       std::common_type_t<double, P> a,
                                       std::common_type_t<double, P> log_constant,
                                       std::string_view role)
{
    const auto stored = [role](std::common_type_t<double, P> value, char name) {
        return finite_parameter<P>(
          value,
          std::format("cubic coefficient {} of ln Delta from {}", name, role));
    };
    return {
      .a = stored(fit.A, 'a'),
      .b = stored(fit.B + (3.0 * fit.A * a), 'b'),
      .c = stored(fit.C + (2.0 * fit.B * a) + (3.0 * fit.A * a * a) - 0.5, 'c'),
      .d = stored(fit.D + (fit.C * a) + (fit.B * a * a) + (fit.A * a * a * a) - (0.5 * a) +
                    log_constant,
                  'd'),
    };
}

}  // namespace detail

// Gupta-Yos transport of the N species of a thermo model, named in its
// species order, with their particle masses from species and their
// rotational degrees of freedom, from the collision integrals of the pairs
// (each pair named in either order, in any order), in the unit system of the
// scales. Every pair of species, i = j included, needs exactly one entry.
// The fits are folded with their constant factors and the temperature scale
// into the cubics of CollisionTermFit, in the model units of Delta (those of
// a particle mass over a dynamic viscosity). Every stored number is computed
// in at least double precision and checked after its conversion to P. Throws
// ParameterError when:
// - species.n is not N, or a name is listed twice;
// - with the species name in the message, a particle mass is not finite and
//   positive, a number of rotational degrees of freedom is not 0 or 2, or
//   the gas constant is not finite and positive in P;
// - with "collision pair <index> (<first>, <second>)" in the message, a
//   pair names a species that is null or unknown, is listed twice, has a
//   coefficient that is not finite, or gives a cubic coefficient that is not
//   finite in P;
// - with "species <i> and <j>" in the message, alpha_ij is not finite in P;
// - a pair of species has no entry.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t K>
    requires(N > 0 && K > 0)
GuptaYosTransport<P, U, N> make_gupta_yos_transport(const std::array<CollisionPairInput, K>& pairs,
                                                    const std::array<std::string_view, N>& names,
                                                    const SpeciesData<P, U, N>& species,
                                                    const std::array<int, N>& rotational_degrees,
                                                    const UnitScales<P, U>& scales)
{
    using W = std::common_type_t<double, P>;
    detail::check_transport_species(names, species);
    for (const SpeciesIndex i : species_indices(0, N)) {
        const int d = rotational_degrees[i.value()];
        if (d != 0 && d != 2) {
            throw ParameterError(
              std::format("species {}: rotational degrees of freedom must be 0 or 2, got {}",
                          names[i.value()],
                          d));
        }
    }
    std::array<std::array<std::size_t, N>, N> entry {};
    for (auto& row : entry) { row.fill(K); }
    const auto name_or_empty = [](const char* name) { return name != nullptr ? name : ""; };
    for (std::size_t e = 0; e < K; ++e) {
        const CollisionPairInput& pair = pairs[e];
        try {
            const SpeciesIndex i = detail::pair_species(names, pair.first);
            const SpeciesIndex j = detail::pair_species(names, pair.second);
            if (entry[i.value()][j.value()] != K) {
                throw ParameterError(
                  std::format("the pair is listed twice, first as collision pair {}",
                              entry[i.value()][j.value()]));
            }
            detail::check_collision_fit(pair.pi_omega_11, "pi Omega(1,1)");
            detail::check_collision_fit(pair.pi_omega_22, "pi Omega(2,2)");
            entry[i.value()][j.value()] = e;
            entry[j.value()][i.value()] = e;
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("collision pair {} ({}, {}): {}",
                                             e,
                                             name_or_empty(pair.first),
                                             name_or_empty(pair.second),
                                             error.what()));
        }
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(i.value(), N)) {
            if (entry[i.value()][j.value()] == K) {
                throw ParameterError(std::format("no collision integrals for the pair ({}, {})",
                                                 names[i.value()],
                                                 names[j.value()]));
            }
        }
    }

    const auto& s = scales;
    const auto ln_tau = std::log(scale_factor<kind::temperature>(s));
    // Delta has the dimension of a particle mass over a dynamic viscosity.
    const Dimension delta_dimension {.length = 1, .time = 1, .count = -1};
    const auto ln_delta_scale = std::log(scale_factor(delta_dimension, s));
    const double boltzmann = constants::boltzmann<double>;

    GuptaYosTransport<P, U, N> model {};
    model.mass = species.mass;
    model.rotational_degrees = rotational_degrees;
    model.boltzmann = boltzmann_constant(scales);
    for (const SpeciesIndex i : species_indices(0, N)) {
        const ParticleMass<W, U> m_i {species.mass[i.value()]};
        model.gas_constant[i.value()] = detail::positive_thermodynamic_parameter<P>(
          specific_gas_constant(BoltzmannConstant<W, U> {model.boltzmann}, m_i),
          std::format("species {}: gas constant", names[i.value()]));
        for (const SpeciesIndex j : species_indices(0, N)) {
            const ParticleMass<W, U> m_j {species.mass[j.value()]};
            const W r = m_i.value() / m_j.value();
            model.alpha[i.value()][j.value()] = detail::finite_parameter<P>(
              1.0 + ((1.0 - r) * (0.45 - (2.54 * r)) / ((1.0 + r) * (1.0 + r))),
              std::format("species {} and {}: alpha", names[i.value()], names[j.value()]));
            const ParticleMass<W, si> mu_ij = detail::reduced_mass(m_i, m_j, scales);
            const auto log_root =
              0.5 * std::log(2.0 * mu_ij.value() / (std::numbers::pi * boltzmann));
            const std::size_t e = entry[i.value()][j.value()];
            const CollisionPairInput& pair = pairs[e];
            try {
                model.collision[i.value()][j.value()] = {
                  .delta_1 =
                    detail::collision_term_fit<P>(pair.pi_omega_11,
                                                  ln_tau,
                                                  std::log(8.0 / 3.0) + log_root - ln_delta_scale,
                                                  "pi Omega(1,1)"),
                  .delta_2 =
                    detail::collision_term_fit<P>(pair.pi_omega_22,
                                                  ln_tau,
                                                  std::log(16.0 / 5.0) + log_root - ln_delta_scale,
                                                  "pi Omega(2,2)"),
                };
            } catch (const ParameterError& error) {
                throw ParameterError(std::format("collision pair {} ({}, {}): {}",
                                                 e,
                                                 pair.first,
                                                 pair.second,
                                                 error.what()));
            }
        }
    }
    return model;
}

// A mixture model whose per-cell mixture value is the NormalisedComposition
// of N species with value type T, from which the mixture transport models
// take the mole fractions.
template <typename G, typename T, std::size_t N>
concept HasNormalisedComposition =
  HasMixture<G, T> && std::same_as<mixture_t<G, T>, NormalisedComposition<T, typename G::units, N>>;

// A thermo model of N species whose species_isobaric_specific_heat(gas,
// temp) gives cp_i of every species at a temperature with value type T.
template <typename G, typename T, std::size_t N>
concept HasSpeciesIsobaricSpecificHeat =
  GasModelType<G> && Scalar<T> && requires(const G& gas, Temperature<T, typename G::units> temp) {
      {
          species_isobaric_specific_heat(gas, temp)
      } -> std::convertible_to<std::array<SpecificHeat<T, typename G::units>, N>>;
  };

// A thermo model whose single internal energy mode holds the vibrational and
// electronic energy of every species at T_ve = mode_temp[0], without free
// electrons (Park's two-temperature model). The model declares this with a
// static constexpr bool vibrational_electronic_mode equal to true.
template <typename G>
concept HasVibrationalElectronicMode = detail::HasModeCount<G> && (G::mode_count == 1) && requires {
    requires std::bool_constant<G::vibrational_electronic_mode>::value;
};

// Binary diffusion coefficients D_ij of N species (symmetric).
template <Scalar T, UnitSystem U, std::size_t N>
using BinaryDiffusionCoefficients = std::array<std::array<DiffusionCoefficient<T, U>, N>, N>;

namespace detail
{

// Mole fractions x_i = (y_i / m_i) / sum_j (y_j / m_j).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE std::array<T, N> mole_fractions(const std::array<ParticleMass<P, U>, N>& mass,
                                                   const NormalisedComposition<T, U, N>& mix)
{
    using V = ScalarTraits<T>::value_type;
    std::array<T, N> x {};
    T sum = T(V(0));
    for (const SpeciesIndex i : species_indices(0, N)) {
        x[i.value()] = T(mix.fractions[i.value()].value() / mass[i.value()].value());
        sum = sum + x[i.value()];
    }
    const T inv_sum = T(V(1)) / sum;
    for (const SpeciesIndex i : species_indices(0, N)) { x[i.value()] = x[i.value()] * inv_sum; }
    return x;
}

// sum_j x_j phi_ij for every i, where phi_ii = 1 and, for j != i,
// phi_ij = [1 + (root_i / root_j) (m_j / m_i)^(1/4)]^2 [8 (1 + m_i / m_j)]^(-1/2)
// with root_i the square root of the mixed species value (mu_i or k_i).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE std::array<T, N> wilke_denominators(const WilkeTransport<P, U, N, S>& model,
                                                       const std::array<T, N>& x,
                                                       const std::array<T, N>& root)
{
    using V = ScalarTraits<T>::value_type;
    std::array<T, N> inverse_root {};
    for (const SpeciesIndex j : species_indices(0, N)) {
        inverse_root[j.value()] = T(V(1)) / root[j.value()];
    }
    std::array<T, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        T sum = x[i.value()];
        for (const SpeciesIndex j : species_indices(0, N)) {
            if (j.value() == i.value()) { continue; }
            const T w = T(V(1) + T(T(root[i.value()] * inverse_root[j.value()]) *
                                   model.mass_ratio[i.value()][j.value()]));
            sum = sum + T(x[j.value()] * T(T(w * w) * model.mixing_factor[i.value()][j.value()]));
        }
        out[i.value()] = sum;
    }
    return out;
}

// Viscosity and conductivity of a mixture, as values in model units.
template <Scalar T> struct ViscosityAndConductivity {
    T viscosity;
    T conductivity;
};

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE ViscosityAndConductivity<T>
  wilke_values(const WilkeTransport<P, U, N, S>& model,
               const NormalisedComposition<T, U, N>& mix,
               Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const std::array<T, N> x = mole_fractions(model.mass, mix);
    const T t = temp.value();
    const T ln_t = math::log(t);
    const T inv_t = T(V(1)) / t;
    std::array<T, N> mu {};
    std::array<T, N> k {};
    std::array<T, N> root {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        mu[i.value()] = transport_curve_value(model.viscosity[i.value()], t, ln_t, inv_t);
        k[i.value()] = transport_curve_value(model.conductivity[i.value()], t, ln_t, inv_t);
        root[i.value()] = math::sqrt(mu[i.value()]);
    }
    const std::array<T, N> mu_sums = wilke_denominators(model, x, root);
    std::array<T, N> k_sums = mu_sums;
    if (model.mixing == ConductivityMixing::conductivity_ratio) {
        for (const SpeciesIndex i : species_indices(0, N)) {
            root[i.value()] = math::sqrt(k[i.value()]);
        }
        k_sums = wilke_denominators(model, x, root);
    }
    T mu_mix = T(V(0));
    T k_mix = T(V(0));
    for (const SpeciesIndex i : species_indices(0, N)) {
        mu_mix = mu_mix + T(T(x[i.value()] * mu[i.value()]) / mu_sums[i.value()]);
        k_mix = k_mix + T(T(x[i.value()] * k[i.value()]) / k_sums[i.value()]);
    }
    return {.viscosity = mu_mix, .conductivity = k_mix};
}

// Omega_D of a WilkeTransport pair (Neufeld, Janzen and Aziz) at t, with
// ln_t = ln t.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T wilke_omega(const BinaryDiffusionPair<P, U>& pair, T t, T ln_t)
{
    using V = ScalarTraits<T>::value_type;
    const V a = V(1.06036);
    const V b = V(0.15610);
    const V c = V(0.19300);
    const V d = V(0.47635);
    const V e = V(1.03587);
    const V f = V(1.52996);
    const V g = V(1.76474);
    const V h = V(3.89411);
    const T t_star = T(t / pair.well_depth.value());
    const T ln_t_star = T(ln_t - pair.log_well_depth);
    // The exponentials of -D T* etc. underflow to zero at large T*.
    return T(T(T(a * math::exp(T(-(b * ln_t_star)))) + T(c * math::exp(T(-(d * t_star))))) +
             T(T(e * math::exp(T(-(f * t_star)))) + T(g * math::exp(T(-(h * t_star))))));
}

// D_ij of WilkeTransport for j >= i (j > i without Diagonal), mirrored to
// j < i; with Diagonal false, D_ii is zero.
template <bool Diagonal,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE std::array<std::array<T, N>, N> wilke_binary_diffusion(
  const WilkeTransport<P, U, N, S>& model, Temperature<T, U> temp, Pressure<T, U> p)
{
    const T t = temp.value();
    const T ln_t = math::log(t);
    // T^(3/2) / p.
    const T t_over_p = T(T(t * math::sqrt(t)) / p.value());
    std::array<std::array<T, N>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(Diagonal ? i.value() : i.value() + 1, N)) {
            const BinaryDiffusionPair<P, U>& pair = model.diffusion[i.value()][j.value()];
            out[i.value()][j.value()] = T(T(pair.factor * t_over_p) / wilke_omega(pair, t, ln_t));
            out[j.value()][i.value()] = out[i.value()][j.value()];
        }
    }
    return out;
}

// 1 / D_ij of WilkeTransport for j != i, mirrored; the diagonal is zero.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE std::array<std::array<T, N>, N> wilke_inverse_binary_diffusion(
  const WilkeTransport<P, U, N, S>& model, Temperature<T, U> temp, Pressure<T, U> p)
{
    const T t = temp.value();
    const T ln_t = math::log(t);
    // T^(3/2) / p.
    const T t_over_p = T(T(t * math::sqrt(t)) / p.value());
    std::array<std::array<T, N>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(i.value() + 1, N)) {
            const BinaryDiffusionPair<P, U>& pair = model.diffusion[i.value()][j.value()];
            out[i.value()][j.value()] = T(wilke_omega(pair, t, ln_t) / T(pair.factor * t_over_p));
            out[j.value()][i.value()] = out[i.value()][j.value()];
        }
    }
    return out;
}

// exp of the cubic of a collision term at ln_t.
template <Scalar T, std::floating_point P>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T collision_term(const CollisionTermFit<P>& fit, T ln_t)
{
    T r = T(fit.a * ln_t);
    r = T(T(r + fit.b) * ln_t);
    r = T(T(r + fit.c) * ln_t);
    return math::exp(T(r + fit.d));
}

// Delta(1)_ij and, with Viscous, Delta(2)_ij of every pair at temperature
// temp, from one logarithm; each pair is evaluated once.
template <Scalar T, std::size_t N> struct CollisionTermValues {
    std::array<std::array<T, N>, N> delta_1;
    std::array<std::array<T, N>, N> delta_2;
};

template <bool Viscous, Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE CollisionTermValues<T, N>
  collision_term_values(const GuptaYosTransport<P, U, N>& model, Temperature<T, U> temp)
{
    const T ln_t = math::log(temp.value());
    CollisionTermValues<T, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(i.value(), N)) {
            const CollisionTerms<P>& terms = model.collision[i.value()][j.value()];
            out.delta_1[i.value()][j.value()] = collision_term(terms.delta_1, ln_t);
            out.delta_1[j.value()][i.value()] = out.delta_1[i.value()][j.value()];
            if constexpr (Viscous) {
                out.delta_2[i.value()][j.value()] = collision_term(terms.delta_2, ln_t);
                out.delta_2[j.value()][i.value()] = out.delta_2[i.value()][j.value()];
            }
        }
    }
    return out;
}

// mu, k_tr + k_rot and k_ve of GuptaYosTransport, as values in model units,
// with cp_i of the thermo model at temp_ve. An absent species (zero mole
// fraction) contributes exactly zero to every sum, also when its collision
// terms are not finite at temp; its composition derivative is then dropped
// (safe_product, safe_divide).
template <Scalar T> struct GuptaYosValues {
    T viscosity;
    T conductivity;
    T vibrational_electronic_conductivity;
};

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, typename G>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasMin<T> && HasMax<T> &&
             HasSpeciesIsobaricSpecificHeat<G, T, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE GuptaYosValues<T> gupta_yos_values(const GuptaYosTransport<P, U, N>& model,
                                                      const G& gas,
                                                      const NormalisedComposition<T, U, N>& mix,
                                                      Temperature<T, U> temp,
                                                      Temperature<T, U> temp_ve)
{
    using V = ScalarTraits<T>::value_type;
    const std::array<T, N> x = mole_fractions(model.mass, mix);
    const CollisionTermValues<T, N> delta = collision_term_values<true>(model, temp);
    const std::array<SpecificHeat<T, U>, N> cp = species_isobaric_specific_heat(gas, temp_ve);
    T mu = T(V(0));
    T translational = T(V(0));
    T rotational = T(V(0));
    T vibrational = T(V(0));
    for (const SpeciesIndex i : species_indices(0, N)) {
        T sum_1 = safe_product(x[0], delta.delta_1[i.value()][0]);
        T sum_2 = safe_product(x[0], delta.delta_2[i.value()][0]);
        T sum_alpha = T(model.alpha[i.value()][0] * sum_2);
        for (const SpeciesIndex j : species_indices(1, N)) {
            const T term_2 = safe_product(x[j.value()], delta.delta_2[i.value()][j.value()]);
            sum_1 = sum_1 + safe_product(x[j.value()], delta.delta_1[i.value()][j.value()]);
            sum_2 = sum_2 + term_2;
            sum_alpha = sum_alpha + T(model.alpha[i.value()][j.value()] * term_2);
        }
        mu = mu + safe_divide(T(model.mass[i.value()].value() * x[i.value()]), sum_2);
        translational = translational + safe_divide(x[i.value()], sum_alpha);
        if (model.rotational_degrees[i.value()] > 0) {
            const T q = safe_divide(x[i.value()], sum_1);
            const T ratio = T(cp[i.value()].value() / model.gas_constant[i.value()].value());
            const T c = math::max(T(math::min(ratio, T(V(4.5))) - V(3.5)), T(V(0)));
            rotational = rotational + q;
            vibrational = vibrational + T(c * q);
        }
    }
    const P k_b = model.boltzmann.value();
    return {
      .viscosity = mu,
      .conductivity = T(T((P(15) / P(4)) * k_b * translational) + T(k_b * rotational)),
      .vibrational_electronic_conductivity = T(k_b * vibrational),
    };
}

// D_ij = k_B T / (p Delta(1)_ij) of GuptaYosTransport, i = j included.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE std::array<std::array<T, N>, N> gupta_yos_binary_diffusion(
  const GuptaYosTransport<P, U, N>& model, Temperature<T, U> temp, Pressure<T, U> p)
{
    const CollisionTermValues<T, N> delta = collision_term_values<false>(model, temp);
    const T kt_over_p = T(T(model.boltzmann.value() * temp.value()) / p.value());
    std::array<std::array<T, N>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(0, N)) {
            out[i.value()][j.value()] = T(kt_over_p / delta.delta_1[i.value()][j.value()]);
        }
    }
    return out;
}

// D_i = sum_{j != i} x_j / sum_{j != i} x_j / D_ij, from inv_d_ij = 1 / D_ij;
// 0 when the second sum is zero.
template <Scalar T, UnitSystem U, std::size_t N>
CANARY_HOST_DEVICE std::array<DiffusionCoefficient<T, U>, N>
  mixture_diffusion(const std::array<T, N>& x, const std::array<std::array<T, N>, N>& inv_d)
{
    using V = ScalarTraits<T>::value_type;
    std::array<DiffusionCoefficient<T, U>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        T others = T(V(0));
        T sum = T(V(0));
        for (const SpeciesIndex j : species_indices(0, N)) {
            if (j.value() == i.value()) { continue; }
            others = others + x[j.value()];
            sum = sum + T(x[j.value()] * inv_d[i.value()][j.value()]);
        }
        out[i.value()] = DiffusionCoefficient<T, U> {
          ScalarTraits<T>::value(sum) == V(0) ? T(V(0)) : T(others / sum),
        };
    }
    return out;
}

template <Scalar T, UnitSystem U, std::size_t N>
CANARY_HOST_DEVICE BinaryDiffusionCoefficients<T, U, N>
  as_diffusion_coefficients(const std::array<std::array<T, N>, N>& d)
{
    BinaryDiffusionCoefficients<T, U, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(0, N)) {
            out[i.value()][j.value()] = DiffusionCoefficient<T, U> {d[i.value()][j.value()]};
        }
    }
    return out;
}

// Complete one diffusion row after accumulating partners in species order.
// A zero value-part denominator retains the constant-zero result convention.
template <Scalar T, UnitSystem U, std::size_t N>
CANARY_HOST_DEVICE DiffusionCoefficient<T, U>
  mixture_diffusion_row(const std::array<T, N>& x, SpeciesIndex i, T sum)
{
    using V = ScalarTraits<T>::value_type;
    T others = T(V(0));
    for (const SpeciesIndex j : species_indices(0, N)) {
        if (j == i) { continue; }
        others = others + x[j.value()];
    }
    return DiffusionCoefficient<T, U> {
      ScalarTraits<T>::value(sum) == V(0) ? T(V(0)) : T(others / sum)};
}

// Accumulate direct p Delta(1)_ij / (k_B T) values without pair matrices.
// An absent species contributes exactly zero to every sum, also when its
// collision term is not finite at temp (safe_product).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE std::array<DiffusionCoefficient<T, U>, N>
  mixture_diffusion(const GuptaYosTransport<P, U, N>& model,
                    const std::array<T, N>& x,
                    Temperature<T, U> temp,
                    Pressure<T, U> p)
{
    using V = ScalarTraits<T>::value_type;
    const T ln_t = math::log(temp.value());
    const T p_over_kt = T(p.value() / T(model.boltzmann.value() * temp.value()));
    std::array<T, N> sums {};
    for (auto& sum : sums) { sum = T(V(0)); }
    std::array<DiffusionCoefficient<T, U>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        for (const SpeciesIndex j : species_indices(i.value() + 1, N)) {
            const T delta = collision_term(model.collision[i.value()][j.value()].delta_1, ln_t);
            const T inverse = T(delta * p_over_kt);
            sums[i.value()] = sums[i.value()] + safe_product(x[j.value()], inverse);
            sums[j.value()] = sums[j.value()] + safe_product(x[i.value()], inverse);
        }
        out[i.value()] = mixture_diffusion_row<T, U>(x, i, sums[i.value()]);
    }
    return out;
}

}  // namespace detail

// Viscosity and conductivity of WilkeTransport at a state of a
// single-temperature mixture model, at the state temperature.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          GasModelType G,
          Scalar T>
    requires HasNormalisedComposition<G, T, N> && std::same_as<state_t<G, T>, ThermoState<T, U>> &&
             ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE TransportCoefficients<T, U> transport(const WilkeTransport<P, U, N, S>& model,
                                                         const G& /*gas*/,
                                                         const NormalisedComposition<T, U, N>& mix,
                                                         const ThermoState<T, U>& state)
{
    const auto values = detail::wilke_values(model, mix, state.temp);
    return {
      .viscosity = DynamicViscosity<T, U> {values.viscosity},
      .conductivity = ThermalConductivity<T, U> {values.conductivity},
    };
}

// As above at a state with M internal energy modes, at the translational
// temperature; every mode conductivity is zero.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          GasModelType G,
          Scalar T,
          std::size_t M>
    requires HasNormalisedComposition<G, T, N> &&
             std::same_as<state_t<G, T>, MultiTemperatureState<T, U, M>> && ScalarWith<T, P> &&
             HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE MultiTemperatureTransportCoefficients<T, U, M>
  transport(const WilkeTransport<P, U, N, S>& model,
            const G& /*gas*/,
            const NormalisedComposition<T, U, N>& mix,
            const MultiTemperatureState<T, U, M>& state)
{
    using V = ScalarTraits<T>::value_type;
    const auto values = detail::wilke_values(model, mix, state.temp);
    MultiTemperatureTransportCoefficients<T, U, M> out {
      .viscosity = DynamicViscosity<T, U> {values.viscosity},
      .conductivity = ThermalConductivity<T, U> {values.conductivity},
      .mode_conductivity = {},
    };
    out.mode_conductivity.fill(ThermalConductivity<T, U> {T(V(0))});
    return out;
}

// Viscosity and conductivity k_tr + k_rot + k_ve of GuptaYosTransport at a
// state of a single-temperature mixture model, with T_ve = T.
template <std::floating_point P, UnitSystem U, std::size_t N, GasModelType G, Scalar T>
    requires HasNormalisedComposition<G, T, N> && std::same_as<state_t<G, T>, ThermoState<T, U>> &&
             HasSpeciesIsobaricSpecificHeat<G, T, N> && ScalarWith<T, P> && HasLog<T> &&
             HasExp<T> && HasMin<T> && HasMax<T>
CANARY_HOST_DEVICE TransportCoefficients<T, U> transport(const GuptaYosTransport<P, U, N>& model,
                                                         const G& gas,
                                                         const NormalisedComposition<T, U, N>& mix,
                                                         const ThermoState<T, U>& state)
{
    const auto values = detail::gupta_yos_values(model, gas, mix, state.temp, state.temp);
    return {
      .viscosity = DynamicViscosity<T, U> {values.viscosity},
      .conductivity =
        ThermalConductivity<T, U> {
          T(values.conductivity + values.vibrational_electronic_conductivity),
        },
    };
}

// Viscosity, conductivity k_tr + k_rot and the mode conductivity k_ve of
// GuptaYosTransport at a state of a mixture model with one
// vibrational-electronic mode (HasVibrationalElectronicMode), with
// T_ve = mode_temp[0].
template <std::floating_point P, UnitSystem U, std::size_t N, GasModelType G, Scalar T>
    requires HasNormalisedComposition<G, T, N> && HasVibrationalElectronicMode<G> &&
             std::same_as<state_t<G, T>, MultiTemperatureState<T, U, 1>> &&
             HasSpeciesIsobaricSpecificHeat<G, T, N> && ScalarWith<T, P> && HasLog<T> &&
             HasExp<T> && HasMin<T> && HasMax<T>
CANARY_HOST_DEVICE MultiTemperatureTransportCoefficients<T, U, 1>
  transport(const GuptaYosTransport<P, U, N>& model,
            const G& gas,
            const NormalisedComposition<T, U, N>& mix,
            const MultiTemperatureState<T, U, 1>& state)
{
    const auto values = detail::gupta_yos_values(model,
                                                 gas,
                                                 mix,
                                                 state.temp,
                                                 Temperature<T, U> {state.mode_temp[0].value()});
    return {
      .viscosity = DynamicViscosity<T, U> {values.viscosity},
      .conductivity = ThermalConductivity<T, U> {values.conductivity},
      .mode_conductivity = {ThermalConductivity<T, U> {values.vibrational_electronic_conductivity}},
    };
}

// Binary diffusion coefficients D_ij of WilkeTransport, i = j included, at
// the temperature and pressure of a state of a mixture model.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          GasModelType G,
          Scalar T,
          ThermoStateType State>
    requires HasNormalisedComposition<G, T, N> && std::same_as<state_t<G, T>, State> &&
             ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE BinaryDiffusionCoefficients<T, U, N>
  binary_diffusion_coefficients(const WilkeTransport<P, U, N, S>& model,
                                const G& /*gas*/,
                                const NormalisedComposition<T, U, N>& /*mix*/,
                                const State& state)
{
    return detail::as_diffusion_coefficients<T, U>(
      detail::wilke_binary_diffusion<true>(model, state.temp, state.p));
}

// Mixture-averaged diffusion coefficients D_i of WilkeTransport at a state
// of a mixture model.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          GasModelType G,
          Scalar T,
          ThermoStateType State>
    requires HasNormalisedComposition<G, T, N> && std::same_as<state_t<G, T>, State> &&
             ScalarWith<T, P> && HasLog<T> && HasExp<T> && HasSqrt<T>
CANARY_HOST_DEVICE std::array<DiffusionCoefficient<T, U>, N>
  mixture_diffusion_coefficients(const WilkeTransport<P, U, N, S>& model,
                                 const G& /*gas*/,
                                 const NormalisedComposition<T, U, N>& mix,
                                 const State& state)
{
    return detail::mixture_diffusion<T, U>(
      detail::mole_fractions(model.mass, mix),
      detail::wilke_inverse_binary_diffusion(model, state.temp, state.p));
}

// Binary diffusion coefficients D_ij = k_B T / (p Delta(1)_ij) of
// GuptaYosTransport, i = j included, at the temperature and pressure of a
// state of a mixture model.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          GasModelType G,
          Scalar T,
          ThermoStateType State>
    requires HasNormalisedComposition<G, T, N> && std::same_as<state_t<G, T>, State> &&
             ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE BinaryDiffusionCoefficients<T, U, N>
  binary_diffusion_coefficients(const GuptaYosTransport<P, U, N>& model,
                                const G& /*gas*/,
                                const NormalisedComposition<T, U, N>& /*mix*/,
                                const State& state)
{
    return detail::as_diffusion_coefficients<T, U>(
      detail::gupta_yos_binary_diffusion(model, state.temp, state.p));
}

// Mixture-averaged diffusion coefficients D_i of GuptaYosTransport at a
// state of a mixture model.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          GasModelType G,
          Scalar T,
          ThermoStateType State>
    requires HasNormalisedComposition<G, T, N> && std::same_as<state_t<G, T>, State> &&
             ScalarWith<T, P> && HasLog<T> && HasExp<T>
CANARY_HOST_DEVICE std::array<DiffusionCoefficient<T, U>, N>
  mixture_diffusion_coefficients(const GuptaYosTransport<P, U, N>& model,
                                 const G& /*gas*/,
                                 const NormalisedComposition<T, U, N>& mix,
                                 const State& state)
{
    return detail::mixture_diffusion(
      model, detail::mole_fractions(model.mass, mix), state.temp, state.p);
}

namespace detail
{

template <std::size_t S> struct InputStorage<SpeciesTransportInput<S>> {
    OwnedSpeciesName name {};
    TransportCurveInput<S> viscosity {};
    TransportCurveInput<S> conductivity {};
    LennardJonesInput lennard_jones {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const SpeciesTransportInput<S>& input) :
        name(input.name), viscosity(input.viscosity), conductivity(input.conductivity),
        lennard_jones(input.lennard_jones)
    {
    }
    SpeciesTransportInput<S> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .name = name.c_str(),
          .viscosity = viscosity,
          .conductivity = conductivity,
          .lennard_jones = lennard_jones,
        };
    }
};

template <> struct InputStorage<CollisionPairInput> {
    OwnedSpeciesName first {};
    OwnedSpeciesName second {};
    CollisionFitInput pi_omega_11 {};
    CollisionFitInput pi_omega_22 {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const CollisionPairInput& input) :
        first(input.first), second(input.second), pi_omega_11(input.pi_omega_11),
        pi_omega_22(input.pi_omega_22)
    {
    }
    CollisionPairInput borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .first = first.c_str(),
          .second = second.c_str(),
          .pi_omega_11 = pi_omega_11,
          .pi_omega_22 = pi_omega_22,
        };
    }
};

}  // namespace detail

// Builds from owned and direct inputs, keeping borrowed views inside this call.
template <typename... Explicit, typename... Inputs>
    requires detail::HasOwnedInput<Inputs...>
auto make_wilke_transport(const Inputs&... inputs)
{
    return detail::with_inputs(
      [](const auto&... borrowed) { return make_wilke_transport<Explicit...>(borrowed...); },
      inputs...);
}

// Builds from owned and direct inputs, keeping borrowed views inside this call.
template <typename... Explicit, typename... Inputs>
    requires detail::HasOwnedInput<Inputs...>
auto make_gupta_yos_transport(const Inputs&... inputs)
{
    return detail::with_inputs(
      [](const auto&... borrowed) { return make_gupta_yos_transport<Explicit...>(borrowed...); },
      inputs...);
}

}  // namespace canary

#endif  // CANARY_TRANSPORT_MIXTURE_HPP
