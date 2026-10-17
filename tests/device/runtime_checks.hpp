#ifndef CANARY_TESTS_DEVICE_RUNTIME_CHECKS_HPP
#define CANARY_TESTS_DEVICE_RUNTIME_CHECKS_HPP

// Host checking of the complete results of the numerical checks.
//
// Exact fields (status, reason, diagnostic mode, species count) are compared
// with the values a fixture declares. Every numerical field is visited by name
// with each scalar component (the value, and the derivative of a Dual). A
// numerical field without an acceptance interval is open: it is counted, and an
// open field fails the run as a missing criterion, not as a pass.

#include "../dual.hpp"

#include <array>
#include <canary/energy_exchange.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/transport/transport.hpp>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <map>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace canary::test
{

// The result of running the checks of one subsystem.
class Report
{
  public:
    // One evaluation of operation, with open_numbers scalar components that have
    // no acceptance interval.
    void evaluated(const std::string& operation, std::size_t open_numbers)
    {
        Open& open = open_[operation];
        ++open.evaluations;
        open.numbers += open_numbers;
        ++evaluations_;
    }

    // An exact field of an evaluation that differs from the declared value.
    void mismatch(std::string_view operation,
                  std::string_view fixture,
                  std::string_view field,
                  std::string_view expected,
                  std::string_view observed)
    {
        mismatches_.push_back(std::string(operation) + "\n    fixture " + std::string(fixture) +
                              ", " + std::string(field) + ": expected " + std::string(expected) +
                              ", observed " + std::string(observed));
    }

    // A failure of the backend or of the group runner itself.
    void error(std::string message) { errors_.push_back(std::move(message)); }

    void merge(const Report& other)
    {
        for (const auto& [operation, open] : other.open_) {
            open_[operation].evaluations += open.evaluations;
            open_[operation].numbers += open.numbers;
        }
        evaluations_ += other.evaluations_;
        mismatches_.insert(mismatches_.end(), other.mismatches_.begin(), other.mismatches_.end());
        errors_.insert(errors_.end(), other.errors_.begin(), other.errors_.end());
    }

    // True when at least one evaluation ran and every evaluated field met its
    // criterion, with no error.
    [[nodiscard]] bool passed() const
    {
        std::size_t open_numbers = 0;
        for (const auto& entry : open_) { open_numbers += entry.second.numbers; }
        return evaluations_ > 0 && open_numbers == 0 && mismatches_.empty() && errors_.empty();
    }

    void print(std::ostream& out) const
    {
        for (const std::string& error : errors_) { out << "error: " << error << '\n'; }
        for (const std::string& mismatch : mismatches_) { out << "mismatch: " << mismatch << '\n'; }
        std::size_t open_numbers = 0;
        for (const auto& [operation, open] : open_) {
            if (open.numbers > 0) {
                out << "open (no acceptance interval): " << operation << "\n    "
                    << open.evaluations << " evaluations, " << open.numbers
                    << " scalar components\n";
            }
            open_numbers += open.numbers;
        }
        out << evaluations_ << " evaluations, " << mismatches_.size() << " mismatches, "
            << open_numbers << " open scalar components, " << errors_.size() << " errors\n";
    }

  private:
    struct Open {
        std::size_t evaluations = 0;
        std::size_t numbers = 0;
    };
    std::map<std::string, Open> open_;
    std::size_t evaluations_ = 0;
    std::vector<std::string> mismatches_;
    std::vector<std::string> errors_;
};

// The name of the check entry Entry with its template arguments, taken from
// the compiler's function signature text; the whole text if it has no
// "Entry = " part.
template <auto Entry> std::string entry_name()
{
    const std::string_view signature = __PRETTY_FUNCTION__;
    constexpr std::string_view marker = "Entry = ";
    const std::size_t start = signature.find(marker);
    if (start == std::string_view::npos) { return std::string(signature); }
    int depth = 0;
    std::size_t end = start + marker.size();
    for (; end < signature.size(); ++end) {
        const char c = signature[end];
        if (c == '<' || c == '(' || c == '[') { ++depth; }
        if (c == '>' || c == ')') { --depth; }
        if ((c == ';' || c == ']') && depth == 0) { break; }
        if (c == ']') { --depth; }
    }
    return std::string(signature.substr(start + marker.size(), end - start - marker.size()));
}

// Numerical fields. for_each_number(result, name, visit) calls
// visit(field_name, component, value) for every scalar component of result,
// with component "value" or "derivative".

template <std::floating_point V, typename F>
void for_each_number(V value, const std::string& name, F& visit)
{ visit(name, "value", static_cast<double>(value)); }

template <std::floating_point V, typename F>
void for_each_number(const Dual<V>& value, const std::string& name, F& visit)
{
    visit(name, "value", static_cast<double>(value.real()));
    visit(name, "derivative", static_cast<double>(value.dual()));
}

template <KindTag K, typename T, typename U, typename F>
void for_each_number(const Quantity<K, T, U>& quantity, const std::string& name, F& visit);
template <typename V, std::size_t N, typename F>
void for_each_number(const std::array<V, N>& values, const std::string& name, F& visit);
template <typename T, typename U, typename F>
void for_each_number(const ThermoState<T, U>& state, const std::string& name, F& visit);
template <typename T, typename U, std::size_t M, typename F>
void for_each_number(const MultiTemperatureState<T, U, M>& state,
                     const std::string& name,
                     F& visit);
template <typename S, typename F>
void for_each_number(const ClosureOf<S>& closure, const std::string& name, F& visit);
template <typename T, typename U, typename F>
void for_each_number(const Mixture<T, U>& mixture, const std::string& name, F& visit);
template <typename T, typename U, std::size_t N, typename F>
void for_each_number(const NormalisedComposition<T, U, N>& mixture,
                     const std::string& name,
                     F& visit);
template <typename T, typename U, typename F>
void for_each_number(const AcousticProperties<T, U>& acoustics, const std::string& name, F& visit);
template <typename T, typename U, typename F>
void for_each_number(const SingleTemperatureProperties<T, U>& properties,
                     const std::string& name,
                     F& visit);
template <typename T, typename U, std::size_t M, typename F>
void for_each_number(const MultiTemperatureProperties<T, U, M>& properties,
                     const std::string& name,
                     F& visit);
template <typename T, typename U, std::size_t N, typename F>
void for_each_number(const SpeciesProperties<T, U, N>& species, const std::string& name, F& visit);
template <typename T, typename F>
void for_each_number(const Nasa9CpH<T>& curve, const std::string& name, F& visit);
template <typename T, typename U, std::size_t R, typename F>
void for_each_number(const RatesOfProgress<T, U, R>& rates, const std::string& name, F& visit);
template <typename T, typename U, std::size_t N, typename F>
void for_each_number(const TwoTemperatureSources<T, U, N>& sources,
                     const std::string& name,
                     F& visit);
template <typename T, typename U, std::size_t N, std::size_t M, typename F>
void for_each_number(const MultiTemperatureSources<T, U, N, M>& sources,
                     const std::string& name,
                     F& visit);
template <typename T, typename U, typename F>
void for_each_number(const TransportCoefficients<T, U>& transport,
                     const std::string& name,
                     F& visit);
template <typename T, typename U, std::size_t M, typename F>
void for_each_number(const MultiTemperatureTransportCoefficients<T, U, M>& transport,
                     const std::string& name,
                     F& visit);

template <KindTag K, typename T, typename U, typename F>
void for_each_number(const Quantity<K, T, U>& quantity, const std::string& name, F& visit)
{ for_each_number(quantity.value(), name, visit); }

template <typename V, std::size_t N, typename F>
void for_each_number(const std::array<V, N>& values, const std::string& name, F& visit)
{
    for (std::size_t i = 0; i < N; ++i) {
        for_each_number(values[i], name + "[" + std::to_string(i) + "]", visit);
    }
}

template <typename T, typename U, typename F>
void for_each_number(const ThermoState<T, U>& state, const std::string& name, F& visit)
{
    for_each_number(state.rho, name + ".rho", visit);
    for_each_number(state.e, name + ".e", visit);
    for_each_number(state.p, name + ".p", visit);
    for_each_number(state.temp, name + ".temp", visit);
}

template <typename T, typename U, std::size_t M, typename F>
void for_each_number(const MultiTemperatureState<T, U, M>& state, const std::string& name, F& visit)
{
    for_each_number(state.rho, name + ".rho", visit);
    for_each_number(state.e, name + ".e", visit);
    for_each_number(state.p, name + ".p", visit);
    for_each_number(state.temp, name + ".temp", visit);
    for_each_number(state.mode_e, name + ".mode_e", visit);
    for_each_number(state.mode_temp, name + ".mode_temp", visit);
}

template <typename S, typename F>
void for_each_number(const ClosureOf<S>& closure, const std::string& name, F& visit)
{ for_each_number(closure.state, name + ".state", visit); }

template <typename T, typename U, typename F>
void for_each_number(const Mixture<T, U>& mixture, const std::string& name, F& visit)
{
    for_each_number(mixture.gas_constant, name + ".gas_constant", visit);
    for_each_number(mixture.cv, name + ".cv", visit);
    for_each_number(mixture.cp, name + ".cp", visit);
}

template <typename T, typename U, std::size_t N, typename F>
void for_each_number(const NormalisedComposition<T, U, N>& mixture,
                     const std::string& name,
                     F& visit)
{
    for_each_number(mixture.fractions, name + ".fractions", visit);
    for_each_number(mixture.gas_constant, name + ".gas_constant", visit);
}

template <typename T, typename U, typename F>
void for_each_number(const AcousticProperties<T, U>& acoustics, const std::string& name, F& visit)
{
    for_each_number(acoustics.sound_speed, name + ".sound_speed", visit);
    for_each_number(acoustics.acoustic_exponent, name + ".acoustic_exponent", visit);
}

template <typename T, typename U, typename F>
void for_each_number(const SingleTemperatureProperties<T, U>& properties,
                     const std::string& name,
                     F& visit)
{
    for_each_number(properties.cv, name + ".cv", visit);
    for_each_number(properties.cp, name + ".cp", visit);
    for_each_number(properties.capacity_ratio, name + ".capacity_ratio", visit);
    for_each_number(properties.gas_constant, name + ".gas_constant", visit);
    for_each_number(properties.acoustics, name + ".acoustics", visit);
}

template <typename T, typename U, std::size_t M, typename F>
void for_each_number(const MultiTemperatureProperties<T, U, M>& properties,
                     const std::string& name,
                     F& visit)
{
    for_each_number(properties.cv_heavy, name + ".cv_heavy", visit);
    for_each_number(properties.mode_cv, name + ".mode_cv", visit);
    for_each_number(properties.gas_constant_heavy, name + ".gas_constant_heavy", visit);
    for_each_number(properties.gas_constant_electron, name + ".gas_constant_electron", visit);
    for_each_number(properties.acoustics, name + ".acoustics", visit);
}

template <typename T, typename U, std::size_t N, typename F>
void for_each_number(const SpeciesProperties<T, U, N>& species, const std::string& name, F& visit)
{
    for_each_number(species.energy, name + ".energy", visit);
    for_each_number(species.cv, name + ".cv", visit);
    for_each_number(species.enthalpy, name + ".enthalpy", visit);
    for_each_number(species.entropy, name + ".entropy", visit);
    for_each_number(species.gibbs_energy, name + ".gibbs_energy", visit);
}

template <typename T, typename F>
void for_each_number(const Nasa9CpH<T>& curve, const std::string& name, F& visit)
{
    for_each_number(curve.cp_over_r, name + ".cp_over_r", visit);
    for_each_number(curve.h_over_r, name + ".h_over_r", visit);
}

template <typename T, typename U, std::size_t R, typename F>
void for_each_number(const RatesOfProgress<T, U, R>& rates, const std::string& name, F& visit)
{
    for (const ReactionIndex r : reaction_indices(0, R)) {
        const std::string index = "[" + std::to_string(r.value()) + "]";
        for_each_number(rates.forward()[r], name + ".forward" + index, visit);
        for_each_number(rates.backward()[r], name + ".backward" + index, visit);
    }
}

template <typename T, typename U, std::size_t N, typename F>
void for_each_number(const TwoTemperatureSources<T, U, N>& sources,
                     const std::string& name,
                     F& visit)
{
    for_each_number(sources.species, name + ".species", visit);
    for_each_number(sources.vibrational_electronic, name + ".vibrational_electronic", visit);
}

template <typename T, typename U, std::size_t N, std::size_t M, typename F>
void for_each_number(const MultiTemperatureSources<T, U, N, M>& sources,
                     const std::string& name,
                     F& visit)
{
    for_each_number(sources.species, name + ".species", visit);
    for_each_number(sources.modes, name + ".modes", visit);
}

template <typename T, typename U, typename F>
void for_each_number(const TransportCoefficients<T, U>& transport,
                     const std::string& name,
                     F& visit)
{
    for_each_number(transport.viscosity, name + ".viscosity", visit);
    for_each_number(transport.conductivity, name + ".conductivity", visit);
}

template <typename T, typename U, std::size_t M, typename F>
void for_each_number(const MultiTemperatureTransportCoefficients<T, U, M>& transport,
                     const std::string& name,
                     F& visit)
{
    for_each_number(transport.viscosity, name + ".viscosity", visit);
    for_each_number(transport.conductivity, name + ".conductivity", visit);
    for_each_number(transport.mode_conductivity, name + ".mode_conductivity", visit);
}

// Declared values of the exact fields of a result.
struct ClosureExpectation {
    Status status = Status::Valid;
    ClosureReason reason = ClosureReason::None;
    OptionalModeIndex mode {};
};

struct CompositionExpectation {
    Status status = Status::Valid;
    CompositionReason reason = CompositionReason::None;
};

// Declares every numerical component of a result exactly zero, value and
// derivative (either sign of zero).
struct ZeroExpectation {};

// Compares the exact fields of result with expected. Results without exact
// fields take no expectation.
template <typename S>
void check_exact(Report& report,
                 const std::string& operation,
                 std::string_view fixture,
                 const ClosureOf<S>& closure,
                 const ClosureExpectation& expected)
{
    const auto compare = [&](std::string_view field, auto observed, auto declared) {
        if (observed != declared) {
            report.mismatch(operation,
                            fixture,
                            field,
                            std::to_string(static_cast<long long>(declared)),
                            std::to_string(static_cast<long long>(observed)));
        }
    };
    compare("status", closure.status, expected.status);
    compare("diagnostic.reason", closure.diagnostic.reason, expected.reason);
    const auto mode_text = [](OptionalModeIndex mode) {
        return mode.has_value() ? std::to_string((*mode).value()) : std::string("none");
    };
    if (closure.diagnostic.mode != expected.mode) {
        report.mismatch(operation,
                        fixture,
                        "diagnostic.mode",
                        mode_text(expected.mode),
                        mode_text(closure.diagnostic.mode));
    }
}

template <typename Composition>
    requires requires(const Composition& c) {
        { c.status } -> std::convertible_to<Status>;
        { c.reason } -> std::convertible_to<CompositionReason>;
    }
void check_exact(Report& report,
                 const std::string& operation,
                 std::string_view fixture,
                 const Composition& mixture,
                 const CompositionExpectation& expected)
{
    if (mixture.status != expected.status) {
        report.mismatch(operation,
                        fixture,
                        "status",
                        std::to_string(static_cast<int>(expected.status)),
                        std::to_string(static_cast<int>(mixture.status)));
    }
    if (mixture.reason != expected.reason) {
        report.mismatch(operation,
                        fixture,
                        "reason",
                        std::to_string(static_cast<int>(expected.reason)),
                        std::to_string(static_cast<int>(mixture.reason)));
    }
}

template <typename Result>
void check_exact(Report& report,
                 const std::string& operation,
                 std::string_view fixture,
                 const Result& result,
                 ZeroExpectation /*expected*/)
{
    auto zero = [&](const std::string& field, std::string_view component, double value) {
        if (value != 0.0) {
            char observed[32];
            std::snprintf(observed, sizeof(observed), "%a", value);
            report.mismatch(operation,
                            fixture,
                            field + " " + std::string(component),
                            "0",
                            observed);
        }
    };
    for_each_number(result, "result", zero);
}

// Records one evaluation of result: its exact fields against expected (when
// the result has any) and its numerical fields, which are all open until
// acceptance intervals exist, unless an expectation declares them exactly.
template <typename Result, typename... Expected>
void check_result(Report& report,
                  const std::string& operation,
                  [[maybe_unused]] std::string_view fixture,
                  const Result& result,
                  const Expected&... expected)
{
    (check_exact(report, operation, fixture, result, expected), ...);
    std::size_t open = 0;
    if constexpr (!(std::same_as<Expected, ZeroExpectation> || ...)) {
        auto count = [&open](const std::string&, std::string_view, double) { ++open; };
        for_each_number(result, "result", count);
    }
    report.evaluated(operation, open);
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_RUNTIME_CHECKS_HPP
