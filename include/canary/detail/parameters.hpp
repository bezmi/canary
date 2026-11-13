#ifndef CANARY_DETAIL_PARAMETERS_HPP
#define CANARY_DETAIL_PARAMETERS_HPP

// Host-side checks of the numbers that a builder stores in a model. Each
// function receives a value in a host precision W, after unit conversion,
// and returns it in the parameter type P of the model, or throws
// ParameterError naming the value.

#include <canary/build.hpp>
#include <canary/quantity.hpp>
#include <cmath>
#include <concepts>
#include <format>
#include <limits>
#include <string_view>
#include <type_traits>

namespace canary::detail
{

// A quantity that is finite in P, with its kind and units; throws
// ParameterError otherwise.
template <std::floating_point P, KindTag K, std::floating_point W, UnitTag U>
Quantity<K, P, U> thermodynamic_parameter(Quantity<K, W, U> value, std::string_view name)
{
    if (!std::isfinite(value.value()) || std::abs(value.value()) > std::numeric_limits<P>::max()) {
        throw ParameterError(std::format("{} is outside the finite parameter range", name));
    }
    return Quantity<K, P, U> {value};
}

// A quantity that is finite and positive in P, with its kind and units;
// throws ParameterError otherwise.
template <std::floating_point P, KindTag K, std::floating_point W, UnitTag U>
Quantity<K, P, U> positive_thermodynamic_parameter(Quantity<K, W, U> value, std::string_view name)
{
    const auto stored = thermodynamic_parameter<P>(value, name);
    if (!(stored.value() > P(0))) {
        throw ParameterError(std::format("{} must remain positive after conversion", name));
    }
    return stored;
}

// A finite value converted to P; throws ParameterError otherwise, naming
// what.
template <std::floating_point P>
P finite_parameter(std::common_type_t<double, P> value, std::string_view what)
{
    const auto converted = static_cast<P>(value);
    if (!std::isfinite(value) || !std::isfinite(converted)) {
        throw ParameterError(std::format("{} must be finite, got {}", what, value));
    }
    return converted;
}

// A finite, positive value converted to P; throws ParameterError otherwise,
// naming what.
template <std::floating_point P>
P positive_parameter(std::common_type_t<double, P> value, std::string_view what)
{
    const auto converted = static_cast<P>(value);
    if (!std::isfinite(value) || !(value > 0.0) || !std::isfinite(converted) ||
        !(converted > P(0))) {
        throw ParameterError(std::format("{} must be finite and positive, got {}", what, value));
    }
    return converted;
}

}  // namespace canary::detail

#endif  // CANARY_DETAIL_PARAMETERS_HPP
