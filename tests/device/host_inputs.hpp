#ifndef CANARY_TESTS_DEVICE_HOST_INPUTS_HPP
#define CANARY_TESTS_DEVICE_HOST_INPUTS_HPP

// Host preparation of fixture inputs shared by the subsystem fixtures
// ({subsystem}_inputs.hpp).
//
// A fixture declares each state input as a binary64 value in SI units. An
// input of state scalar T with value type V is that value rounded to V,
// converted to the unit system U with the fixture scales in parameter type P,
// and stored in T. For a Dual T, each evaluation seeds at most one input: its
// derivative is one SI unit of its kind, converted the same way (one for a raw
// mass-fraction weight); every other input has derivative zero.

#include "../dual.hpp"
#include "cell_inputs.hpp"

#include <array>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace canary::test
{

// The unit scales of every fixture in parameter type P.
template <std::floating_point P, UnitSystem U> UnitScales<P, U> fixture_scales()
{
    if constexpr (std::same_as<U, si>) {
        return {};
    } else {
        return UnitScales<P, U> {{
          .mass = P(4.65e-26),
          .density = P(1.0e-3),
          .temperature = P(1000.0),
          .length = P(1.0),
        }};
    }
}

// The scalar inputs of one evaluation (see above). Inputs are numbered in the
// order they are requested; seeded selects the input that carries the seed.
template <typename T> class Seeds
{
  public:
    using value_type = ScalarTraits<T>::value_type;

    explicit Seeds(std::optional<std::size_t> seeded) : seeded_(seeded) {}

    // The input name, given as value_si in SI units, as a quantity of kind K in U.
    template <KindTag K, std::floating_point P, UnitSystem U>
    Quantity<K, T, U>
      quantity(const std::string& name, double value_si, const UnitScales<P, U>& scales)
    {
        const auto value =
          convert<U>(Quantity<K, value_type, si> {static_cast<value_type>(value_si)}, scales)
            .value();
        const auto unit = convert<U>(Quantity<K, value_type, si> {value_type(1)}, scales).value();
        return Quantity<K, T, U> {scalar(name, value, unit)};
    }

    // The values of the inputs name[i] of kind K, given in SI units, in U.
    template <KindTag K, std::size_t C, std::floating_point P, UnitSystem U>
    std::array<T, C> values(const std::string& name,
                            const std::array<double, C>& values_si,
                            const UnitScales<P, U>& scales)
    {
        std::array<T, C> out {};
        for (std::size_t i = 0; i < C; ++i) {
            out[i] =
              quantity<K>(name + "[" + std::to_string(i) + "]", values_si[i], scales).value();
        }
        return out;
    }

    // The raw mass-fraction weights name[i], each one input.
    template <std::size_t N>
    std::array<T, N> fractions(const std::string& name, const std::array<double, N>& weights)
    {
        std::array<T, N> out {};
        for (std::size_t i = 0; i < N; ++i) {
            out[i] = scalar(name + "[" + std::to_string(i) + "]",
                            static_cast<value_type>(weights[i]),
                            value_type(1));
        }
        return out;
    }

    // The number of inputs requested so far.
    [[nodiscard]] std::size_t count() const { return count_; }

    // "seed <input name>" or "no seed".
    [[nodiscard]] std::string label() const
    { return seeded_name_.empty() ? std::string("no seed") : "seed " + seeded_name_; }

  private:
    T scalar(const std::string& name, value_type value, value_type unit)
    {
        const bool seeded = seeded_ == count_;
        ++count_;
        if constexpr (std::floating_point<T>) {
            return value;
        } else {
            if (seeded) { seeded_name_ = name; }
            return T(value, seeded ? unit : value_type(0));
        }
    }

    std::optional<std::size_t> seeded_;
    std::size_t count_ = 0;
    std::string seeded_name_;
};

// Calls evaluate(seeds) once without a seed and, for a Dual T, once more for
// each input seeded in turn.
template <typename T, typename F> void for_each_seed(F&& evaluate)
{
    Seeds<T> unseeded(std::nullopt);
    evaluate(unseeded);
    if constexpr (!std::floating_point<T>) {
        for (std::size_t seeded = 0; seeded < unseeded.count(); ++seeded) {
            Seeds<T> seeds(seeded);
            evaluate(seeds);
        }
    }
}

// Storage of the C component values of one cell as storage S reads it (see
// cell_inputs.hpp): the values themselves for Contiguous; otherwise the middle
// row of a [3, C] grid after one padding element, with the padding, the
// other rows and one trailing element holding 12345.
template <Storage S, typename T, std::size_t C>
std::vector<T> cell_storage(const std::array<T, C>& values)
{
    if constexpr (S == Storage::Contiguous) {
        return std::vector<T>(values.begin(), values.end());
    } else {
        using Layout =
          std::conditional_t<S == Storage::LayoutRight, Kokkos::layout_right, Kokkos::layout_left>;
        std::vector<T> out((3 * C) + 2, T(12345));
        const typename Layout::template mapping<Kokkos::extents<std::size_t, 3, C>> grid {};
        for (std::size_t c = 0; c < C; ++c) { out[1 + grid(1, c)] = values[c]; }
        return out;
    }
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_HOST_INPUTS_HPP
