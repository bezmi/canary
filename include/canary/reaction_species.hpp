// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_REACTION_SPECIES_HPP
#define CANARY_REACTION_SPECIES_HPP

#include <array>
#include <canary/species.hpp>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <string_view>
#include <utility>

namespace canary
{

// Signed number of elementary charges on a chemical species.
class ChargeNumber
{
  public:
    CANARY_HOST_DEVICE constexpr explicit ChargeNumber(std::int32_t value) noexcept : value_(value)
    {
    }
    template <std::floating_point T> ChargeNumber(T) = delete;
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::int32_t value() const noexcept
    { return value_; }
    friend CANARY_HOST_DEVICE constexpr bool operator==(ChargeNumber,
                                                        ChargeNumber) noexcept = default;

  private:
    std::int32_t value_;
};

// Host view of species definitions in reaction/model order. The definitions,
// their names and elemental metadata must outlive the view and construction.
// No thermodynamic fit is accessed. Source records provide name, molar_mass
// as MolarMass<double> or kg/mol, charge as ChargeNumber or elementary charges,
// and SpeciesMetadataView metadata.
template <typename Source, std::size_t N> class ReactionSpeciesView
{
  public:
    constexpr explicit ReactionSpeciesView(const std::array<Source, N>& input) noexcept :
        input_(&input)
    {
    }
    ReactionSpeciesView(const std::array<Source, N>&&) = delete;
    [[nodiscard]] std::string_view name(SpeciesIndex i) const
    {
        const auto& name = (*input_)[i.value()].name;
        if constexpr (std::is_pointer_v<std::remove_cvref_t<decltype(name)>>) {
            if (name == nullptr) { throw ParameterError("species name must not be null"); }
        }
        const std::string_view value {name};
        if (value.empty()) { throw ParameterError("species name must not be empty"); }
        return value;
    }
    [[nodiscard]] MolarMass<double> molar_mass(SpeciesIndex i) const
    { return MolarMass<double> {(*input_)[i.value()].molar_mass}; }
    [[nodiscard]] SpeciesMetadataView metadata(SpeciesIndex i) const
    { return (*input_)[i.value()].metadata; }
    // Only reaction participants require integral charge numbers.
    [[nodiscard]] ChargeNumber charge_number(SpeciesIndex i) const
    {
        const auto charge = (*input_)[i.value()].charge;
        if constexpr (std::same_as<std::remove_cvref_t<decltype(charge)>, ChargeNumber>) {
            return charge;
        } else {
            const double value = charge;
            if (!std::isfinite(value) || value != std::floor(value) ||
                value < std::numeric_limits<std::int32_t>::min() ||
                value > std::numeric_limits<std::int32_t>::max()) {
                throw ParameterError(std::format(
                  R"(species "{}": charge number must be an integer in the signed 32-bit range)",
                  name(i)));
            }
            return ChargeNumber {static_cast<std::int32_t>(value)};
        }
    }

  private:
    const std::array<Source, N>* input_;
};

template <typename Source, std::size_t N>
ReactionSpeciesView<Source, N> reaction_species(const std::array<Source, N>& input) noexcept
{ return ReactionSpeciesView<Source, N> {input}; }

template <typename Source, std::size_t N>
void reaction_species(const std::array<Source, N>&&) = delete;

namespace detail
{
template <typename Input> struct OwnedReactionSpecies {
    OwnedInput<Input> input;
};
}  // namespace detail

// Retain owned inputs, including temporaries, until reaction construction ends.
template <typename Source, std::size_t N>
auto reaction_species(OwnedInput<std::array<Source, N>> input)
{ return detail::OwnedReactionSpecies<std::array<Source, N>> {std::move(input)}; }

}  // namespace canary
#endif
