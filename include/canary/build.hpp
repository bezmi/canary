#ifndef CANARY_BUILD_HPP
#define CANARY_BUILD_HPP

// Host-side support for building models. Not used in per-cell code.

#include <array>
#include <concepts>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace canary
{

// Thrown during checked host construction when a parameter is out of range.
// The message names the parameter and its value.
class ParameterError : public std::invalid_argument
{
  public:
    using std::invalid_argument::invalid_argument;
};

namespace detail
{

// Storage is specialized by each supported host input domain. Its scratch
// storage owns temporary view descriptors, never the model's numeric data.
template <typename Input> struct InputStorage;
struct InputAccess;

struct SpeciesNameTag {};
struct ReactionEquationTag {};
struct TemperatureLabelTag {};
struct ElementNameTag {};

template <typename Tag> struct OwnedText {
    std::optional<std::string> value;
    OwnedText() = default;
    explicit OwnedText(const char* text)
    {
        if (text != nullptr) { value.emplace(text); }
    }
    explicit OwnedText(std::string_view text) : value(std::string(text)) {}
    [[nodiscard]] const char* c_str() const noexcept { return value ? value->c_str() : nullptr; }
};

using OwnedSpeciesName = OwnedText<SpeciesNameTag>;
using OwnedReactionEquation = OwnedText<ReactionEquationTag>;
using OwnedTemperatureLabel = OwnedText<TemperatureLabelTag>;
using OwnedElementName = OwnedText<ElementNameTag>;

template <typename T, std::size_t N> struct InputStorage<std::array<T, N>> {
    std::array<InputStorage<T>, N> entries {};
    using Scratch = std::array<typename InputStorage<T>::Scratch, N>;
    InputStorage() = default;
    explicit InputStorage(const std::array<T, N>& input)
    {
        for (std::size_t i = 0; i < N; ++i) { entries[i] = InputStorage<T> {input[i]}; }
    }
    std::array<T, N> borrow(Scratch& scratch) const
    {
        std::array<T, N> input {};
        for (std::size_t i = 0; i < N; ++i) { input[i] = entries[i].borrow(scratch[i]); }
        return input;
    }
};

}  // namespace detail

// Independent host-owned copy of a supported construction input. All numeric
// fields must be initialized and borrowed strings/spans valid during copying.
// Null text remains null; no validation is substituted for the model builder.
// Copying is independent, and moved-to owners retain their inputs after the
// source is destroyed. A moved-from owner may be assigned or destroyed.
// Builders accept this owner directly and bind borrowed views only for their
// call. There is no borrowed conversion. Owning storage is never device data.
template <typename Input> class OwnedInput
{
  public:
    explicit OwnedInput(const Input& input) : storage_(input) {}

  private:
    detail::InputStorage<Input> storage_;
    friend struct detail::InputAccess;
};

namespace detail
{

struct InputAccess {
    template <typename Input, typename Function>
    static decltype(auto) apply(const OwnedInput<Input>& owner, Function&& function)
    {
        typename InputStorage<Input>::Scratch scratch {};
        const Input input = owner.storage_.borrow(scratch);
        return std::invoke(std::forward<Function>(function), input);
    }
};

template <typename T> inline constexpr bool is_owned_input = false;
template <typename Input> inline constexpr bool is_owned_input<OwnedInput<Input>> = true;
template <typename... Inputs>
concept HasOwnedInput = (is_owned_input<std::remove_cvref_t<Inputs>> || ...);

template <typename Input, typename Function>
decltype(auto) with_input(const Input& input, Function&& function)
{
    if constexpr (is_owned_input<Input>) {
        return InputAccess::apply(input, std::forward<Function>(function));
    } else {
        return std::invoke(std::forward<Function>(function), input);
    }
}

template <typename Function> decltype(auto) with_inputs(Function&& function)
{ return std::invoke(std::forward<Function>(function)); }

template <typename Function, typename First, typename... Rest>
decltype(auto) with_inputs(Function&& function, const First& first, const Rest&... rest)
{
    return with_input(first, [&](const auto& input) -> decltype(auto) {
        return with_inputs(
          [&](const auto&... tail) -> decltype(auto) {
              return std::invoke(std::forward<Function>(function), input, tail...);
          },
          rest...);
    });
}

}  // namespace detail

// A model together with the species data (particle masses and charges) that
// were given to build it.
template <typename Model, typename Species>
    requires std::is_trivially_copyable_v<Model> && std::is_trivially_copyable_v<Species>
struct Built {
    Model model;
    Species species;
};

}  // namespace canary

#endif  // CANARY_BUILD_HPP
