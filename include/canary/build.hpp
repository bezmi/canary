#ifndef CANARY_BUILD_HPP
#define CANARY_BUILD_HPP

// Host-side support for building models. Not used in per-cell code.

#include <stdexcept>
#include <type_traits>

namespace canary
{

// Thrown during checked host construction when a parameter is out of range.
// The message names the parameter and its value.
class ParameterError : public std::invalid_argument
{
  public:
    using std::invalid_argument::invalid_argument;
};

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
