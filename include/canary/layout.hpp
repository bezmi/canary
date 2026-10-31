#ifndef CANARY_LAYOUT_HPP
#define CANARY_LAYOUT_HPP

// Compile-time checks that an object consists only of its members' bytes.
// A type without padding has one byte representation for each set of
// member values, so two objects with equal members are equal byte for byte
// and a byte copy of the object identifies its values.
//
// A class takes part by listing the types of all its non-static data
// members in a static member function:
//   static consteval bool padding_free()
//   { return layout::members_fill<C, decltype(a), decltype(b)>(); }
// A member left out of the list makes the sum of the sizes smaller than
// sizeof(C), so the check is false. The list holds types, not members:
// listing one member twice in place of another of the same type is not
// detected, so each member is listed once, in declaration order.
// The check applies to classes without base classes, unions, bit-fields
// and [[no_unique_address]] members; a class with any of these is not
// described by the sizes of its members.

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace canary
{

// An arithmetic type other than long double (whose storage holds unused
// bytes on common targets), or an enumeration.
template <typename T>
concept PaddingFreeScalar =
  (std::is_arithmetic_v<T> && !std::same_as<std::remove_cv_t<T>, long double>) ||
  std::is_enum_v<T>;

// A class that lists its members for the check through padding_free().
template <typename T>
concept DeclaresMembers = requires {
    { T::padding_free() } -> std::same_as<bool>;
};

namespace detail
{
template <typename T> inline constexpr bool padding_free_v = false;
template <PaddingFreeScalar T> inline constexpr bool padding_free_v<T> = true;
template <DeclaresMembers T> inline constexpr bool padding_free_v<T> = T::padding_free();
template <typename E, std::size_t N>
inline constexpr bool padding_free_v<std::array<E, N>> =
  sizeof(std::array<E, N>) == N * sizeof(E) && padding_free_v<std::remove_cv_t<E>>;
}  // namespace detail

// A type without padding bytes: a PaddingFreeScalar, a std::array of
// PaddingFree elements that holds nothing besides them, or a
// DeclaresMembers class whose padding_free() is true.
template <typename T>
concept PaddingFree = detail::padding_free_v<std::remove_cv_t<T>>;

namespace layout
{

// True when the sizes of Members, the types of all non-static data members
// of C, add up to sizeof(C) and every member type is PaddingFree.
template <typename C, typename... Members> consteval bool members_fill()
{ return sizeof(C) == (std::size_t {0} + ... + sizeof(Members)) && (PaddingFree<Members> && ...); }

}  // namespace layout

}  // namespace canary

#endif  // CANARY_LAYOUT_HPP
