#ifndef CANARY_TESTS_DEVICE_NATIVE_RUNTIME_HPP
#define CANARY_TESTS_DEVICE_NATIVE_RUNTIME_HPP

// Evaluation of one check entry in the memory of the selected backend. The
// including unit first includes one backend header (host_backend.hpp,
// cuda_backend.hpp or hip_backend.hpp), which defines canary::test::backend.
//
// run<entry>(args...) copies each argument for a pointer parameter into backend
// storage, passes every other argument by value, launches the entry once,
// copies its complete result back and releases the storage. Arguments are
// host values only; views over backend storage are built inside the entry.
// run_on_host<entry>(args...) calls the same entry as an ordinary host
// function on the host values.

#include "cell_inputs.hpp"
#include "runtime_checks.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <ranges>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>

namespace canary::test
{

// Storage for count values of V in backend memory. Values are copied in and
// out as bytes, so V must be trivially copyable.
template <typename V> class BackendArray
{
    static_assert(std::is_trivially_copyable_v<V>);

  public:
    explicit BackendArray(std::size_t count) :
        data_(static_cast<V*>(backend::allocate_bytes(count * sizeof(V))))
    {
    }

    explicit BackendArray(std::span<const V> values) : BackendArray(values.size())
    { backend::copy_to_device(data_, values.data(), values.size_bytes()); }

    BackendArray(const BackendArray&) = delete;
    BackendArray& operator=(const BackendArray&) = delete;
    ~BackendArray() { backend::release(data_); }

    [[nodiscard]] V* data() const { return data_; }

    // The first value, copied to the host.
    [[nodiscard]] V front() const
    {
        std::array<std::byte, sizeof(V)> bytes {};
        backend::copy_to_host(bytes.data(), data_, sizeof(V));
        return std::bit_cast<V>(bytes);
    }

  private:
    V* data_;
};

namespace detail
{

// The argument passed to an entry parameter of type Param. A pointer parameter
// const V* receives backend storage holding either one V or every value of a
// contiguous range of V; any other parameter receives the argument unchanged,
// which must already have type Param.
template <typename Param> class Staged
{
  public:
    explicit Staged(const Param& value) : value_(value) {}
    [[nodiscard]] Param argument() const { return value_; }

  private:
    Param value_;
};

template <typename V> class Staged<const V*>
{
  public:
    explicit Staged(const V& value) : storage_(std::span<const V>(&value, 1)) {}

    template <std::ranges::contiguous_range R>
        requires std::same_as<std::ranges::range_value_t<R>, V>
    explicit Staged(const R& values) : storage_(std::span<const V>(values))
    {
    }

    [[nodiscard]] const V* argument() const { return storage_.data(); }

  private:
    BackendArray<V> storage_;
};

template <auto Entry, typename Result, typename... Params, typename... Args>
Result run(void (* /*entry*/)(Output<Result>, Params...), const Args&... args)
{
    static_assert(sizeof...(Params) == sizeof...(Args), "one argument per entry parameter");
    const BackendArray<Result> result(1);
    const std::tuple<Staged<Params>...> staged(args...);
    std::apply(
      [&](const auto&... argument) {
          backend::launch<Entry>(Output<Result>(result.data()), argument.argument()...);
      },
      staged);
    return result.front();
}

// The argument a host call passes to an entry parameter: for a pointer
// parameter const V* the address of the host value, or of the first value of
// a contiguous range of V; for any other parameter the argument unchanged.
template <typename V, typename Arg>
const V* host_argument(std::type_identity<const V*> /*parameter*/, const Arg& argument)
{
    if constexpr (std::same_as<Arg, V>) {
        return &argument;
    } else {
        return std::ranges::data(argument);
    }
}

template <typename Param, typename Arg>
Param host_argument(std::type_identity<Param> /*parameter*/, const Arg& argument)
{ return argument; }

template <auto Entry, typename Result, typename... Params, typename... Args>
Result run_on_host(void (* /*entry*/)(Output<Result>, Params...), const Args&... args)
{
    Result result {};
    Entry(Output<Result>(&result), host_argument(std::type_identity<Params> {}, args)...);
    return result;
}

}  // namespace detail

// The complete result of one evaluation of Entry with args on the backend
// (see above).
template <auto Entry, typename... Args> auto run(const Args&... args)
{ return detail::run<Entry>(Entry, args...); }

// The complete result of one call of Entry with args on the host, reading the
// arguments from host memory.
template <auto Entry, typename... Args> auto run_on_host(const Args&... args)
{ return detail::run_on_host<Entry>(Entry, args...); }

// Evaluates Entry with args on the backend and on the host and records both
// results in report under fixture (see check_result), with the limits of
// agreement limits. The result has no exact fields.
template <auto Entry, typename... Args>
void evaluate_within(Report& report,
                     const std::string& fixture,
                     const Agreement& limits,
                     const Args&... args)
{
    check_result(report,
                 entry_name<Entry>(),
                 fixture,
                 run<Entry>(args...),
                 run_on_host<Entry>(args...),
                 limits);
}

// As evaluate_within, with the default limits.
template <auto Entry, typename... Args>
void evaluate(Report& report, const std::string& fixture, const Args&... args)
{ evaluate_within<Entry>(report, fixture, Agreement {}, args...); }

// As evaluate, for a result with exact fields declared by expected.
template <auto Entry, typename Expected, typename... Args>
void evaluate_expecting(Report& report,
                        const std::string& fixture,
                        const Expected& expected,
                        const Args&... args)
{
    check_result(report,
                 entry_name<Entry>(),
                 fixture,
                 run<Entry>(args...),
                 run_on_host<Entry>(args...),
                 Agreement {},
                 expected);
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_NATIVE_RUNTIME_HPP
