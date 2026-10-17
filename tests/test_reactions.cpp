#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_kim_jo.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <random>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;
using Air = canary::Air5<double, si>;
using AirScaled = canary::Air5<double, scaled>;
using Mechanism = canary::Air5KimJo<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;

constexpr std::size_t species_count = 5;
constexpr std::size_t reaction_count = 13;
using Fractions = std::array<double, species_count>;
using PerReaction = std::array<double, reaction_count>;

constexpr double atmosphere = 101325.0;
constexpr double avogadro = canary::constants::avogadro<double>;
// Universal gas constant of the reference implementation [J/(mol K)].
constexpr double reference_gas_constant = 8.31451;

const canary::UnitScales<double, scaled> reference_scales {{
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
}};

template <typename T> auto fractions(const std::array<T, species_count>& y)
{
    return canary::as_kind<canary::kind::mass_fraction>(
      Kokkos::mdspan<const T, Extents>(y.data(), species_count));
}

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

// A state of the reference implementation of the mechanism: its rates of
// progress in mol/(m^3 s) and its mass production rates in kg/(m^3 s),
// computed with R_u = 8.31451 J/(mol K) and the standard pressure 1 atm.
struct ReferenceState {
    Fractions y;
    double temp;
    double rho;
    PerReaction forward;
    PerReaction backward;
    Fractions source;
};

// clang-format off
const std::array<ReferenceState, 5> reference_states {
 {
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .temp = 5000.0, .rho = 0.006939829108798786,
    .forward = {0.02086957258136873, 0.01193076483970515, 0.0, 0.0, 0.0, 21.452051125326943, 53.11942775274426, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .backward = {},
    .source = {-0.0009188489723115111, -2.386197838323625, 0.0009188489723115111, 2.386197838323625, 0.0},
  },
  {
    .y = {0.6, 0.1, 0.1, 0.15, 0.05}, .temp = 8000.0, .rho = 0.007037619261282653,
    .forward = {36.11917338782714, 9.449099298657492, 5.0383061883730536, 31.919883547238125, 121.48841955416773, 564.9614172603348, 258.8716541637941, 43.95353546523503, 9092.468950324308, 3582.0073365831445, 362.6249273227603, 571.4791326314848, 6547.884340817446},
    .backward = {0.0465333278431809, 0.012173535389804535, 0.006490988902777724, 0.04112326685528125, 0.15651688358320112, 0.014252179409990634, 0.006530508361429085, 0.001108807883185507, 0.22937406838895197, 0.09036265070356325, 0.34483357978179385, 21542.18651615003, 8871.057891464243},
    .source = {59.372202435036435, -1104.3627016168407, -315.4882659828913, 811.8096129430239, 548.6691522216718},
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .temp = 12000.0, .rho = 0.0009537284289470741,
    .forward = {11.3324353641388, 5.070311373421975, 2.703514944892124, 59.79775804381441, 108.64952943047089, 55.42078097247583, 19.88725849679765, 8.623382721593023, 3253.5127831374225, 712.0732261767776, 165.10652060828053, 39.337479941326585, 330.90062474954647},
    .backward = {0.00028215932058449665, 0.00012624255654735965, 6.731315163329544e-05, 0.0014888675064065548, 0.0027051976403005633, 5.7203579795290745e-05, 2.0527000131884224e-05, 8.90078329760392e-06, 0.0033581731408229513, 0.0007349794950367644, 0.004831933515966059, 1115.3539068743562, 281.5217883242919},
    .source = {-6.637150873842631, -164.01079908131675, -6.813393138677023, 148.6466780625252, 28.81466503131116},
  },
  // Below the lower temperature limit: rates at 300 K.
  {
    .y = {0.6, 0.1, 0.1, 0.15, 0.05}, .temp = 250.0, .rho = 1.1260190818052245,
    .forward = {9.302211786680948e-151, 8.642965192338661e-151, 4.608471520734223e-151, 8.551087362336059e-151, 1.1112383818721133e-149, 4.405252399611495e-76, 4.0035624018478834e-73, 3.4272502805339773e-77, 3.98903507851722e-73, 1.5714931769578992e-73, 1.5378558164223767e-98, 1.1597380101060866e-20, 1.397106543456571e-44},
    .backward = {24152302.333758928, 22440631.666159894, 11965455.100114929, 22202079.67671047, 288522407.1363414, 186304.83039025188, 169316748.86525002, 14494.363189635562, 168702366.24567926, 66460838.842832245, 4497691.490702478, 8020.115369279894, 1240303008.3030038},
    .source = {45089973.208899, 12949041.84974347, -27780531.21325408, 6822973.419029424, -37081457.26441781},
  },
  // Above the upper temperature limit: rates at 30000 K.
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .temp = 40000.0, .rho = 0.00028611852868412217,
    .forward = {96.28826408351549, 30.24683197649958, 16.12775947306739, 502.52737869697563, 648.1463994964196, 76.11369194089052, 3.2822268167112036, 11.843165766388294, 1451.145961431145, 317.60200598112453, 638.195866106026, 9.362685332361984, 81.2541390796447},
    .backward = {2.371158872684948e-07, 7.44847201207008e-08, 3.9715618860802984e-08, 1.2375051768831707e-06, 1.5961011454435885e-06, 2.1313739978806003e-07, 9.191057106673343e-09, 3.316383021687169e-08, 4.063572125401226e-06, 8.893651588318524e-07, 1.6100815331684422e-05, 84.35239825157483, 7.931200106984255},
    .source = {-38.28478118462953, -61.91693433659379, 45.14642817564361, 69.75477157648831, -14.699484230908613},
  },
 },
};
// clang-format on

// Change of the coefficient of species i in reaction k.
int coefficient_change(const Mechanism& mech, std::size_t k, std::size_t i)
{
    return mech.reactions()[canary::ReactionIndex {k}].products()[canary::SpeciesIndex {i}] -
           mech.reactions()[canary::ReactionIndex {k}].reactants()[canary::SpeciesIndex {i}];
}

struct Rates {
    PerReaction forward;
    PerReaction backward;
    Fractions production;
};

template <typename M, typename G>
Rates evaluate(const M& mech, const G& gas, const Fractions& y, double rho, double temp)
{
    using U = M::units;
    const auto mix = canary::mixture(gas, fractions(y));
    const canary::Density<double, U> rho_q {rho};
    const canary::Temperature<double, U> temp_q {temp};
    const auto q = canary::rates_of_progress(mech, gas, mix, rho_q, temp_q);
    const auto omega = canary::production_rates(mech, gas, mix, rho_q, temp_q);
    Rates out {};
    for (const auto k : canary::reaction_indices(0, reaction_count)) {
        out.forward[k.value()] = q.forward()[k].value();
        out.backward[k.value()] = q.backward()[k].value();
    }
    for (std::size_t i = 0; i < species_count; ++i) { out.production[i] = omega[i].value(); }
    return out;
}

// Equilibrium constant of reaction k per particle in SI, from the
// definition: exp(-sum dnu_i g_i m_i / (k_B T)) (p0 / (k_B T))^sum(dnu).
double
  equilibrium_constant(const Mechanism& mech, const Air& gas, std::size_t k, double temp, double p0)
{
    const double kb = canary::constants::boltzmann<double>;
    const auto sp = canary::species_properties(gas, canary::Temperature<double, si> {temp});
    double exponent = 0;
    int change = 0;
    for (std::size_t i = 0; i < species_count; ++i) {
        const int dnu = coefficient_change(mech, k, i);
        exponent -= dnu * sp.gibbs_energy[i].value() *
                    mech.mass()[canary::SpeciesIndex {i}].value() / (kb * temp);
        change += dnu;
    }
    return std::exp(exponent) * std::pow(p0 / (kb * temp), change);
}

template <typename F> std::string error_message(const F& f)
{
    try {
        f();
    } catch (const canary::ParameterError& error) {
        return error.what();
    }
    return "no exception";
}

}  // namespace

static_assert(canary::HasSpeciesProperties<Air, double, 5>);
static_assert(canary::HasSpeciesProperties<Air, float, 5>);
static_assert(canary::HasSpeciesProperties<Air, D, 5>);
static_assert(canary::HasSpeciesProperties<AirScaled, D, 5>);
static_assert(!canary::HasSpeciesProperties<Air, double, 4>);
static_assert(!canary::HasSpeciesProperties<canary::CaloricallyPerfect<double, si>, double, 1>);
static_assert(canary::ReactingModel<Air, double, double, 5>);
static_assert(canary::ReactingModel<Air, float, double, 5>);
static_assert(canary::ReactingModel<Air, D, double, 5>);
static_assert(std::is_trivially_copyable_v<Mechanism>);
static_assert(std::is_trivially_copyable_v<canary::Air5KimJo<float, scaled>>);
static_assert(std::is_trivially_copyable_v<canary::RatesOfProgress<D, si, 13>>);
static_assert(canary::air5_kim_jo.reactions.size() == 13);

static_assert(Mechanism::reaction_count == reaction_count);
static_assert(canary::RatesOfProgress<D, si, 13>::reaction_count == reaction_count);
static_assert(sizeof(canary::ReactionIndex) == sizeof(std::size_t));
static_assert(std::is_trivially_copyable_v<canary::ReactionIndex>);
static_assert(std::is_standard_layout_v<canary::ReactionIndex>);
static_assert(std::is_nothrow_default_constructible_v<canary::ReactionIndex>);
static_assert(std::is_nothrow_constructible_v<canary::ReactionIndex, std::size_t>);
static_assert(!std::convertible_to<std::size_t, canary::ReactionIndex>);
static_assert(!std::convertible_to<canary::ReactionIndex, std::size_t>);
static_assert(!std::constructible_from<canary::ReactionIndex, canary::SpeciesIndex>);
static_assert(!std::constructible_from<canary::ReactionIndex, canary::ModeIndex>);
static_assert(!std::constructible_from<canary::SpeciesIndex, canary::ReactionIndex>);
static_assert(!std::constructible_from<canary::ModeIndex, canary::ReactionIndex>);
static_assert(canary::ReactionIndex {} == canary::ReactionIndex {0});
static_assert(canary::ReactionIndex {std::numeric_limits<std::size_t>::max()}.value() ==
              std::numeric_limits<std::size_t>::max());
static_assert(std::ranges::forward_range<canary::ReactionIndexRange>);
static_assert(std::ranges::sized_range<canary::ReactionIndexRange>);
static_assert(
  std::same_as<std::ranges::range_value_t<canary::ReactionIndexRange>, canary::ReactionIndex>);
static_assert(noexcept(canary::reaction_indices(0, 1)));
static_assert(noexcept(++std::declval<canary::ReactionIndexRange::iterator&>()));
static_assert(canary::reaction_indices(4, 4).size() == 0);
static_assert(canary::reaction_indices(4, 7).size() == 3);

namespace
{
template <typename Owner>
concept BorrowsReactions = requires(Owner&& owner) { std::forward<Owner>(owner).reactions(); };
template <typename Owner>
concept BorrowsReactionMass = requires(Owner&& owner) { std::forward<Owner>(owner).mass(); };
template <typename Owner>
concept BorrowsReactants = requires(Owner&& owner) { std::forward<Owner>(owner).reactants(); };
template <typename Owner>
concept BorrowsProducts = requires(Owner&& owner) { std::forward<Owner>(owner).products(); };
template <typename Owner>
concept BorrowsEfficiencies =
  requires(Owner&& owner) { std::forward<Owner>(owner).efficiencies(); };
template <typename Owner>
concept BorrowsForward = requires(Owner&& owner) { std::forward<Owner>(owner).forward(); };
template <typename Owner>
concept BorrowsBackward = requires(Owner&& owner) { std::forward<Owner>(owner).backward(); };
template <typename View, typename Index>
concept ReadsReactionView = requires(View v, Index i) { v[i]; };

template <std::floating_point P, canary::UnitSystem U, typename T>
consteval bool reaction_view_contracts()
{
    using Record = canary::Reaction<P, U, 5>;
    using Model = canary::ReactionMechanism<P, U, 5, 13>;
    using Output = canary::RatesOfProgress<T, U, 13>;
    static_assert(std::is_trivially_copyable_v<Record> && std::is_trivially_copyable_v<Model> &&
                  std::is_trivially_copyable_v<Output>);
    static_assert(std::is_copy_assignable_v<Model> && std::is_move_assignable_v<Model>);
    static_assert(!std::is_aggregate_v<Model> && !std::is_default_constructible_v<Model>);
    static_assert(BorrowsReactions<Model&> && BorrowsReactions<const Model&>);
    static_assert(!BorrowsReactions<Model> && !BorrowsReactions<const Model>);
    static_assert(BorrowsReactionMass<Model&> && BorrowsReactionMass<const Model&>);
    static_assert(!BorrowsReactionMass<Model> && !BorrowsReactionMass<const Model>);
    static_assert(BorrowsReactants<Record&> && BorrowsReactants<const Record&>);
    static_assert(!BorrowsReactants<Record> && !BorrowsReactants<const Record>);
    static_assert(BorrowsProducts<Record&> && BorrowsProducts<const Record&>);
    static_assert(!BorrowsProducts<Record> && !BorrowsProducts<const Record>);
    static_assert(BorrowsEfficiencies<Record&> && BorrowsEfficiencies<const Record&>);
    static_assert(!BorrowsEfficiencies<Record> && !BorrowsEfficiencies<const Record>);
    static_assert(BorrowsForward<Output&> && BorrowsForward<const Output&>);
    static_assert(!BorrowsForward<Output> && !BorrowsForward<const Output>);
    static_assert(BorrowsBackward<Output&> && BorrowsBackward<const Output&>);
    static_assert(!BorrowsBackward<Output> && !BorrowsBackward<const Output>);
    using Records = decltype(std::declval<Model&>().reactions());
    using Masses = decltype(std::declval<Model&>().mass());
    using Reactants = decltype(std::declval<Record&>().reactants());
    using Products = decltype(std::declval<Record&>().products());
    using Efficiencies = decltype(std::declval<Record&>().efficiencies());
    using Forward = decltype(std::declval<Output&>().forward());
    using Backward = decltype(std::declval<Output&>().backward());
    static_assert(std::same_as<typename Records::element_type, const Record>);
    static_assert(std::same_as<typename Masses::element_type, const canary::ParticleMass<P, U>>);
    static_assert(std::same_as<typename Reactants::element_type, int> &&
                  std::same_as<typename Products::element_type, int>);
    static_assert(std::same_as<typename Efficiencies::element_type, P>);
    static_assert(std::same_as<typename Forward::element_type, canary::NumberRateDensity<T, U>>);
    static_assert(std::same_as<typename Backward::element_type, canary::NumberRateDensity<T, U>>);
    static_assert(
      std::is_const_v<typename decltype(std::declval<const Record&>().reactants())::element_type>);
    static_assert(
      std::is_const_v<typename decltype(std::declval<const Record&>().products())::element_type>);
    static_assert(std::is_const_v<
                  typename decltype(std::declval<const Record&>().efficiencies())::element_type>);
    static_assert(
      std::is_const_v<typename decltype(std::declval<const Output&>().forward())::element_type>);
    static_assert(
      std::is_const_v<typename decltype(std::declval<const Output&>().backward())::element_type>);
    static_assert(Records::extents_type::static_extent(0) == 13 &&
                  Forward::extents_type::static_extent(0) == 13 &&
                  Backward::extents_type::static_extent(0) == 13);
    static_assert(Masses::extents_type::static_extent(0) == 5 &&
                  Reactants::extents_type::static_extent(0) == 5 &&
                  Products::extents_type::static_extent(0) == 5 &&
                  Efficiencies::extents_type::static_extent(0) == 5);
    static_assert(ReadsReactionView<Records, canary::ReactionIndex>);
    static_assert(!ReadsReactionView<Records, canary::SpeciesIndex> &&
                  !ReadsReactionView<Records, canary::ModeIndex> &&
                  !ReadsReactionView<Records, std::size_t>);
    static_assert(!ReadsReactionView<Masses, canary::ReactionIndex> &&
                  !ReadsReactionView<Reactants, canary::ReactionIndex> &&
                  !ReadsReactionView<Products, canary::ModeIndex> &&
                  !ReadsReactionView<Efficiencies, std::size_t>);
    static_assert(!ReadsReactionView<Forward, canary::SpeciesIndex> &&
                  !ReadsReactionView<Backward, canary::ModeIndex>);
    static_assert(noexcept(std::declval<const Model&>().reactions()) &&
                  noexcept(std::declval<const Model&>().mass()) &&
                  noexcept(std::declval<Record&>().reactants()) &&
                  noexcept(std::declval<const Record&>().reactants()) &&
                  noexcept(std::declval<Record&>().products()) &&
                  noexcept(std::declval<const Record&>().products()) &&
                  noexcept(std::declval<Record&>().efficiencies()) &&
                  noexcept(std::declval<const Record&>().efficiencies()) &&
                  noexcept(std::declval<Output&>().forward()) &&
                  noexcept(std::declval<const Output&>().forward()) &&
                  noexcept(std::declval<Output&>().backward()) &&
                  noexcept(std::declval<const Output&>().backward()));
    return true;
}
static_assert(reaction_view_contracts<double, si, double>());
static_assert(reaction_view_contracts<float, scaled, float>());
static_assert(reaction_view_contracts<double, scaled, D>());
// Reference traversal retains original reaction order and independent source
// caches, permitting exact comparison with grouped execution.
namespace schedule_reference
{
using namespace canary;
using namespace canary::detail;
template <std::size_t M,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  original_rates_of_progress_at(const ReactionMechanism<P, U, N, R>& mech,
                                const G& gas,
                                const NormalisedComposition<T, U, N>& mix,
                                Density<T, U> rho,
                                T t,
                                const std::array<T, M>& other)
{
    using V = ScalarTraits<T>::value_type;
    const auto records = mech.reactions();
    const T ln_t = math::log(t);
    const T inv_t = T(V(1)) / t;
    const EquilibriumTerms<T, N> at_t = equilibrium_terms(mech, gas, t);

    // Number densities n_i.
    std::array<T, N> n {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        n[i.value()] = T(rho.value() * mix.fractions[i.value()].value() / mech.mass()[i].value());
    }

    // ln(T_o / T_l) of every other temperature and, at each other
    // temperature that has equilibrium backward rates, ln T_o, 1 / T_o and
    // the equilibrium terms.
    std::array<T, M> ln_ratio {};
    std::array<T, M> ln_other {};
    std::array<T, M> inv_other {};
    std::array<EquilibriumTerms<T, N>, M> at_other {};
    if constexpr (M > 0) {
        std::array<bool, M> used {};
        for (const ReactionIndex k : reaction_indices(0, R)) {
            const auto& r = records[k];
            if (r.backward_model == BackwardRateModel::equilibrium &&
                r.equilibrium_at_other_temperature) {
                used[(*r.backward.other_mode).value()] = true;
            }
        }
        for (const ModeIndex m : mode_indices(0, M)) {
            ln_ratio[m.value()] = math::log(T(other[m.value()] / t));
            if (!used[m.value()]) { continue; }
            ln_other[m.value()] = math::log(other[m.value()]);
            inv_other[m.value()] = T(V(1)) / other[m.value()];
            at_other[m.value()] = equilibrium_terms(mech, gas, other[m.value()]);
        }
    }

    RatesOfProgress<T, U, R> out {};
    for (const ReactionIndex k : reaction_indices(0, R)) {
        const Reaction<P, U, N>& r = records[k];
        // ln(k_f / A) at T_l, and at the rate-controlling temperature.
        const T forward_exponent = arrhenius_exponent(r.forward, ln_t, inv_t);
        T controlled_exponent = forward_exponent;
        if constexpr (M > 0) {
            if (r.forward.other_mode && r.forward.translational_exponent != P(1)) {
                controlled_exponent = arrhenius_exponent(r.forward,
                                                         ln_t,
                                                         inv_t,
                                                         ln_ratio[(*r.forward.other_mode).value()]);
            }
        }
        const T k_f = T(r.forward.A * math::exp(controlled_exponent));
        T k_b = k_f;
        if (r.backward_model == BackwardRateModel::equilibrium) {
            bool at_other_temperature = false;
            if constexpr (M > 0) { at_other_temperature = r.equilibrium_at_other_temperature; }
            T exponent = forward_exponent;
            if (at_other_temperature) {
                const ModeIndex m = (*r.backward.other_mode);
                exponent = equilibrium_exponent(
                  r,
                  arrhenius_exponent(r.forward, ln_other[m.value()], inv_other[m.value()]),
                  at_other[m.value()]);
            } else {
                exponent = equilibrium_exponent(r, forward_exponent, at_t);
            }
            k_b = T(r.forward.A * math::exp(exponent));
        } else {
            k_b = T(r.backward.A * math::exp(rate_exponent<M>(r.backward, ln_t, inv_t, ln_ratio)));
        }
        T q_f = k_f;
        T q_b = k_b;
        if (r.third_body) {
            T m = T(r.efficiencies()[SpeciesIndex {0}] * n[0]);
            for (const SpeciesIndex i : species_indices(1, N)) {
                m = m + T(r.efficiencies()[i] * n[i.value()]);
            }
            q_f = q_f * m;
            q_b = q_b * m;
        }
        for (const SpeciesIndex i : species_indices(0, N)) {
            for (int c = 0; c < r.reactants()[i]; ++c) { q_f = q_f * n[i.value()]; }
            for (int c = 0; c < r.products()[i]; ++c) { q_b = q_b * n[i.value()]; }
        }
        out.forward()[k] = NumberRateDensity<T, U> {q_f};
        out.backward()[k] = NumberRateDensity<T, U> {q_b};
    }
    return out;
}

}  // namespace schedule_reference

template <typename G> struct CountedReactionGas {
    using parameter_type = typename G::parameter_type;
    using units = typename G::units;
    static constexpr std::size_t mode_count = G::mode_count;
    G model;
    std::vector<double>* calls;
};
template <typename G, canary::Scalar T>
auto species_properties(const CountedReactionGas<G>& gas,
                        canary::Temperature<T, typename G::units> t)
{
    gas.calls->push_back(static_cast<double>(canary::ScalarTraits<T>::value(t.value())));
    return species_properties(gas.model, t);
}

auto interleaved_reaction_inputs()
{
    using namespace canary;
    MechanismInput<22, 1> input {.temperature_lower = 300.0,
                                 .temperature_upper = 30000.0,
                                 .reactions = {}};
    for (std::size_t k = 0; k < input.reactions.size(); ++k) {
        const std::size_t behavior = (k * 7) % 20;
        const bool equilibrium = behavior >= 8;
        const bool f = equilibrium ? (behavior % 4 >= 2) : (behavior >= 4);
        const bool b = !equilibrium && (behavior % 4 >= 2);
        auto& r = input.reactions[k];
        r.equation = "N2 <=> N + N";
        r.reactants = {.count = 1, .entries = {{{.name = "N2", .value = 1.0}}}};
        r.products = {.count = 1, .entries = {{{.name = "N", .value = 2.0}}}};
        r.forward = {.A = 1e7 + 1e5 * static_cast<double>(k),
                     .n = 0.5,
                     .activation_temperature = 10000.0 + 20.0 * static_cast<double>(k),
                     .translational_exponent = 0.5,
                     .other_temperature =
                       f ? PhysicalTemperatureInput {VibrationTemperatureInput {"N2"}}
                         : PhysicalTemperatureInput {TranslationTemperatureInput {}}};
        r.third_body = behavior % 2 != 0;
        if (equilibrium) {
            r.backward_model = BackwardRateModel::equilibrium;
            r.equilibrium_temperature =
              behavior < 12   ? PhysicalTemperatureInput {TranslationTemperatureInput {}}
              : behavior < 16 ? PhysicalTemperatureInput {VibrationTemperatureInput {"N2"}}
                              : PhysicalTemperatureInput {FreeElectronTemperatureInput {}};
        } else {
            r.backward_model = BackwardRateModel::arrhenius;
            r.backward = {.A = 2e6 + 1e4 * static_cast<double>(k),
                          .n = -0.25,
                          .activation_temperature = 3000.0,
                          .translational_exponent = 0.25,
                          .other_temperature =
                            b ? PhysicalTemperatureInput {FreeElectronTemperatureInput {}}
                              : PhysicalTemperatureInput {TranslationTemperatureInput {}}};
        }
    }
    return input;
}

void check_prepared_reactions()
{
    using namespace boost::ut;
    using namespace canary;
    "prepared rates preserve outputs for interleaved behaviors and original identities"_test = [] {
        const auto input = interleaved_reaction_inputs();
        const auto check = [&]<typename P, typename U>(const UnitScales<P, U>& scales) {
            const auto built = make_air11_three_temperature(scales);
            const auto& gas = built.model;
            const auto mechanism =
              make_reaction_mechanism(input,
                                      reaction_species(air11_species),
                                      scales,
                                      convert<U>(Pressure<P, si> {P(101325)}, scales),
                                      gas);
            static_assert(std::is_trivially_copyable_v<decltype(mechanism)>);
            auto copy = mechanism;
            copy = mechanism;
            const auto evaluate = [&]<typename T>() {
                std::array<T, 11> y {};
                y[0] = T(0.5);
                y[2] = T(0.5);
                const auto mix = mixture(gas,
                                         as_kind<kind::mass_fraction>(
                                           Kokkos::mdspan<const T, Extents>(y.data(), y.size())));
                for (const double density : {1e-3, 1e-2, 1e-1, 1.0}) {
                    const auto rho = convert<U>(Density<T, si> {T(density)}, scales);
                    const T t = convert<U>(Temperature<T, si> {T(8000)}, scales).value();
                    std::array<T, 2> other {
                      convert<U>(Temperature<T, si> {T(5000)}, scales).value(),
                      convert<U>(Temperature<T, si> {T(11000)}, scales).value()};
                    if constexpr (std::same_as<T, D>) {
                        other[0] = D(t.real(), 1.0);
                        other[1] = D(t.real(), 2.0);
                    }
                    const auto same = [&](const auto& a, const auto& b) {
                        for (const auto k : reaction_indices(0, input.reactions.size())) {
                            if constexpr (std::same_as<T, D>) {
                                expect(a.forward()[k].value().real() ==
                                       b.forward()[k].value().real());
                                expect(a.forward()[k].value().dual() ==
                                       b.forward()[k].value().dual());
                                expect(a.backward()[k].value().real() ==
                                       b.backward()[k].value().real());
                                expect(a.backward()[k].value().dual() ==
                                       b.backward()[k].value().dual());
                            } else {
                                expect(a.forward()[k] == b.forward()[k]) << k.value();
                                expect(a.backward()[k] == b.backward()[k]) << k.value();
                            }
                        }
                    };
                    const std::array temperatures {ModeTemperature<T, U> {other[0]},
                                                   ModeTemperature<T, U> {other[1]}};
                    same(
                      rates_of_progress(copy, gas, mix, rho, Temperature<T, U> {t}, temperatures),
                      schedule_reference::original_rates_of_progress_at<2>(mechanism,
                                                                           gas,
                                                                           mix,
                                                                           rho,
                                                                           t,
                                                                           other));
                    same(rates_of_progress(copy, gas, mix, rho, Temperature<T, U> {t}),
                         schedule_reference::original_rates_of_progress_at<0>(mechanism,
                                                                              gas,
                                                                              mix,
                                                                              rho,
                                                                              t,
                                                                              std::array<T, 0> {}));
                }
            };
            evaluate.template operator()<double>();
            if constexpr (std::same_as<P, double>) { evaluate.template operator()<D>(); }
        };
        check(UnitScales<double, si> {});
        check(UnitScales<float, si> {});
        check(UnitScales<double, scaled> {
          {.mass = 4.65e-26, .density = 1e-3, .temperature = 1000, .length = 1}});
    };
    "prepared rates request only their equilibrium temperature sources"_test = [] {
        const UnitScales<double, si> scales {};
        const auto built = make_air11_three_temperature(scales);
        const auto& gas = built.model;
        const std::array<double, 11> y {0.5, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0};
        const auto mix = mixture(
          gas,
          as_kind<kind::mass_fraction>(Kokkos::mdspan<const double, Extents>(y.data(), y.size())));
        std::vector<double> calls;
        const CountedReactionGas counted {gas, &calls};
        auto input = interleaved_reaction_inputs();
        const auto run = [&](bool collapse) {
            const auto mechanism = make_reaction_mechanism(input,
                                                           reaction_species(air11_species),
                                                           scales,
                                                           Pressure<double, si> {101325},
                                                           gas);
            calls.clear();
            if (collapse) {
                (void)rates_of_progress(mechanism,
                                        counted,
                                        mix,
                                        Density<double, si> {0.01},
                                        Temperature<double, si> {8000});
            } else {
                const std::array other {ModeTemperature<double, si> {5000},
                                        ModeTemperature<double, si> {11000}};
                (void)rates_of_progress(mechanism,
                                        counted,
                                        mix,
                                        Density<double, si> {0.01},
                                        Temperature<double, si> {8000},
                                        other);
            }
        };
        run(false);
        expect(calls == std::vector<double> {8000, 5000, 11000});
        run(true);
        expect(calls == std::vector<double> {8000});
        for (auto& r : input.reactions) {
            if (r.backward_model == BackwardRateModel::equilibrium) {
                r.equilibrium_temperature = FreeElectronTemperatureInput {};
            }
        }
        run(false);
        expect(calls == std::vector<double> {11000});
        run(true);
        expect(calls == std::vector<double> {8000});
        for (auto& r : input.reactions) {
            r.backward_model = BackwardRateModel::arrhenius;
            r.backward = r.forward;
            r.equilibrium_temperature.reset();
        }
        run(false);
        expect(calls.empty());
        run(true);
        expect(calls.empty());
    };
    "production from stored net rates equals the reduction of both directional rates"_test = [] {
        const auto input = interleaved_reaction_inputs();
        const auto check = [&]<typename P, typename U>(const UnitScales<P, U>& scales) {
            const auto built = make_air11_three_temperature(scales);
            const auto& gas = built.model;
            const auto mechanism =
              make_reaction_mechanism(input,
                                      reaction_species(air11_species),
                                      scales,
                                      convert<U>(Pressure<P, si> {P(101325)}, scales),
                                      gas);
            const auto evaluate = [&]<typename T>() {
                std::array<T, 11> y {};
                y[0] = T(0.5);
                y[2] = T(0.5);
                const auto mix = mixture(gas,
                                         as_kind<kind::mass_fraction>(
                                           Kokkos::mdspan<const T, Extents>(y.data(), y.size())));
                // Bit patterns, so that equal non-finite results also compare equal.
                const auto bits = [](double x) { return std::bit_cast<std::uint64_t>(x); };
                const auto same = [&](const auto& a, const auto& b) {
                    for (std::size_t i = 0; i < a.size(); ++i) {
                        if constexpr (std::same_as<T, D>) {
                            expect(bits(a[i].value().real()) == bits(b[i].value().real()) &&
                                   bits(a[i].value().dual()) == bits(b[i].value().dual()))
                              << i;
                        } else {
                            expect(bits(a[i].value()) == bits(b[i].value())) << i;
                        }
                    }
                };
                for (const double density : {1e-3, 1e-2, 1e-1, 1.0}) {
                    for (const double temperature : {100.0, 8000.0, 90000.0}) {
                        const auto rho = convert<U>(Density<T, si> {T(density)}, scales);
                        const auto t = convert<U>(Temperature<T, si> {T(temperature)}, scales);
                        std::array modes {
                          ModeTemperature<T, U> {
                            convert<U>(Temperature<T, si> {T(5000)}, scales).value()},
                          ModeTemperature<T, U> {
                            convert<U>(Temperature<T, si> {T(11000)}, scales).value()}};
                        if constexpr (std::same_as<T, D>) {
                            modes[0] = ModeTemperature<T, U> {D(t.value().real(), 1.0)};
                            modes[1] = ModeTemperature<T, U> {D(t.value().real(), 2.0)};
                        }
                        const auto view = as_indexed<ModeIndex>(
                          Kokkos::mdspan<const ModeTemperature<T, U>, Extents>(modes.data(), 2));
                        same(production_rates(mechanism, gas, mix, rho, t, modes),
                             canary::detail::mass_production_rates(
                               mechanism,
                               rates_of_progress(mechanism, gas, mix, rho, t, modes)));
                        same(production_rates(mechanism, gas, mix, rho, t, view),
                             canary::detail::mass_production_rates(
                               mechanism,
                               rates_of_progress(mechanism, gas, mix, rho, t, view)));
                        same(production_rates(mechanism, gas, mix, rho, t),
                             canary::detail::mass_production_rates(
                               mechanism,
                               rates_of_progress(mechanism, gas, mix, rho, t)));
                    }
                }
            };
            evaluate.template operator()<double>();
            if constexpr (std::same_as<P, double>) { evaluate.template operator()<D>(); }
        };
        check(UnitScales<double, si> {});
        check(UnitScales<float, si> {});
        check(UnitScales<double, scaled> {
          {.mass = 4.65e-26, .density = 1e-3, .temperature = 1000, .length = 1}});
    };
}

// Net-rate storage is owned by the producing call and read through a const view.
template <typename A>
concept BorrowsNetRates =
  requires(A&& net) { canary::detail::net_rate_view(std::forward<A>(net)); };
using NetRates = std::array<canary::NumberRateDensity<double, canary::si>, 3>;
static_assert(BorrowsNetRates<const NetRates&>);
static_assert(!BorrowsNetRates<NetRates&&>);
static_assert(std::is_const_v<decltype(canary::detail::net_rate_view(
                std::declval<const NetRates&>()))::element_type>);
static_assert(
  std::same_as<canary::detail::RateOutputOf<canary::detail::RateOutput::Net, double, canary::si, 3>,
               NetRates>);
static_assert(
  std::same_as<
    canary::detail::RateOutputOf<canary::detail::RateOutput::Directional, double, canary::si, 3>,
    canary::RatesOfProgress<double, canary::si, 3>>);

}  // namespace

void test_reactions()
{
    check_prepared_reactions();
    using namespace boost::ut;
    const canary::UnitScales<double, si> si_scales {};
    const auto built = canary::make_air5<double>(si_scales);
    const Air& gas = built.model;
    const Mechanism mech =
      canary::make_air5_kim_jo(si_scales, canary::Pressure<double, si> {atmosphere});

    "reaction views retain ownership across copies and assignments"_test = [&] {
        auto copied = mech;
        const auto original_records = mech.reactions();
        const auto copied_records = copied.reactions();
        expect(original_records.base().data_handle() != copied_records.base().data_handle());
        expect(mech.mass().base().data_handle() != copied.mass().base().data_handle());
        for (const auto k : canary::reaction_indices(0, Mechanism::reaction_count)) {
            expect(&copied_records[k] == &copied_records.base()[k.value()]);
            for (const auto i : canary::species_indices(0, species_count)) {
                expect(copied_records[k].reactants()[i] == original_records[k].reactants()[i]);
                expect(copied_records[k].products()[i] == original_records[k].products()[i]);
                expect(copied_records[k].efficiencies()[i] ==
                       original_records[k].efficiencies()[i]);
            }
        }
        const canary::ReactionIndex first {0};
        const canary::SpeciesIndex nitrogen {0};
        auto detached = copied_records[first];
        const auto reactants = detached.reactants();
        const auto products = detached.products();
        const auto efficiencies = detached.efficiencies();
        reactants[nitrogen] = 7;
        products[nitrogen] = 8;
        efficiencies[nitrogen] = 9;
        expect(detached.reactants()[nitrogen] == 7);
        expect(detached.products()[nitrogen] == 8);
        expect(detached.efficiencies()[nitrogen] == 9.0);
        expect(copied_records[first].reactants()[nitrogen] ==
               original_records[first].reactants()[nitrogen]);
        auto moved = std::move(copied);
        copied = canary::make_air5_kim_jo(si_scales);
        expect(copied.standard_pressure.value() == 100000.0);
        expect(moved.standard_pressure == mech.standard_pressure);
        expect(moved.reactions().base().data_handle() != copied.reactions().base().data_handle());
        copied = moved;
        expect(copied.standard_pressure == mech.standard_pressure);
        expect(copied.reactions()[first].reactants()[nitrogen] ==
               original_records[first].reactants()[nitrogen]);
    };

    "rate views alias their owner and preserve copied outputs"_test = [] {
        const auto check = []<canary::Scalar T, canary::UnitSystem U>() {
            using Rate = canary::NumberRateDensity<T, U>;
            canary::RatesOfProgress<T, U, 3> rates {};
            const auto forward = rates.forward();
            const auto backward = rates.backward();
            for (const auto k : canary::reaction_indices(0, rates.reaction_count)) {
                expect(canary::ScalarTraits<T>::value(forward[k].value()) == 0);
                expect(canary::ScalarTraits<T>::value(backward[k].value()) == 0);
                forward[k] = Rate {T(k.value() + 1)};
                backward[k] = Rate {T(k.value() + 4)};
                expect(&forward[k] == &forward.base()[k.value()]);
                expect(&backward[k] == &backward.base()[k.value()]);
            }
            auto copied = rates;
            const canary::ReactionIndex last {2};
            copied.forward()[last] = Rate {T(10)};
            copied.backward()[last] = Rate {T(20)};
            expect(canary::ScalarTraits<T>::value(rates.forward()[last].value()) == 3);
            expect(canary::ScalarTraits<T>::value(rates.backward()[last].value()) == 6);
            auto moved = std::move(copied);
            rates = moved;
            expect(canary::ScalarTraits<T>::value(rates.forward()[last].value()) == 10);
            expect(canary::ScalarTraits<T>::value(rates.backward()[last].value()) == 20);
            const auto& read_only = rates;
            expect(read_only.forward().base().data_handle() ==
                   rates.forward().base().data_handle());
            expect(read_only.backward().base().data_handle() ==
                   rates.backward().base().data_handle());
        };
        check.template operator()<double, si>();
        check.template operator()<float, scaled>();
        check.template operator()<D, si>();
    };

    "reaction indices preserve identity and half-open endpoints"_test = [] {
        std::size_t expected = 3;
        for (const auto k : canary::reaction_indices(3, 6)) { expect(k.value() == expected++); }
        expect(expected == std::size_t {6});
        const auto empty = canary::reaction_indices(6, 6);
        expect(empty.begin() == empty.end());
        const auto last = std::numeric_limits<std::size_t>::max();
        const auto range = canary::reaction_indices(last - 2, last);
        auto it = range.begin();
        const auto original = it++;
        expect((*original).value() == last - 2);
        expect((*it).value() == last - 1);
        ++it;
        expect(it == range.end());
    };

    "the builder converts rate constants to per-particle units"_test = [&] {
        // N2 + N2 <=> N + N + N2: order 2. NO + M <=> N + O + M: order 2.
        expect(relative(mech.reactions()[canary::ReactionIndex {0}].forward.A,
                        1.216e14 / avogadro,
                        1.0e-15));
        expect(mech.reactions()[canary::ReactionIndex {0}].reactants()[canary::SpeciesIndex {0}] ==
                 2_i &&
               mech.reactions()[canary::ReactionIndex {0}].products()[canary::SpeciesIndex {2}] ==
                 2_i);
        expect(mech.reactions()[canary::ReactionIndex {0}].coefficient_change == 1_i);
        expect(mech.reactions()[canary::ReactionIndex {10}].third_body);
        expect(relative(mech.reactions()[canary::ReactionIndex {10}].forward.A,
                        1.45e9 / avogadro,
                        1.0e-15));
        expect(
          mech.reactions()[canary::ReactionIndex {10}].efficiencies()[canary::SpeciesIndex {0}] ==
          1.0_d);
        expect(
          mech.reactions()[canary::ReactionIndex {10}].efficiencies()[canary::SpeciesIndex {4}] ==
          0.664827586_d);
        expect(mech.reactions()[canary::ReactionIndex {11}].coefficient_change == 0_i);
        // The translational exponent defaults to 1.
        expect(mech.reactions()[canary::ReactionIndex {0}].forward.translational_exponent == 1.0_d);
        expect(mech.temperature_lower.value() == 300.0_d);
        expect(mech.temperature_upper.value() == 30000.0_d);
        expect(mech.standard_pressure.value() == atmosphere);
    };

    "rates of progress and production rates equal the reference at five states"_test = [&] {
        // The reference values come from GDTk (gdtk-uq/gdtk commit bf559c3f,
        // libgas.so through its C interface), whose mechanism preparation
        // writes the third-body efficiency of NO, N and O in
        // NO + M <=> N + O + M to 7 significant digits (0.6648276); the file
        // keeps the source value 0.664827586.
        auto input = canary::air5_kim_jo;
        for (std::size_t e = 2; e < 5; ++e) {
            input.reactions[10].efficiencies.entries[e].value = 0.6648276;
        }
        const Mechanism reference_mech =
          canary::make_reaction_mechanism(input,
                                          canary::reaction_species(canary::air5_species),
                                          si_scales,
                                          canary::Pressure<double, si> {atmosphere});
        const double r = canary::constants::universal_gas_constant<double> / reference_gas_constant;
        for (const ReferenceState& s : reference_states) {
            const Rates got = evaluate(reference_mech, gas, s.y, s.rho, s.temp);
            PerReaction forward {};
            PerReaction backward {};
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                forward[k.value()] = avogadro * s.forward[k.value()];
                backward[k.value()] = avogadro * s.backward[k.value()] *
                                      std::pow(r, mech.reactions()[k].coefficient_change);
                expect(relative(got.forward[k.value()], forward[k.value()], 1.0e-12))
                  << s.temp << k.value() << got.forward[k.value()] << forward[k.value()];
                expect(relative(got.backward[k.value()], backward[k.value()], 1.0e-10))
                  << s.temp << k.value() << got.backward[k.value()] << backward[k.value()];
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                double expected = 0;
                double scale = 0;
                for (std::size_t k = 0; k < reaction_count; ++k) {
                    const double term = built.species.mass[i].value() *
                                        coefficient_change(mech, k, i) * (forward[k] - backward[k]);
                    expected += term;
                    scale += std::abs(term);
                }
                expect(std::abs(got.production[i] - expected) <= 1.0e-10 * scale)
                  << s.temp << i << got.production[i] << expected;
            }
        }
    };

    "production rates conserve mass"_test = [&] {
        const auto check = [&](const Fractions& y, double rho, double temp) {
            const Rates got = evaluate(mech, gas, y, rho, temp);
            double sum = 0;
            double scale = 0;
            for (const double w : got.production) {
                sum += w;
                scale += std::abs(w);
            }
            expect(std::abs(sum) <= 1.0e-13 * scale) << temp << sum << scale;
        };
        for (const ReferenceState& s : reference_states) { check(s.y, s.rho, s.temp); }
        // A fixed seed gives the same compositions in every run.
        std::mt19937 engine(12345);  // NOLINT(bugprone-random-generator-seed)
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (int n = 0; n < 200; ++n) {
            Fractions y {};
            double total = 0;
            for (double& v : y) {
                v = unit(engine);
                total += v;
            }
            for (double& v : y) { v /= total; }
            const double rho = std::pow(10.0, -4.0 + (4.0 * unit(engine)));
            const double temp = 300.0 + (29700.0 * unit(engine));
            check(y, rho, temp);
        }
    };

    "forward and backward rates balance at the equilibrium composition of each reaction"_test =
      [&] {
          for (const double temp : {3000.0, 8000.0, 15000.0}) {
              for (std::size_t k = 0; k < reaction_count; ++k) {
                  const double kc = equilibrium_constant(mech, gas, k, temp, atmosphere);
                  // Every species at n_base, except the first one whose
                  // coefficient increases, which is set so that
                  // prod n_i^dnu_i = K_c.
                  const double n_base = 1.0e22;
                  std::size_t j = species_count;
                  double others = 1.0;
                  for (std::size_t i = 0; i < species_count; ++i) {
                      const int dnu = coefficient_change(mech, k, i);
                      if (dnu > 0 && j == species_count) {
                          j = i;
                      } else {
                          others *= std::pow(n_base, dnu);
                      }
                  }
                  std::array<double, species_count> n {};
                  n.fill(n_base);
                  n[j] = std::pow(kc / others, 1.0 / coefficient_change(mech, k, j));
                  double rho = 0;
                  for (std::size_t i = 0; i < species_count; ++i) {
                      rho += n[i] * built.species.mass[i].value();
                  }
                  Fractions y {};
                  for (std::size_t i = 0; i < species_count; ++i) {
                      y[i] = n[i] * built.species.mass[i].value() / rho;
                  }
                  const Rates got = evaluate(mech, gas, y, rho, temp);
                  expect(relative(got.backward[k], got.forward[k], 1.0e-12))
                    << temp << k << got.forward[k] << got.backward[k];
              }
          }
      };

    "the default standard pressure is 1 bar"_test = [&] {
        const Mechanism bar = canary::make_air5_kim_jo(si_scales);
        expect(bar.standard_pressure.value() == 1.0e5);
        const ReferenceState& s = reference_states[1];
        const Rates a = evaluate(bar, gas, s.y, s.rho, s.temp);
        const Rates b = evaluate(mech, gas, s.y, s.rho, s.temp);
        for (const auto k : canary::reaction_indices(0, reaction_count)) {
            expect(a.forward[k.value()] == b.forward[k.value()]);
            // K_c(1 bar) / K_c(1 atm) = (q_f / q_b)(1 bar) / (q_f / q_b)(1 atm).
            const double ratio = b.backward[k.value()] / a.backward[k.value()];
            const double expected =
              std::pow(1.0e5 / atmosphere, bar.reactions()[k].coefficient_change);
            expect(relative(ratio, expected, 1.0e-13)) << k.value() << ratio << expected;
        }
    };

    "rates outside the temperature limits are those at the limits"_test = [&] {
        const canary::Density<double, si> rho {0.01};
        for (const auto& [outside, limit] :
             {std::array {250.0, 300.0}, std::array {4.0e4, 3.0e4}}) {
            const Rates a = evaluate(mech, gas, reference_states[1].y, rho.value(), outside);
            const Rates b = evaluate(mech, gas, reference_states[1].y, rho.value(), limit);
            expect(a.forward == b.forward && a.backward == b.backward &&
                   a.production == b.production)
              << outside;
            const std::array<D, species_count> y_d {0.6, 0.1, 0.1, 0.15, 0.05};
            const auto mix_d = canary::mixture(gas, fractions(y_d));
            const canary::Density<D, si> rho_d {D(rho.value())};
            const auto q = canary::rates_of_progress(mech,
                                                     gas,
                                                     mix_d,
                                                     rho_d,
                                                     canary::Temperature<D, si> {D(outside, 1.0)});
            const auto omega =
              canary::production_rates(mech,
                                       gas,
                                       mix_d,
                                       rho_d,
                                       canary::Temperature<D, si> {D(outside, 1.0)});
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                expect(q.forward()[k].value().dual() == 0.0 &&
                       q.backward()[k].value().dual() == 0.0);
                expect(q.forward()[k].value().real() == a.forward[k.value()]);
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(omega[i].value().dual() == 0.0);
            }
        }
    };

    "dual-number derivatives equal central differences"_test = [&] {
        for (const std::size_t state : {std::size_t {1}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            Fractions rho_i {};
            for (std::size_t j = 0; j < species_count; ++j) { rho_i[j] = s.rho * s.y[j]; }
            // omega from partial densities and T, with value type T.
            const auto omega = [&]<typename T>(const std::array<T, species_count>& partial,
                                               T temp) {
                T rho = partial[0];
                for (std::size_t j = 1; j < species_count; ++j) { rho = rho + partial[j]; }
                const auto view = canary::as_kind<canary::kind::density, si>(
                  Kokkos::mdspan<const T, Extents>(partial.data(), species_count));
                const canary::Density<T, si> rho_q {rho};
                const auto mix = canary::mixture(gas, view, rho_q);
                return canary::production_rates(mech,
                                                gas,
                                                mix,
                                                rho_q,
                                                canary::Temperature<T, si> {temp});
            };
            // Columns 0 to 4: d/d(rho_j); column 5: d/dT.
            for (std::size_t c = 0; c <= species_count; ++c) {
                std::array<D, species_count> partial_d {};
                for (std::size_t j = 0; j < species_count; ++j) {
                    partial_d[j] = D(rho_i[j], j == c ? 1.0 : 0.0);
                }
                const auto ad = omega(partial_d, D(s.temp, c == species_count ? 1.0 : 0.0));
                const double h = 1.0e-5 * (c == species_count ? s.temp : rho_i[c]);
                Fractions up = rho_i;
                Fractions down = rho_i;
                double t_up = s.temp;
                double t_down = s.temp;
                if (c == species_count) {
                    t_up += h;
                    t_down -= h;
                } else {
                    up[c] += h;
                    down[c] -= h;
                }
                const auto plus = omega(up, t_up);
                const auto minus = omega(down, t_down);
                for (std::size_t i = 0; i < species_count; ++i) {
                    const double fd = (plus[i].value() - minus[i].value()) / (2.0 * h);
                    const double derivative = ad[i].value().dual();
                    expect(relative(fd, derivative, 1.0e-6)) << state << c << i << derivative << fd;
                }
            }
        }
    };

    "a scaled mechanism gives the SI rates after conversion"_test = [&] {
        const auto built_s = canary::make_air5<double>(reference_scales);
        const auto mech_s = canary::make_air5_kim_jo(
          reference_scales,
          canary::convert<scaled>(canary::Pressure<double, si> {atmosphere}, reference_scales));
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto mix_s = canary::mixture(built_s.model, fractions(s.y));
            const canary::Density<double, si> rho {s.rho};
            const canary::Temperature<double, si> temp {s.temp};
            const auto rho_s = canary::convert<scaled>(rho, reference_scales);
            const auto temp_s = canary::convert<scaled>(temp, reference_scales);
            const auto q = canary::rates_of_progress(mech, gas, mix, rho, temp);
            const auto q_s = canary::rates_of_progress(mech_s, built_s.model, mix_s, rho_s, temp_s);
            // The Gibbs energies of the scaled thermo model differ from the
            // SI ones by up to 6e-13 relative at 30000 K; in K_c they are
            // multiplied by m / (k_B T), about 40.
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                const double f = canary::convert<si>(q_s.forward()[k], reference_scales).value();
                const double b = canary::convert<si>(q_s.backward()[k], reference_scales).value();
                expect(relative(f, q.forward()[k].value(), 1.0e-12)) << s.temp << k.value();
                expect(relative(b, q.backward()[k].value(), 1.0e-10)) << s.temp << k.value();
            }
            const auto omega = canary::production_rates(mech, gas, mix, rho, temp);
            const auto omega_s =
              canary::production_rates(mech_s, built_s.model, mix_s, rho_s, temp_s);
            double scale = 0;
            for (const auto w : omega) { scale = std::max(scale, std::abs(w.value())); }
            for (std::size_t i = 0; i < species_count; ++i) {
                const double w = canary::convert<si>(omega_s[i], reference_scales).value();
                expect(std::abs(w - omega[i].value()) <= 1.0e-12 * scale)
                  << s.temp << i << w << omega[i].value();
            }
        }
    };

    "float states give the double rates within float precision"_test = [&] {
        for (const std::size_t state : {std::size_t {1}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            std::array<float, species_count> y {};
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = static_cast<float>(s.y[i]); }
            const auto mix = canary::mixture(gas, fractions(y));
            const auto omega =
              canary::production_rates(mech,
                                       gas,
                                       mix,
                                       canary::Density<float, si> {static_cast<float>(s.rho)},
                                       canary::Temperature<float, si> {static_cast<float>(s.temp)});
            const Rates ref = evaluate(mech, gas, s.y, s.rho, s.temp);
            double scale = 0;
            for (const double w : ref.production) { scale = std::max(scale, std::abs(w)); }
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(std::abs(omega[i].value() - ref.production[i]) <= 1.0e-4 * scale)
                  << state << i << omega[i].value() << ref.production[i];
            }
        }
    };

    "a float-parameter scaled mechanism gives the SI rates within float precision"_test = [&] {
        const canary::UnitScales<float, scaled> scales_f {{
          .mass = 4.65e-26F,
          .density = 1.0e-3F,
          .temperature = 1000.0F,
          .length = 1.0F,
        }};
        const auto built_f = canary::make_air5<float>(scales_f);
        const auto mech_f = canary::make_air5_kim_jo(
          scales_f,
          canary::convert<scaled>(canary::Pressure<float, si> {101325.0F}, scales_f));
        for (const std::size_t state : {std::size_t {1}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            std::array<float, species_count> y {};
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = static_cast<float>(s.y[i]); }
            const auto omega = canary::production_rates(
              mech_f,
              built_f.model,
              canary::mixture(built_f.model, fractions(y)),
              canary::convert<scaled>(canary::Density<float, si> {static_cast<float>(s.rho)},
                                      scales_f),
              canary::convert<scaled>(canary::Temperature<float, si> {static_cast<float>(s.temp)},
                                      scales_f));
            const Rates ref = evaluate(mech, gas, s.y, s.rho, s.temp);
            double scale = 0;
            for (const double w : ref.production) { scale = std::max(scale, std::abs(w)); }
            for (std::size_t i = 0; i < species_count; ++i) {
                const double w = canary::convert<si>(omega[i], scales_f).value();
                expect(std::abs(w - ref.production[i]) <= 1.0e-3 * scale)
                  << state << i << w << ref.production[i];
            }
        }
    };

    "stored reaction parameters reject narrowing overflow"_test = [] {
        const canary::UnitScales<float, si> scales;
        const std::array species {canary::air5_species[2]};
        canary::MechanismInput<1, 1> base {};
        base.temperature_lower = 300;
        base.temperature_upper = 30000;
        auto& reaction = base.reactions[0];
        reaction.equation = "N <=> N";
        reaction.reactants = {.count = 1, .entries = {{{"N", 1.0}}}};
        reaction.products = reaction.reactants;
        reaction.forward = {.A = 1.0, .n = 0.0, .activation_temperature = 0.0,
                            .translational_exponent = 1.0};
        const auto build = [&](const auto& input) {
            return canary::make_reaction_mechanism(input,
                                                   canary::reaction_species(species),
                                                   scales);
        };
        const double large = 2.0 * std::numeric_limits<float>::max();
        for (const bool backward : {false, true}) {
            const std::string role = backward ? "backward" : "forward";
            const std::array parameters {
              std::pair {&canary::ArrheniusInput::A, "A ="},
              std::pair {&canary::ArrheniusInput::n, "n is"},
              std::pair {&canary::ArrheniusInput::activation_temperature, "activation temperature"}};
            for (const auto& [parameter, diagnostic] : parameters) {
                auto input = base;
                input.reactions[0].backward_model = canary::BackwardRateModel::arrhenius;
                input.reactions[0].backward = input.reactions[0].forward;
                auto& rate = backward ? input.reactions[0].backward : input.reactions[0].forward;
                rate.*parameter = large;
                const auto message = error_message([&] { (void)build(input); });
                expect(message.starts_with("reaction 0 (N <=> N): " + role + " rate:")) << message;
                expect(message.find(diagnostic) != std::string::npos) << message;
            }
        }
        auto input = base;
        input.reactions[0].third_body = true;
        input.reactions[0].forward.A = canary::constants::avogadro<double>;
        input.reactions[0].efficiencies = {.count = 1, .entries = {{{"N", large}}}};
        expect(error_message([&] { (void)build(input); }).starts_with(
          "reaction 0 (N <=> N): third-body efficiency of \"N\""));
        input.reactions[0].efficiencies.entries[0].value = 0.0;
        const auto zero = build(input);
        expect(
          zero.reactions()[canary::ReactionIndex {0}].efficiencies()[canary::SpeciesIndex {0}] ==
          0.0F);
        expect(zero.reactions()[canary::ReactionIndex {0}].forward.activation_temperature.value() ==
               0.0F);
        for (const double a : {1.0e-100, std::numeric_limits<double>::max()}) {
            input = base;
            input.reactions[0].forward.A = a;
            expect(error_message([&] { (void)build(input); }).find("A =") != std::string::npos);
        }
        input = base;
        input.temperature_upper = std::nextafter(input.temperature_lower,
                                                std::numeric_limits<double>::infinity());
        expect(error_message([&] { (void)build(input); }).find("remain ordered") != std::string::npos);
        input = base;
        input.temperature_lower = 1.0e-100;
        expect(error_message([&] { (void)build(input); }).starts_with("temperature lower limit"));
        input = base;
        input.temperature_upper = large;
        expect(error_message([&] { (void)build(input); }).starts_with("temperature upper limit"));
    };

    "stored reaction temperatures are converted before narrowing"_test = [] {
        const auto check = [&]<typename P>(P temperature_scale, double lower) {
            const canary::UnitScales<P, scaled> scales {{
              .mass = P(1), .density = P(1), .temperature = temperature_scale, .length = P(1)}};
            const std::array species {canary::air5_species[2]};
            canary::MechanismInput<1, 1> input {};
            input.temperature_lower = lower;
            input.temperature_upper = 2 * lower;
            auto& reaction = input.reactions[0];
            reaction.equation = "N <=> N";
            reaction.reactants = {.count = 1, .entries = {{{"N", 1.0}}}};
            reaction.products = reaction.reactants;
            reaction.forward = {.A = 1.0, .n = 0.0, .activation_temperature = lower,
                                .translational_exponent = 1.0};
            const auto result = canary::make_reaction_mechanism(input,
                                                                canary::reaction_species(species),
                                                                scales,
                                                                canary::Pressure<P, scaled> {P(1)});
            const P expected = static_cast<P>(lower / static_cast<double>(temperature_scale));
            expect(std::isfinite(result.temperature_lower.value()));
            expect(result.temperature_lower.value() > P(0));
            expect(result.temperature_lower.value() == expected);
            expect(result.temperature_upper.value() == P(2) * expected);
            expect(result.reactions()[canary::ReactionIndex {0}]
                     .forward.activation_temperature.value() == expected);
            input.reactions[0].forward.activation_temperature = std::numeric_limits<double>::max();
            expect(error_message([&] {
                       (void)canary::make_reaction_mechanism(input,
                                                             canary::reaction_species(species),
                                                             scales,
                                                             canary::Pressure<P, scaled> {P(1)});
                   }).starts_with("reaction 0 (N <=> N): forward rate: activation temperature"));
        };
        check(1.0e10F, 1.0e40);
        check(1.0e-10F, 1.0e-46);
        check(1.0e-10, 300.0);
    };

    "stored reaction masses must remain finite and positive"_test = [] {
        const auto check = []<typename P, typename U>(const canary::UnitScales<P, U>& scales) {
            std::array species {canary::air5_species[2]};
            canary::MechanismInput<1, 1> input {};
            input.temperature_lower = 300;
            input.temperature_upper = 30000;
            auto& reaction = input.reactions[0];
            reaction.equation = "N <=> N";
            reaction.reactants = {.count = 1, .entries = {{{"N", 1.0}}}};
            reaction.products = reaction.reactants;
            reaction.forward = {.A = 1.0, .n = 0.0, .activation_temperature = 0.0,
                                .translational_exponent = 1.0};
            for (const P mass : {P(0), P(-1), std::numeric_limits<P>::infinity(),
                                 std::numeric_limits<P>::quiet_NaN()}) {
                species[0].molar_mass = static_cast<double>(mass);
                expect(error_message([&] {
                           (void)canary::make_reaction_mechanism(input,
                                                                 canary::reaction_species(species),
                                                                 scales);
                       }).starts_with("species 0: molar mass"));
            }
            species[0].molar_mass = 0.014;
            const auto result =
              canary::make_reaction_mechanism(input, canary::reaction_species(species), scales);
            expect(result.mass()[canary::SpeciesIndex {0}] ==
                   canary::detail::particle_mass_parameter(
                     canary::MolarMass<double> {species[0].molar_mass},
                     scales));
        };
        check(canary::UnitScales<float, si> {});
        check(canary::UnitScales<double, si> {});
        check(canary::UnitScales<float, scaled> {{
          .mass = 1, .density = 1, .temperature = 1000, .length = 1}});
        check(canary::UnitScales<double, scaled> {{
          .mass = 1, .density = 1, .temperature = 1000, .length = 1}});
    };

    "reaction orders and dimension exponents are checked before narrowing"_test = [&] {
        const std::intmax_t maximum = std::numeric_limits<int>::max();
        expect(canary::detail::reaction_order(std::array {std::numeric_limits<int>::max(), 1},
                                              true,
                                              "reactant") == maximum + 2);
        const auto limit = maximum / 3 + 1;
        const auto dimension = canary::detail::rate_dimension(limit, "forward");
        expect(dimension.length == 3 * (limit - 1));
        expect(dimension.count == -(limit - 1));
        for (const auto order : {std::intmax_t(0), limit + 1, maximum}) {
            expect(error_message([&] {
                       (void)canary::detail::rate_dimension(order, "forward");
                   }).starts_with("forward rate: reaction order"));
        }
        const auto build = [&](const auto& input) {
            (void)canary::make_reaction_mechanism(input,
                                                  canary::reaction_species(canary::air5_species),
                                                  si_scales);
        };
        auto input = canary::air5_kim_jo;
        auto& reaction = input.reactions[0];
        reaction.reactants = {.count = 2, .entries = {{{"N2", double(maximum)}, {"O2", 1.0}}}};
        reaction.products = reaction.reactants;
        expect(error_message([&] { build(input); }).find("forward rate: reaction order") !=
               std::string::npos);
        reaction.third_body = true;
        expect(error_message([&] { build(input); }).find("forward rate: reaction order") !=
               std::string::npos);
        reaction.third_body = false;
        reaction.reactants = {.count = 1, .entries = {{{"N2", 1.0}}}};
        reaction.products = {.count = 2,
                             .entries = {{{"N", double(maximum)}, {"O", double(maximum)}}}};
        expect(error_message([&] { build(input); }).find("net stoichiometric coefficient") !=
               std::string::npos);
        std::swap(reaction.reactants, reaction.products);
        expect(error_message([&] { build(input); }).find("net stoichiometric coefficient") !=
               std::string::npos);
    };

    "the builder rejects invalid input with the reaction in the message"_test = [&] {
        const auto build = [&](const canary::MechanismInput<13, 5>& input) {
            (void)canary::make_reaction_mechanism(input,
                                                  canary::reaction_species(canary::air5_species),
                                                  si_scales);
        };
        const auto starts = [](const std::string& message, std::string_view prefix) {
            return message.starts_with(prefix);
        };
        auto input = canary::air5_kim_jo;
        input.reactions[3].reactants.entries[1].name = "Ar";
        std::string message = error_message([&] { build(input); });
        expect(starts(message, "reaction 3 (N2 + N <=> N + N + N): unknown species \"Ar\""))
          << message;

        input = canary::air5_kim_jo;
        input.reactions[5].forward.A = -1.0;
        message = error_message([&] { build(input); });
        expect(starts(message, "reaction 5 (O2 + N2 <=> O + O + N2): forward rate: A must be"))
          << message;
        input.reactions[5].forward.A = 0.0;
        expect(starts(error_message([&] { build(input); }), "reaction 5"));
        input.reactions[5].forward.A = std::numeric_limits<double>::infinity();
        expect(starts(error_message([&] { build(input); }), "reaction 5"));

        input = canary::air5_kim_jo;
        input.reactions[1].forward.activation_temperature = -1.0;
        expect(starts(error_message([&] { build(input); }), "reaction 1"));
        input = canary::air5_kim_jo;
        input.reactions[1].forward.n = std::numeric_limits<double>::quiet_NaN();
        expect(starts(error_message([&] { build(input); }), "reaction 1"));
        input = canary::air5_kim_jo;
        input.reactions[0].reactants.entries[0].value = 1.5;
        expect(starts(error_message([&] { build(input); }), "reaction 0"));
        input = canary::air5_kim_jo;
        input.reactions[0].products.count = 0;
        expect(starts(error_message([&] { build(input); }), "reaction 0"));
        input = canary::air5_kim_jo;
        input.reactions[1].reactants.entries[1].name = "N2";
        expect(starts(error_message([&] { build(input); }), "reaction 1"));
        input = canary::air5_kim_jo;
        input.reactions[10].efficiencies.entries[0].value = -1.0;
        expect(starts(error_message([&] { build(input); }), "reaction 10"));
        input = canary::air5_kim_jo;
        input.reactions[10].third_body = false;
        expect(starts(error_message([&] { build(input); }), "reaction 10"));
        input = canary::air5_kim_jo;
        input.reactions[11].backward_model = canary::BackwardRateModel::arrhenius;
        input.reactions[11].equilibrium_temperature.reset();
        expect(starts(error_message([&] { build(input); }),
                      "reaction 11 (NO + O <=> O2 + N): backward rate"));

        for (const double s : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
            input = canary::air5_kim_jo;
            input.reactions[2].forward.translational_exponent = s;
            message = error_message([&] { build(input); });
            expect(starts(message,
                          "reaction 2 (N2 + NO <=> N + N + NO): forward rate: translational "
                          "exponent must be finite and between 0 and 1"))
              << message;
            input = canary::air5_kim_jo;
            input.reactions[11].backward_model = canary::BackwardRateModel::arrhenius;
            input.reactions[11].equilibrium_temperature.reset();
            input.reactions[11].backward = input.reactions[11].forward;
            input.reactions[11].backward.translational_exponent = s;
            message = error_message([&] { build(input); });
            expect(starts(message,
                          "reaction 11 (NO + O <=> O2 + N): backward rate: translational exponent"))
              << message;
        }
        // An equilibrium backward rate ignores its Arrhenius input.
        input = canary::air5_kim_jo;
        input.reactions[11].backward.translational_exponent = 2.0;
        build(input);

        input = canary::air5_kim_jo;
        input.temperature_lower = 40000.0;
        expect(starts(error_message([&] { build(input); }), "temperature limits"));
        expect(starts(error_message([&] {
                          (void)canary::make_air5_kim_jo(si_scales,
                                                         canary::Pressure<double, si> {0.0});
                      }),
                      "standard pressure"));
        auto invalid_species = canary::air5_species;
        invalid_species[0].name = nullptr;
        expect(starts(error_message([&] {
                          (void)canary::make_reaction_mechanism(
                            canary::air5_kim_jo,
                            canary::reaction_species(invalid_species),
                            si_scales);
                      }),
                      "species 0: species name must not be null"));
        auto duplicate_species = canary::air5_species;
        duplicate_species[3].name = "N2";
        expect(starts(error_message([&] {
                          (void)canary::make_reaction_mechanism(
                            canary::air5_kim_jo,
                            canary::reaction_species(duplicate_species),
                            si_scales);
                      }),
                      "species \"N2\" is named twice"));
    };
    "reaction construction conserves elements and integer charge numbers"_test = [&] {
        static_assert(sizeof(canary::ChargeNumber) == sizeof(std::int32_t));
        static_assert(!std::is_convertible_v<canary::ChargeNumber, int>);
        static_assert(!std::is_convertible_v<int, canary::ChargeNumber>);
        static_assert(std::is_trivially_copyable_v<canary::ChargeNumber>);
        const auto view = canary::reaction_species(canary::air11_species);
        const auto source =
          *canary::detail::species_index(canary::detail::air11_kim_jo_species_names, "O");
        const auto ion =
          *canary::detail::species_index(canary::detail::air11_kim_jo_species_names, "O+");
        const auto electron =
          *canary::detail::species_index(canary::detail::air11_kim_jo_species_names, "e-");
        canary::MechanismInput<1, 11> input {.temperature_lower = 300.0,
                                             .temperature_upper = 30000.0,
                                             .reactions = {canary::air11_kim_jo.reactions[38]}};
        input.reactions[0].forward.other_temperature = std::nullopt;
        input.reactions[0].forward.translational_exponent = 1.0;
        const auto check = [&](const auto& definitions, const auto& mechanism) {
            return error_message([&] {
                (void)canary::make_reaction_mechanism(mechanism,
                                                      canary::reaction_species(definitions),
                                                      si_scales);
            });
        };
        expect(view.charge_number(electron).value() == -1);
        expect(check(canary::air11_species, input) == "no exception");
        auto changed = input;
        changed.reactions[0].products.entries[1].value = 1.0;
        const auto charge_error = check(canary::air11_species, changed);
        expect(charge_error.starts_with("reaction 0 (")) << charge_error;
        expect(charge_error.find("charge is not conserved") != std::string::npos) << charge_error;

        auto species = canary::air11_species;
        for (const double charge : {0.5,
                                    -0.5,
                                    2147483648.0,
                                    -2147483649.0,
                                    std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::quiet_NaN()}) {
            species[source.value()].charge = charge;
            const auto message = check(species, input);
            expect(message.find("charge number must be an integer") != std::string::npos)
              << message;
        }
        species = canary::air11_species;
        species[ion.value()].charge = 2;
        changed = input;
        changed.reactions[0].products.entries[1].value = 3.0;
        expect(check(species, changed) == "no exception");
        species[ion.value()].charge = -1;
        changed.reactions[0].products.count = 1;
        expect(check(species, changed) == "no exception");

        species = canary::air11_species;
        species[source.value()].metadata.elements.reset();
        expect(check(species, input).find("elemental composition is required") !=
               std::string::npos);
        species = canary::air11_species;
        species[source.value()].metadata.role = canary::SpeciesRole::unspecified;
        expect(check(species, input).find("species role must be specified") != std::string::npos);
        species = canary::air11_species;
        // Unused species do not participate in conservation checks.
        species[0].metadata.elements.reset();
        species[0].charge = 0.5;
        expect(check(species, input) == "no exception");

        auto dissociation = canary::air5_kim_jo;
        dissociation.reactions[0].products.entries[0].value = 1;
        const auto element_error = check(canary::air5_species, dissociation);
        expect(element_error.starts_with("reaction 0 (")) << element_error;
        expect(element_error.find("element N is not conserved") != std::string::npos)
          << element_error;
        const std::array<canary::ElementCount, 2> duplicate {{{"N", 1}, {"N", 1}}};
        auto nitrogen = canary::air5_species;
        nitrogen[0].metadata.elements = duplicate;
        expect(check(nitrogen, canary::air5_kim_jo).find("element N: given twice") !=
               std::string::npos);
        nitrogen = canary::air5_species;
        const std::array<canary::ElementCount, 1> zero {{{"N", 0}}};
        nitrogen[0].metadata.elements = zero;
        expect(check(nitrogen, canary::air5_kim_jo).find("count must be positive") !=
               std::string::npos);
    };

    "reaction species support typed definitions without thermodynamic fits"_test = [&] {
        struct Definition {
            std::string_view name;
            canary::MolarMass<double> molar_mass;
            canary::ChargeNumber charge;
            canary::SpeciesMetadataView metadata;
        };
        static_assert(!std::is_constructible_v<canary::ChargeNumber, double>);
        const std::array<canary::ElementCount, 1> molecule {{{"X", 2}}};
        const std::array<canary::ElementCount, 1> atom {{{"X", 1}}};
        const std::array definitions {Definition {"molecule",
                                                  canary::MolarMass<double> {0.028},
                                                  canary::ChargeNumber {0},
                                                  {canary::SpeciesRole::heavy, molecule}},
                                      Definition {"atom",
                                                  canary::MolarMass<double> {0.014},
                                                  canary::ChargeNumber {0},
                                                  {canary::SpeciesRole::heavy, atom}}};
        using View = decltype(canary::reaction_species(definitions));
        static_assert(!std::is_constructible_v<View, decltype(definitions)&&>);
        canary::MechanismInput<1, 2> input {};
        input.temperature_lower = 300;
        input.temperature_upper = 30000;
        input.reactions[0].equation = "molecule -> 2 atom";
        input.reactions[0].reactants = {.count = 1, .entries = {{{"molecule", 1}}}};
        input.reactions[0].products = {.count = 1, .entries = {{{"atom", 2}}}};
        input.reactions[0].forward.A = 1;
        const auto mechanism =
          canary::make_reaction_mechanism(input, canary::reaction_species(definitions), si_scales);
        expect(mechanism.mass()[canary::SpeciesIndex {0}].value() ==
               2 * mechanism.mass()[canary::SpeciesIndex {1}].value());
        auto broken = input;
        broken.reactions[0].products.entries[0].value = 1;
        expect(error_message([&] {
                   (void)canary::make_reaction_mechanism(broken,
                                                         canary::reaction_species(definitions),
                                                         si_scales);
               }).find("element X is not conserved") != std::string::npos);
    };
}
