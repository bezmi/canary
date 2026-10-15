#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air11.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/data/air11_kim_jo_exchange.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <numbers>
#include <random>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;
using Gas = canary::Air11ThreeTemperature<double, si>;
using Mechanism = canary::Air11KimJo<double, si>;
using Exchange = canary::Air11KimJoExchange<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;
using Rho = canary::Density<double, si>;
using Temp = canary::Temperature<double, si>;

// A composed model provides the species query without exposing an electron field.
struct SpeciesQueryModel {
    using parameter_type = double;
    using units = si;
    static constexpr std::size_t mode_count = Gas::mode_count;
    Gas thermo;
};

canary::OptionalModeIndex species_vibration_mode(const SpeciesQueryModel& gas,
                                                 canary::SpeciesIndex index)
{ return canary::species_vibration_mode(gas.thermo, index); }

canary::OptionalModeIndex electron_energy_mode(const SpeciesQueryModel& gas)
{ return canary::electron_energy_mode(gas.thermo); }

canary::OptionalSpeciesIndex electron_species_index(const SpeciesQueryModel& gas)
{ return canary::electron_species_index(gas.thermo); }

// Raw integer query results do not satisfy the typed mode contract.
struct RawModeQueryModel {
    using parameter_type = double;
    using units = si;
    static constexpr std::size_t mode_count = Gas::mode_count;
    [[maybe_unused]] friend int species_vibration_mode(const RawModeQueryModel&,
                                                       canary::SpeciesIndex)
    { return -1; }
    [[maybe_unused]] friend int electron_energy_mode(const RawModeQueryModel&) { return -1; }
};

template <typename Index>
concept AcceptsModeEnergyIndex = requires(const Gas& gas, Index index, Temp temp) {
    canary::species_mode_energy(gas, index, temp);
};

static_assert(canary::HasEnergyModes<SpeciesQueryModel>);
static_assert(!canary::HasEnergyModes<RawModeQueryModel>);
static_assert(AcceptsModeEnergyIndex<canary::ModeIndex>);
static_assert(!AcceptsModeEnergyIndex<canary::SpeciesIndex>);
static_assert(!AcceptsModeEnergyIndex<std::size_t>);
static_assert(!AcceptsModeEnergyIndex<canary::OptionalModeIndex>);

constexpr std::size_t species_count = 11;
constexpr std::size_t reaction_count = 40;
constexpr std::size_t mode_count = 2;
constexpr std::size_t electron = 10;
using Fractions = std::array<double, species_count>;
using PerReaction = std::array<double, reaction_count>;
using PerMode = std::array<double, mode_count>;

constexpr double atmosphere = 101325.0;
constexpr double avogadro = canary::constants::avogadro<double>;
constexpr double boltzmann = canary::constants::boltzmann<double>;
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

// A state of a reference implementation of the eleven-species
// three-temperature air model, the mechanism air11_kim_jo and the exchange
// air11_kim_jo_exchange, computed from (rho, T, T_v, T_e, y) with
// R_u = 8.31451 J/(mol K), k_B = 1.380658e-23 J/K, N_A = 6.02214e23 /mol and
// the standard pressure 1 atm: rates of progress in mol/(m^3 s), mass
// production rates in kg/(m^3 s), and the mode sources in W/m^3 (vibration,
// electron) in total and per exchange (vibration-translation,
// chemistry-vibration, electron-translation, vibration-electron). It has no
// electron-impact ionisation source.
struct ReferenceState {
    const char* name;
    Fractions y;
    double rho;
    double temp;
    double t_v;
    double t_e;
    PerReaction forward;
    PerReaction backward;
    Fractions source;
    PerMode total;
    PerMode vt;
    PerMode cv;
    PerMode et;
    PerMode ve;
};

// clang-format off
const std::array<ReferenceState, 4> reference_states {
 {
   {
     .name = "post-shock, T_e between",
     .y = {0.4999995398459613, 0.049999953984596134, 0.19999981593838453, 0.19999981593838453, 0.01999998159383845, 0.009999990796919225, 0.009999990796919225, 0.003999996318767691, 0.0019999981593838453, 0.003999996318767691, 9.203080774345904e-07},
     .rho = 0.001, .temp = 12000.0, .t_v = 5000.0, .t_e = 8000.0,
     .forward = {0.3311303189768377, 0.052627452969181254, 0.022448973269704985, 0.7025898268797712, 1.8043698160862143, 0.002649094430315865, 0.002105134210015285, 0.00448987674212909, 0.035130867313617975, 0.09022158439289515, 3.757021967031462, 1.1129291152244, 0.14030075107560014, 150.97993731443114, 39.65268737543988, 0.030056764352446936, 0.04451792784182268, 0.02806066324679883, 7.549292548637031, 1.9827023532283474, 2.281594223857229, 96.86190041881451, 13.839049306984272, 485.04950711670284, 132.0638911465566, 13.863990185844353, 42.5310028348855, 0.29336840211948706, 505.0594083052646, 26.676769158815784, 0.3577808907733399, 0.3716723115506606, 0.2994478545082911, 0.0008755228249258909, 15.54304492867512, 3.185563637673142, 44.658250254821446, 7.113714033138841, 1.4569416353480358e-05, 1.0065880537400348e-05},
     .backward = {0.00024092765228365137, 3.233845620010255e-05, 1.379441903154147e-05, 0.0005085202957490465, 0.0011087470697178016, 1.9274589646330224e-06, 1.2935604253122927e-06, 2.7589342477628097e-06, 2.5427010686419434e-05, 5.5439254430611795e-05, 7.034232353912862e-05, 7.572516736183032e-06, 2.6268360716401334e-06, 0.0016517966352282494, 0.0004338203919644147, 5.6274960891158e-07, 3.0290586258444556e-07, 5.253768197656986e-07, 8.259306668158507e-05, 2.1691763382449616e-05, 0.0029160286972155993, 0.059519588112821036, 408.70612822019496, 82.53188397798903, 4806.771921048249, 8856.440429922985, 4900.148994199068, 31.059970408170955, 75.88039514271566, 160.92347260393558, 6.311094298776658, 0.3368476267356278, 3.024299066658695, 0.030883887029004877, 9.23703834804622, 55.90978505118614, 4.0074808581059544, 22.51804272893099, 0.5108842557802363, 0.6613783687259123},
     .source = {-27.2363192470121, -24.962422780482374, 211.87048102852359, 366.0184737646054, 24.225885243959357, -8.41499380872576, -0.6073673595083087, -123.34500269987552, -280.0505497386352, -137.48810328826903, -0.010080745275893907},
     .total = {2997132240.282368, -1983546836.5308173},
     .vt = {915035979.6116308, 0.0},
     .cv = {-12993446.727560824, 0.0},
     .et = {0.0, 111542870.86748034},
     .ve = {2095089707.3982973, -2095089707.3982973},
   },
   {
     .name = "strongly ionised, T_e above T",
     .y = {0.10100951469793618, 0.010100951469793618, 0.4040380587917447, 0.3030285440938085, 0.010100951469793618, 0.08080761175834894, 0.0606057088187617, 0.010100951469793618, 0.005050475734896809, 0.015151427204690425, 5.8044904318440994e-06},
     .rho = 0.0001, .temp = 10000.0, .t_v = 9000.0, .t_e = 14000.0,
     .forward = {0.0015436797865539034, 0.00022687414659030279, 0.00024194082009703965, 0.032674095910462116, 0.058339066266077856, 0.00015437100177759523, 0.00011343901813825023, 0.0003629178653544314, 0.006535075141185032, 0.01166821334208387, 0.0059401130670267204, 0.0011010873358985613, 0.0005545631168057373, 1.8612427210022535, 0.36662127295767394, 0.0005940229397420298, 0.000550553106848065, 0.0008318598840605525, 0.37226312461911926, 0.07332676887740838, 0.07597563169821263, 13.036584635310339, 0.07657889474983046, 0.9394322571411249, 1.9794729291111564, 0.05076968986777312, 0.4286845401078314, 0.004245006643104379, 6.1417289925808, 0.8236250103810452, 0.004971439702997369, 0.0009621862755506518, 0.008628112706562255, 6.590530295421265e-05, 0.06410497790492585, 0.0771423898198302, 0.3501637348847107, 0.3012488975273424, 0.1983917438380295, 0.43808762302903503},
     .backward = {2.640476979964928e-07, 3.8025854426307823e-08, 4.055114495441056e-08, 5.585398998099e-06, 9.778076852479154e-06, 2.6405286907190576e-08, 1.9013253183837084e-08, 6.082782954372929e-08, 1.117123553964628e-06, 1.955682428471771e-06, 3.066548742884365e-08, 5.032364819835478e-09, 2.8628997621789197e-09, 9.006930278200952e-06, 1.7741545510286254e-06, 3.0666087977718997e-09, 2.5162255490772146e-09, 4.294428158079952e-09, 1.8014566132374844e-06, 3.548430773440072e-07, 9.725645809047831e-06, 0.0021850319968530997, 1.7194071326865983, 1.0895604578132572, 1260.6689412925534, 1340.1850791009358, 894.3528681530463, 0.4852479774358675, 4.5833221147423515, 2.2714652875516053, 0.07300373293963283, 0.0019802748306447806, 0.03638305194972803, 0.0004423407642706451, 0.1183398204868185, 3.1974360726776045, 0.17090475991796752, 0.9868361769248392, 0.32753568988842613, 0.5269828380156827},
     .source = {-0.41272524296158186, -0.20150460538791815, 43.0254188221097, 63.18874878052897, 0.042924673698658634, -0.05104090367064871, -0.004295723638359731, -25.011922117240974, -42.93628486038126, -37.63740258190869, -0.0019161709497631859},
     .total = {78088448.06756195, -102083005.06012443},
     .vt = {395085.6172234663, 0.0},
     .cv = {-629840.3181166119, 0.0},
     .et = {0.0, -23759802.291669354},
     .ve = {78323202.76845509, -78323202.76845509},
   },
   {
     .name = "equal temperatures",
     .y = {0.4999995398459613, 0.049999953984596134, 0.19999981593838453, 0.19999981593838453, 0.01999998159383845, 0.009999990796919225, 0.009999990796919225, 0.003999996318767691, 0.0019999981593838453, 0.003999996318767691, 9.203080774345904e-07},
     .rho = 0.01, .temp = 7000.0, .t_v = 7000.0, .t_e = 7000.0,
     .forward = {7.889144766709588, 1.303823119453573, 0.5561639164672647, 16.75948102778312, 44.70250695269393, 0.06311439413347045, 0.05215381892370655, 0.11123481699395397, 0.838006873579311, 2.235201989997017, 170.3436430983337, 63.773972308542156, 6.361246028735305, 7751.3888796605115, 2035.789691252314, 1.362775832686811, 2.5510026277330624, 1.2722724666398533, 387.5846244984217, 101.79297492013166, 81.0823834347926, 2399.713040154466, 434.8174756767775, 8456.607464613888, 1153.608333910327, 2.6686038443610265, 34.091619120509314, 0.2263003286065131, 18660.374091866433, 450.84465803486006, 2.6377043760877337, 3.296325761356547, 10.766466117761034, 0.07314527909117138, 178.9879879766242, 15.09907219039365, 946.7429829013997, 85.978685529195, 3.503848696933863, 4.2536308623646715},
     .backward = {0.5312153652619788, 0.08779289709072192, 0.037449283384732995, 1.1284992352204282, 3.0100421859676088, 0.004249806148110983, 0.003511776090898008, 0.007489993616118897, 0.05642717184237337, 0.15050727000992756, 0.05765297463681445, 0.021584364060306126, 0.002152970015684728, 2.6234652397399887, 0.6890150363163513, 0.0004612328296407379, 0.000863389364698752, 0.000430601875809138, 0.13117839984746307, 0.03445193313077409, 2.7488219663526094, 161.5846174516668, 43554.32705916809, 16796.42166225658, 637445.63524915, 673755.8295727252, 735642.24936929, 2911.9603050855117, 11266.76570650307, 12696.778185758749, 252.580141278417, 56.05000376509645, 317.1692081479662, 4.948720076173864, 816.4568369791594, 6898.964705850211, 546.4617599826125, 3357.569956613575, 8820.179049524282, 10236.24377282351},
     .source = {-48.009907730825034, -2200.4764214543184, 29232.077159515953, 33157.48764169304, 1029.6243767644144, -462.66811943265355, -148.4420198121975, -20481.27732131959, -21396.25174079867, -18680.930951901937, -1.1326540289057105},
     .total = {-2577242091.622133, 0.0},
     .vt = {0.0, 0.0},
     .cv = {-2577242091.6221323, 0.0},
     .et = {0.0, 0.0},
     .ve = {0.0, 0.0},
   },
   {
     .name = "expansion, frozen modes",
     .y = {0.10100951469793618, 0.010100951469793618, 0.4040380587917447, 0.3030285440938085, 0.010100951469793618, 0.08080761175834894, 0.0606057088187617, 0.010100951469793618, 0.005050475734896809, 0.015151427204690425, 5.8044904318440994e-06},
     .rho = 0.0001, .temp = 3000.0, .t_v = 6000.0, .t_e = 5000.0,
     .forward = {1.6119429237806475e-09, 3.232054402411065e-10, 3.4466945858299213e-10, 3.445004886386648e-08, 8.31099703477131e-08, 1.6119744918590038e-10, 1.6160549074858384e-10, 5.170136404085771e-10, 6.8902796441500495e-09, 1.6622564037080822e-08, 3.178521125121944e-06, 3.7879241777514213e-06, 2.967436077547268e-07, 0.002674184052296462, 0.0005267516967631719, 3.1785833728618953e-07, 1.8939945602626306e-06, 4.451235498033094e-07, 0.0005348577592171515, 0.00010535395181172432, 4.221261225721688e-06, 1.857194898414229e-05, 0.0008187210757294595, 0.00040222184689005854, 0.00034759396528702354, 1.3372122359842202e-11, 1.0178449205757974e-08, 3.494199797631361e-11, 0.19523588761241895, 0.0008797001360215288, 1.017399397211982e-06, 7.226187532637879e-08, 0.00014419097305196522, 1.2204965412864353e-05, 1.945262032123211e-05, 6.467833448559583e-07, 0.001110678010913019, 7.612363316936983e-05, 0.00013770330574934254, 0.00017633476605538178},
     .backward = {1.1429905806240397e-06, 2.6198139040863747e-07, 2.7937953000916177e-07, 2.4529517944916796e-05, 6.736664324793535e-05, 1.1430129648013435e-07, 1.309929409987683e-07, 4.190769569764301e-07, 4.906095745888711e-06, 1.347379064950896e-05, 2.2944264099129995e-08, 6.093982306806669e-08, 2.142053942080338e-09, 2.9538091646649228e-05, 5.818313021744715e-06, 2.2944713436582476e-09, 3.0470433931124946e-08, 3.2131396587894473e-09, 5.907849721903216e-06, 1.1637025062949504e-06, 8.154301180355324e-06, 0.015053908170092325, 1.1151976264789052, 3.9601633233724916, 2647.7476924155444, 491.8640203934629, 1844.6109958379438, 0.4101435460550079, 5.3070553884692595, 2.138368083983717, 0.010457079400294212, 0.004774754006389753, 0.04644968934765906, 0.0014156579418250774, 0.0689608470267417, 3.1231955058862617, 0.22781651622650215, 1.8110828453228263, 165.26129034779856, 304.96261612254784},
     .source = {0.26099539880613537, -0.11706486432372398, 92.90002912533828, 60.85148154715886, -0.0871336214325723, -4.235582152831216, -2.640957593775767, -51.87322624857706, -15.77090099039122, -79.28464829092164, -0.0029921994320760955},
     .total = {-4772335.081999166, -37480681.79179475},
     .vt = {-50428.04397051974, 0.0},
     .cv = {-1626.6585075406067, 0.0},
     .et = {0.0, -42200962.17131586},
     .ve = {-4720280.3795211045, 4720280.3795211045},
   },
 },
};
// clang-format on

// Standard enthalpies over R_u, in K, by which the equilibrium constants of
// the reference differ from those of canary, for the molecules N2, O2 and
// NO and their ions, relative to the atoms, the atomic ions and e- (0): the
// reference's ln K_c is that of canary minus sum_i dnu_i h_i / T at every
// state. They are found from reactions 0, 10 and 22 of the first reference
// state; with them, the equilibrium backward rates of every other reaction
// at T_l in the four states agree to 1e-10.
constexpr Fractions enthalpy_offset {
  0.03869508856,
  1.109733591,
  0.0,
  0.0,
  0.2799083787,
  0.0,
  0.0,
  0.03869508856,
  1.109733591,
  0.2799083787,
  0.0,
};

// Change of the coefficient of species i in reaction k.
template <typename M> int coefficient_change(const M& mech, std::size_t k, std::size_t i)
{ return mech.reactions[k].products[i] - mech.reactions[k].reactants[i]; }

// The state of gas with density rho, translational temperature temp and
// mode temperatures (t_v, t_e), and the composition of y.
template <typename T, typename G>
auto cell_state(const G& gas, const std::array<T, species_count>& y, T rho, T temp, T t_v, T t_e)
{
    using U = G::units;
    const std::array<T, mode_count> t_m {t_v, t_e};
    return canary::thermo_state(gas,
                                canary::mixture(gas, fractions(y)),
                                canary::Density<T, U> {rho},
                                canary::Temperature<T, U> {temp},
                                canary::as_kind<canary::kind::mode_temperature, U>(
                                  Kokkos::mdspan<const T, Extents>(t_m.data(), mode_count)))
      .state;
}

// Number densities n_i = rho y_i / m_i.
template <typename S> Fractions number_densities(const S& species, const Fractions& y, double rho)
{
    Fractions n {};
    for (std::size_t i = 0; i < species_count; ++i) { n[i] = rho * y[i] / species.mass[i].value(); }
    return n;
}

// [M] of reaction k (1 without a third body) times prod_i n_i^nu_i over
// the reactants (forward) or the products (backward).
double concentration_product(const Mechanism& mech, std::size_t k, const Fractions& n, bool forward)
{
    const auto& r = mech.reactions[k];
    double out = 1.0;
    if (r.third_body) {
        out = 0.0;
        for (std::size_t i = 0; i < species_count; ++i) { out += r.efficiencies[i] * n[i]; }
    }
    for (std::size_t i = 0; i < species_count; ++i) {
        const int nu = forward ? r.reactants[i] : r.products[i];
        for (int c = 0; c < nu; ++c) { out *= n[i]; }
    }
    return out;
}

// ln(k_f) of reaction k at the rate-controlling temperature T_c, from the
// formula A T_c^n exp(-theta / T_c).
double forward_log_rate(const Mechanism& mech, std::size_t k, double t_c)
{
    const auto& f = mech.reactions[k].forward;
    return std::log(f.A) + (f.n * std::log(t_c)) - (f.activation_temperature.value() / t_c);
}

// k_b of an equilibrium backward rate at the temperature t_eq from the
// formula k_f(t_eq) / K_c(t_eq), with the Gibbs energies of gas at t_eq.
double equilibrium_backward_rate(const Mechanism& mech, const Gas& gas, std::size_t k, double t_eq)
{
    const auto props = canary::species_properties(gas, Temp {t_eq});
    double ln_k_c =
      mech.reactions[k].coefficient_change * std::log(atmosphere / (boltzmann * t_eq));
    for (std::size_t i = 0; i < species_count; ++i) {
        ln_k_c -= coefficient_change(mech, k, i) * props.gibbs_energy[i].value() *
                  mech.mass[i].value() / (boltzmann * t_eq);
    }
    return std::exp(forward_log_rate(mech, k, t_eq) - ln_k_c);
}

// Rates of progress and the mode sources of canary in one cell, as plain
// values: the sources of each exchange and the total of
// multi_temperature_sources.
struct Values {
    PerReaction forward;
    PerReaction backward;
    Fractions production;
    PerMode total;
    PerMode vt;
    PerMode cv;
    PerMode et;
    PerMode ve;
    PerMode ce;
};

template <typename M, typename X, typename G>
Values evaluate(const M& mech, const X& exchange, const G& gas, const ReferenceState& s)
{
    using U = M::units;
    const auto mix = canary::mixture(gas, fractions(s.y));
    const auto state = cell_state(gas, s.y, s.rho, s.temp, s.t_v, s.t_e);
    const auto q =
      canary::rates_of_progress(mech, gas, mix, state.rho, state.temp, state.mode_temp);
    const auto sources = canary::multi_temperature_sources(mech, exchange, gas, mix, state);
    const auto vt = canary::vibration_translation_source(exchange, gas, mix, state);
    const auto cv = canary::chemistry_vibration_source(exchange, mech, q, state);
    const auto et = canary::electron_translation_source(exchange, mix, state);
    const auto ve = canary::vibration_electron_source(exchange, gas, mix, state);
    const auto ce = canary::electron_impact_ionisation_source(exchange, q);
    static_assert(
      std::same_as<typename decltype(sources.modes)::value_type, canary::PowerDensity<double, U>>);
    Values out {};
    for (std::size_t k = 0; k < reaction_count; ++k) {
        out.forward[k] = q.forward[k].value();
        out.backward[k] = q.backward[k].value();
    }
    for (std::size_t i = 0; i < species_count; ++i) {
        out.production[i] = sources.species[i].value();
    }
    for (const canary::ModeIndex m : canary::mode_indices(0, mode_count)) {
        out.total[m.value()] = sources.modes[m.value()].value();
        out.vt[m.value()] = vt[m.value()].value();
        out.cv[m.value()] = cv[m.value()].value();
        out.et[m.value()] = et[m.value()].value();
        out.ve[m.value()] = ve[m.value()].value();
        out.ce[m.value()] = ce[m.value()].value();
    }
    return out;
}

// The reactions whose forward rate depends on T_e, and those whose
// equilibrium backward rate is evaluated at T_e.
constexpr std::array<std::size_t, 3> forward_at_electron {21, 38, 39};
constexpr std::array<std::size_t, 6> backward_at_electron {21, 24, 25, 26, 38, 39};

bool contains(const auto& list, std::size_t k) { return std::ranges::find(list, k) != list.end(); }

template <typename F> std::string error_message(const F& f)
{
    try {
        f();
    } catch (const canary::ParameterError& error) {
        return error.what();
    }
    return "no exception";
}

// Up to three species coefficients.
canary::SpeciesCoefficients<3> coefficients(std::initializer_list<canary::SpeciesCoefficient> list)
{
    canary::SpeciesCoefficients<3> out {};
    out.count = static_cast<int>(list.size());
    std::ranges::copy(list, out.entries.begin());
    return out;
}

// A reaction without a third body, in SI per mole, for the mode resolution
// tests.
canary::ReactionInput<3> reaction(const char* equation,
                                  canary::SpeciesCoefficients<3> reactants,
                                  canary::SpeciesCoefficients<3> products,
                                  canary::ArrheniusInput forward,
                                  canary::BackwardRateModel backward_model,
                                  canary::ArrheniusInput backward)
{
    canary::ReactionInput<3> r {};
    r.equation = equation;
    r.reactants = reactants;
    r.products = products;
    r.forward = forward;
    r.backward_model = backward_model;
    r.backward = backward;
    r.third_body = false;
    return r;
}

// The eleven species of air with three modes: the vibration of N2; the
// vibration of O2 and NO; the electronic energy of the heavy species and
// e-. The vibrations of the molecular ions are frozen.
const std::array<canary::EnergyModeInput<11>, 3> split_modes {
  {
    {
      .count = 1,
      .components = {{{.species = "N2", .component = canary::ModeComponent::vibration}}},
    },
    {
      .count = 2,
      .components =
        {
          {
            {.species = "O2", .component = canary::ModeComponent::vibration},
            {.species = "NO", .component = canary::ModeComponent::vibration},
          },
        },
    },
    canary::air11_three_temperature_modes[1],
  },
};

// One mode holding every component of the eleven species.
constexpr auto single_mode = [] {
    std::array<canary::EnergyModeInput<17>, 1> modes {};
    std::size_t count = 0;
    for (const char* name : {"N2", "O2", "NO", "N2+", "O2+", "NO+"}) {
        modes[0].components[count++] = {
          .species = name,
          .component = canary::ModeComponent::vibration,
        };
    }
    for (std::size_t c = 0; c < 11; ++c) {
        modes[0].components[count++] = canary::air11_three_temperature_modes[1].components[c];
    }
    modes[0].count = static_cast<int>(count);
    return modes;
}();

}  // namespace

static_assert(canary::HasEnergyModes<Gas>);
static_assert(canary::HasEnergyModes<canary::Air5TwoTemperature<double, si>>);
static_assert(!canary::HasEnergyModes<canary::Air11<double, si>>);
static_assert(canary::HasModeEnergies<Gas, double, species_count>);
static_assert(canary::HasModeEnergies<Gas, D, species_count>);
static_assert(canary::HasModeEnergies<canary::Air5TwoTemperature<double, si>, double, 5>);
static_assert(!canary::HasModeEnergies<Gas, double, 5>);
static_assert(canary::ModeRelaxingModel<Gas, D, double, species_count>);
static_assert(!canary::ModeRelaxingModel<canary::Air11<double, si>, double, double, species_count>);
static_assert(Gas::mode_count == 2 && canary::Air5TwoTemperature<double, si>::mode_count == 1);
static_assert(std::is_trivially_copyable_v<Exchange>);
static_assert(std::is_trivially_copyable_v<canary::Air11KimJoExchange<float, scaled>>);
static_assert(std::is_trivially_copyable_v<canary::MultiTemperatureSources<D, si, 11, 2>>);
static_assert(Exchange::mode_count == 2);

namespace
{
template <typename V>
concept AcceptsReactionModes = requires(
  const Gas& gas, const Mechanism& mechanism, const Gas::Mixture<double>& mix, const V& modes) {
    canary::rates_of_progress(mechanism, gas, mix, Rho {1}, Temp {4000}, modes);
    canary::production_rates(mechanism, gas, mix, Rho {1}, Temp {4000}, modes);
};
using RawModes = canary::ModeTemperatureView<double, si>;
static_assert(AcceptsReactionModes<canary::IndexedView<canary::ModeIndex, RawModes>>);
static_assert(!AcceptsReactionModes<canary::IndexedView<canary::SpeciesIndex, RawModes>>);
static_assert(
  AcceptsReactionModes<std::array<canary::ModeTemperature<double, si>, Gas::mode_count>>);
static_assert(
  !AcceptsReactionModes<std::array<canary::ModeTemperature<double, si>, Gas::mode_count + 1>>);
static_assert(!AcceptsReactionModes<canary::ModeTemperatureView<double, scaled>>);
static_assert(!AcceptsReactionModes<canary::ModeEnergyDensityView<double, si>>);
}  // namespace

void test_multi_temperature_kinetics()
{
    using namespace boost::ut;
    const canary::UnitScales<double, si> si_scales {};
    const auto built = canary::make_air11_three_temperature<double>(si_scales);
    const Gas& gas = built.model;
    const canary::Pressure<double, si> p0 {atmosphere};
    const Mechanism mech = canary::make_air11_kim_jo(si_scales, p0, gas);
    const Exchange exchange =
      canary::make_air11_kim_jo_exchange(built.species, mech, gas, si_scales);
    const std::array<std::string_view, species_count> names =
      canary::detail::air11_kim_jo_species_names;

    "mode views preserve closure and reaction results across cell layouts"_test = [&] {
        const auto run = [&]<typename T, typename Layout>() {
            const auto equal = [](const T& a, const T& b) {
                expect(canary::ScalarTraits<T>::value(a) == canary::ScalarTraits<T>::value(b));
                if constexpr (requires { a.dual(); }) { expect(a.dual() == b.dual()); }
            };
            std::array<T, species_count> y {};
            for (const auto i : canary::species_indices(0, species_count)) {
                y[i.value()] = T(reference_states[0].y[i.value()]);
            }
            const auto mix = canary::mixture(gas, fractions(y));
            const auto seed = [](double v) {
                if constexpr (requires { T(v, 1.0); }) {
                    return T(v, 1.0);
                } else {
                    return T(v);
                }
            };
            const std::array<canary::ModeTemperature<T, si>, 2> modes {
              canary::ModeTemperature<T, si> {seed(5000)},
              canary::ModeTemperature<T, si> {T(8000)}};
            std::array<T, 8> storage {};
            using Grid = Kokkos::mdspan<T, Kokkos::dextents<std::size_t, 2>, Layout>;
            const Grid grid(storage.data() + 1, 3, 2);
            for (std::size_t cell = 0; cell < 3; ++cell) {
                for (const auto m : canary::mode_indices(0, 2)) {
                    storage[1 + grid.mapping()(cell, m.value())] =
                      cell == 1 ? modes[m.value()].value() : T(99);
                }
            }
            using Mapping = Kokkos::layout_stride::mapping<Extents>;
            const Kokkos::mdspan<const T, Extents, Kokkos::layout_stride> raw(
              storage.data() + 1 + grid.mapping()(1, 0),
              Mapping(Extents(2), std::array<std::size_t, 1> {grid.stride(1)}));
            const auto quantities = canary::as_kind<canary::kind::mode_temperature, si>(raw);
            const auto view = canary::as_indexed<canary::ModeIndex>(quantities);
            const canary::Density<T, si> rho {T(0.001)};
            const canary::Temperature<T, si> temp {T(10000)};
            const auto expected = canary::rates_of_progress(mech, gas, mix, rho, temp, modes);
            const auto typed = canary::rates_of_progress(mech, gas, mix, rho, temp, view);
            const auto untyped = canary::rates_of_progress(mech, gas, mix, rho, temp, quantities);
            for (std::size_t k = 0; k < reaction_count; ++k) {
                equal(typed.forward[k].value(), expected.forward[k].value());
                equal(typed.backward[k].value(), expected.backward[k].value());
                equal(untyped.forward[k].value(), expected.forward[k].value());
                equal(untyped.backward[k].value(), expected.backward[k].value());
            }
            const auto omega = canary::production_rates(mech, gas, mix, rho, temp, view);
            const auto omega_expected = canary::production_rates(mech, gas, mix, rho, temp, modes);
            for (const auto i : canary::species_indices(0, species_count)) {
                equal(omega[i.value()].value(), omega_expected[i.value()].value());
            }
            const auto c = canary::thermo_state(gas, mix, view, temp, rho);
            const auto base = canary::thermo_state(gas, mix, rho, temp, quantities);
            expect(c.status == base.status && c.diagnostic.reason == base.diagnostic.reason);
            equal(c.state.e.value(), base.state.e.value());
            equal(c.state.p.value(), base.state.p.value());
            for (const auto m : canary::mode_indices(0, 2)) {
                equal(c.state.mode_e[m.value()].value(), base.state.mode_e[m.value()].value());
            }
        };
        run.template operator()<double, Kokkos::layout_left>();
        run.template operator()<double, Kokkos::layout_right>();
        run.template operator()<D, Kokkos::layout_left>();
        run.template operator()<D, Kokkos::layout_right>();
    };

    "energy exchange resolves electron identity through the typed model query"_test = [&] {
        auto renamed = names;
        renamed[10] = "carrier";
        renamed[2] = "e-";
        const canary::EnergyExchangeInput<0, 0> input {};
        canary::ElectronExchangeInput<1, 0, 0> electrons {};
        electrons.electron_translation[0] =
          canary::air11_kim_jo_exchange_electrons.electron_translation[0];
        electrons.electron_translation[0].partner = "e-";
        const SpeciesQueryModel custom {gas};
        const auto result = canary::make_multi_temperature_exchange(input,
                                                                    electrons,
                                                                    renamed,
                                                                    built.species,
                                                                    mech,
                                                                    custom,
                                                                    si_scales);
        expect(result.electron == canary::OptionalSpeciesIndex {canary::SpeciesIndex {10}});
        expect(result.electron_translation[0].partner == canary::SpeciesIndex {2});
        electrons.electron_translation[0].partner = "carrier";
        expect(error_message([&] {
                   canary::make_multi_temperature_exchange(input,
                                                           electrons,
                                                           renamed,
                                                           built.species,
                                                           mech,
                                                           custom,
                                                           si_scales);
               })
                 .starts_with(
                   "electron-translation pair 0 (carrier): the partner must be a heavy species"));
        auto invalid = custom;
        for (const std::size_t index :
             {std::size_t(11), std::size_t(std::numeric_limits<std::uint32_t>::max()) - 1}) {
            invalid.thermo.electron = canary::OptionalSpeciesIndex {canary::SpeciesIndex {index}};
            expect(error_message([&] {
                       canary::make_multi_temperature_exchange(input,
                                                               electrons,
                                                               renamed,
                                                               built.species,
                                                               mech,
                                                               invalid,
                                                               si_scales);
                   }) == "electron species index is outside the species range");
        }
        const auto no_electrons = canary::make_air5_two_temperature(si_scales).model;
        expect(!canary::electron_species_index(no_electrons));
        expect(canary::species_vibration_mode(gas, canary::SpeciesIndex {0}) ==
               canary::OptionalModeIndex {canary::ModeIndex {0}});
        expect(!canary::species_vibration_mode(gas, canary::SpeciesIndex {11}));
    };

    "multi-temperature kinetics: the builders resolve modes and electrons"_test = [&] {
        for (std::size_t k = 0; k < reaction_count; ++k) {
            const auto& r = mech.reactions[k];
            const bool forward_electron = r.forward.translational_exponent != 1.0 &&
                                          r.forward.other_mode == canary::ModeIndex {1};
            expect(forward_electron == contains(forward_at_electron, k)) << k;
            const bool backward_electron = r.equilibrium_at_other_temperature;
            expect(backward_electron == contains(backward_at_electron, k)) << k;
            if (backward_electron) { expect(r.backward.other_mode == canary::ModeIndex {1}) << k; }
            // Dissociations of N2, O2 and NO: their vibration mode.
            if (k <= 20) { expect(r.forward.other_mode == canary::ModeIndex {0}) << k; }
        }
        expect(exchange.electron == canary::OptionalSpeciesIndex {canary::SpeciesIndex {10}} &&
               exchange.electron_mode == canary::OptionalModeIndex {canary::ModeIndex {1}});
        for (const std::size_t i : std::array<std::size_t, 6> {0, 1, 4, 7, 8, 9}) {
            expect(exchange.vibration_mode[i] == canary::OptionalModeIndex {canary::ModeIndex {0}})
              << i;
        }
        for (const std::size_t i : std::array<std::size_t, 5> {2, 3, 5, 6, 10}) {
            expect(!exchange.vibration_mode[i]) << i;
        }
        const auto& ve = exchange.vibration_electron[0];
        expect(ve.species == canary::SpeciesIndex {0});
        // In SI the fit is that of the input with log10(101325) added to c,
        // for p_e in Pa.
        expect(ve.relaxation_time.a_low == 7.5_d && ve.relaxation_time.b_high == -17.9_d);
        expect(relative(ve.relaxation_time.c_high, 24.35 + std::log10(atmosphere), 1.0e-15));
        expect(ve.relaxation_time.blend_lower.value() == 6950.0_d);
        const auto& ce = exchange.electron_impact_ionisation[1];
        expect(ce.reaction == 39_ul && ce.species == canary::SpeciesIndex {2});
        expect(relative(ce.ionisation_energy, 1402330.360173364 / avogadro, 1.0e-15));
        // (2 pi / 9) C^2 and 9 / (4 pi C^3) with C = e^2 / (4 pi eps_0 k_B),
        // the exact e and k_B and CODATA 2018 eps_0; a few ulp.
        const double e = canary::constants::elementary_charge<double>;
        const double c = e * e / (4.0 * std::numbers::pi * 8.8541878128e-12 * boltzmann);
        expect(
          relative(exchange.coulomb_cross_section, 2.0 * std::numbers::pi / 9.0 * c * c, 1.0e-15));
        expect(relative(exchange.coulomb_logarithm_factor,
                        9.0 / (4.0 * std::numbers::pi * c * c * c),
                        1.0e-15));
        expect(exchange.minimum_electron_density == 1.0_d);
        expect(exchange.electron_translation[0].cross_section ==
               canary::ElectronCrossSection::gnoffo_neutral);
        expect(exchange.electron_translation[5].cross_section ==
               canary::ElectronCrossSection::coulomb);
        // kim_2021 of N2-N2 with omega = 0.55.
        expect(relative(exchange.vibration_translation[0].correction_exponent, -1.05, 1.0e-15));
    };

    "multi-temperature kinetics: rates and mode sources equal the reference"_test = [&] {
        const double r_u =
          canary::constants::universal_gas_constant<double> / reference_gas_constant;
        for (const ReferenceState& s : reference_states) {
            const Values got = evaluate(mech, exchange, gas, s);
            const Fractions n = number_densities(built.species, s.y, s.rho);
            const double t_l = std::clamp(s.temp, 300.0, 30000.0);
            // The reference rates in particles per volume and time, with
            // the reactions that differ in the reference replaced by the
            // formulas at canary's temperatures.
            PerReaction forward {};
            PerReaction backward {};
            for (std::size_t k = 0; k < reaction_count; ++k) {
                const int dnu = mech.reactions[k].coefficient_change;
                forward[k] = avogadro * s.forward[k];
                double offset = 0.0;
                for (std::size_t i = 0; i < species_count; ++i) {
                    offset += coefficient_change(mech, k, i) * enthalpy_offset[i];
                }
                backward[k] =
                  avogadro * s.backward[k] * std::pow(r_u, dnu) * std::exp(offset / t_l);
                // The reference evaluates these forward rates at T_v; the
                // formula at T_c = T^s T_e^(1 - s).
                if (contains(forward_at_electron, k)) {
                    const double s_t = mech.reactions[k].forward.translational_exponent;
                    const double t_c = std::pow(t_l, s_t) * std::pow(s.t_e, 1.0 - s_t);
                    forward[k] = std::exp(forward_log_rate(mech, k, t_c)) *
                                 concentration_product(mech, k, n, true);
                }
                // The reference evaluates these equilibrium backward rates at
                // T; the formula at T_e.
                if (contains(backward_at_electron, k)) {
                    backward[k] = equilibrium_backward_rate(mech, gas, k, s.t_e) *
                                  concentration_product(mech, k, n, false);
                }
                // Forward rates contain no constant that differs between the
                // implementations, and the formulas compute the same
                // expression in another order. Backward rates contain R_u
                // through the Gibbs energies, corrected above, and differ by
                // rounding amplified by the exponent of K_c (up to about
                // 100). The third-body efficiencies of reaction 20 in the
                // reference differ in the eighth digit: 2e-8 in its rates.
                const double bound_f = k == 20 ? 1.0e-7 : 1.0e-12;
                const double bound_b = k == 20 ? 1.0e-7 : 1.0e-10;
                expect(relative(got.forward[k], forward[k], bound_f))
                  << s.name << k << got.forward[k] << forward[k];
                expect(relative(got.backward[k], backward[k], bound_b))
                  << s.name << k << got.backward[k] << backward[k];
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                // The bound is relative to the sum of the magnitudes of the
                // terms, the rounding scale of the backward rates, plus the
                // difference of reaction 20.
                const double m = built.species.mass[i].value();
                double expected = 0;
                double scale = 0;
                for (std::size_t k = 0; k < reaction_count; ++k) {
                    const int dnu = coefficient_change(mech, k, i);
                    expected += m * dnu * (forward[k] - backward[k]);
                    const double terms = m * std::abs(dnu) * (forward[k] + backward[k]);
                    scale += (k == 20 ? 1.0e3 : 1.0) * terms;
                }
                expect(std::abs(got.production[i] - expected) <= 1.0e-10 * scale)
                  << s.name << i << got.production[i] << expected;
            }
            // k_B, R_u and N_A of the reference differ from the exact values
            // by up to 7e-6 relative; they enter the energies, the
            // relaxation times and G linearly or through exponents of order
            // 1. The chemistry-vibration source has the size of the sum of
            // its terms k_B theta_D (q_f + q_b) / 2.
            double cv_scale = 0;
            for (const auto& c : exchange.chemistry_vibration) {
                cv_scale += boltzmann * c.dissociation_temperature.value() *
                            (got.forward[c.reaction] + got.backward[c.reaction]) / 2.0;
            }
            for (const canary::ModeIndex m : canary::mode_indices(0, mode_count)) {
                expect(std::abs(got.vt[m.value()] - s.vt[m.value()]) <=
                       5.0e-5 * std::abs(s.vt[m.value()]))
                  << s.name << m.value() << got.vt[m.value()] << s.vt[m.value()];
                expect(std::abs(got.cv[m.value()] - s.cv[m.value()]) <= 5.0e-5 * cv_scale)
                  << s.name << m.value() << got.cv[m.value()] << s.cv[m.value()];
                expect(std::abs(got.et[m.value()] - s.et[m.value()]) <=
                       5.0e-5 * std::abs(s.et[m.value()]))
                  << s.name << m.value() << got.et[m.value()] << s.et[m.value()];
                expect(std::abs(got.ve[m.value()] - s.ve[m.value()]) <=
                       5.0e-5 * std::abs(s.ve[m.value()]))
                  << s.name << m.value() << got.ve[m.value()] << s.ve[m.value()];
                const double scale = std::abs(s.vt[m.value()]) + cv_scale +
                                     std::abs(s.et[m.value()]) + std::abs(s.ve[m.value()]) +
                                     std::abs(got.ce[m.value()]);
                expect(std::abs(got.total[m.value()] - (s.total[m.value()] + got.ce[m.value()])) <=
                       5.0e-5 * scale)
                  << s.name << m.value() << got.total[m.value()]
                  << s.total[m.value()] + got.ce[m.value()];
            }
            // Q_CE from its definition with canary's rates of progress.
            double ce = 0;
            double ce_scale = 0;
            for (const auto& [k, energy] : {
                   std::pair {std::size_t {38}, 1313942.077121842},
                   std::pair {std::size_t {39}, 1402330.360173364},
                 }) {
                ce -= (got.forward[k] - got.backward[k]) * energy / avogadro;
                ce_scale += (got.forward[k] + got.backward[k]) * energy / avogadro;
            }
            // Two terms and the conversion of I to J per particle.
            expect(got.ce[0] == 0.0_d);
            expect(std::abs(got.ce[1] - ce) <= 1.0e-14 * ce_scale) << s.name << got.ce[1] << ce;
        }
    };

    "multi-temperature kinetics: one mode gives the two-temperature rates"_test = [&] {
        const auto one = canary::make_multi_temperature_mixture(canary::air11_species,
                                                                canary::air11_rotational_degrees,
                                                                canary::air11_electronic_levels,
                                                                single_mode,
                                                                si_scales);
        static_assert(decltype(one.model)::mode_count == 1);
        const Mechanism with_gas = canary::make_air11_kim_jo(si_scales, p0, one.model);
        const Mechanism without = canary::make_air11_kim_jo(si_scales, p0);
        for (std::size_t k = 0; k < reaction_count; ++k) {
            expect(with_gas.reactions[k].forward.other_mode == canary::ModeIndex {0});
            expect(with_gas.reactions[k].backward.other_mode == canary::ModeIndex {0});
            expect(with_gas.reactions[k].equilibrium_at_other_temperature ==
                   without.reactions[k].equilibrium_at_other_temperature);
        }
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(one.model, fractions(s.y));
            const std::array t_m {canary::ModeTemperature<double, si> {s.t_v}};
            const auto modes =
              canary::rates_of_progress(with_gas, one.model, mix, Rho {s.rho}, Temp {s.temp}, t_m);
            const auto two = canary::rates_of_progress(without,
                                                       one.model,
                                                       mix,
                                                       Rho {s.rho},
                                                       Temp {s.temp},
                                                       Temp {s.t_v});
            // The three-temperature mechanism with both modes at T_v.
            const std::array same {
              canary::ModeTemperature<double, si> {s.t_v},
              canary::ModeTemperature<double, si> {s.t_v},
            };
            const auto three =
              canary::rates_of_progress(mech, gas, mix, Rho {s.rho}, Temp {s.temp}, same);
            for (std::size_t k = 0; k < reaction_count; ++k) {
                expect(modes.forward[k] == two.forward[k] && modes.backward[k] == two.backward[k])
                  << s.name << k;
                expect(three.forward[k] == two.forward[k] && three.backward[k] == two.backward[k])
                  << s.name << k;
            }
            const auto omega_modes =
              canary::production_rates(with_gas, one.model, mix, Rho {s.rho}, Temp {s.temp}, t_m);
            const auto omega_two = canary::production_rates(without,
                                                            one.model,
                                                            mix,
                                                            Rho {s.rho},
                                                            Temp {s.temp},
                                                            Temp {s.t_v});
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(omega_modes[i] == omega_two[i]) << s.name << i;
            }
        }
    };

    "multi-temperature kinetics: other temperatures resolve to modes"_test = [&] {
        using canary::BackwardRateModel;
        const auto split = canary::make_multi_temperature_mixture(canary::air11_species,
                                                                  canary::air11_rotational_degrees,
                                                                  canary::air11_electronic_levels,
                                                                  split_modes,
                                                                  si_scales);
        const auto rate = [](double s, const char* other) {
            return canary::ArrheniusInput {
              .A = 1.0e15,
              .n = -1.0,
              .activation_temperature = 50000.0,
              .translational_exponent = s,
              .other_temperature = other,
            };
        };
        const canary::ArrheniusInput none {};
        const auto build = [&](const canary::ReactionInput<3>& r, const auto& model) {
            const canary::MechanismInput<1, 3> input {
              .temperature_lower = 300.0,
              .temperature_upper = 30000.0,
              .reactions = {r},
            };
            return canary::make_reaction_mechanism(input,
                                                   canary::reaction_species(canary::air11_species),
                                                   si_scales,
                                                   p0,
                                                   model)
              .reactions[0];
        };
        const auto o2 = [&](canary::ArrheniusInput f,
                            BackwardRateModel model = BackwardRateModel::equilibrium,
                            canary::ArrheniusInput b = {}) {
            return reaction(
              "O + O2 <=> O + O + O",
              coefficients({{.name = "O", .value = 1.0}, {.name = "O2", .value = 1.0}}),
              coefficients({{.name = "O", .value = 3.0}}),
              f,
              model,
              b);
        };
        const auto ionise = [&](canary::ArrheniusInput f) {
            return reaction(
              "N + O <=> NO+ + e-",
              coefficients({{.name = "N", .value = 1.0}, {.name = "O", .value = 1.0}}),
              coefficients({{.name = "NO+", .value = 1.0}, {.name = "e-", .value = 1.0}}),
              f,
              BackwardRateModel::equilibrium,
              none);
        };
        // Without other_temperature: the unique reactant vibration mode.
        expect(build(o2(rate(0.5, nullptr)), split.model).forward.other_mode ==
               canary::ModeIndex {1});
        const auto dissociation = reaction(
          "N2 + O2 <=> N + N + O2",
          coefficients({{.name = "N2", .value = 1.0}, {.name = "O2", .value = 1.0}}),
          coefficients({{.name = "N", .value = 2.0}, {.name = "O2", .value = 1.0}}),
          rate(0.5, nullptr),
          BackwardRateModel::equilibrium,
          none);
        for (const bool reorder_modes : {false, true}) {
            auto modes = split_modes;
            if (reorder_modes) { std::swap(modes[0], modes[1]); }
            const auto reordered = canary::make_multi_temperature_mixture(
              canary::air11_species, canary::air11_rotational_degrees,
              canary::air11_electronic_levels, modes, si_scales);
            const canary::ModeIndex n2_mode {reorder_modes ? 1U : 0U};
            const canary::ModeIndex o2_mode {reorder_modes ? 0U : 1U};
            for (const bool reorder_reactants : {false, true}) {
                auto input = dissociation;
                if (reorder_reactants) {
                    std::swap(input.reactants.entries[0], input.reactants.entries[1]);
                }
                const auto error = [&](std::string_view role) {
                    const std::string message =
                      error_message([&] { (void)build(input, reordered.model); });
                    expect(message == "reaction 0 (N2 + O2 <=> N + N + O2): " +
                                        std::string(role) + " rate: reactant vibrations occupy "
                                        "different modes; specify other_temperature") << message;
                };
                error("forward");
                // The standard three-temperature model shares one vibration mode.
                expect(build(input, gas).forward.other_mode == canary::ModeIndex {0});
                input.forward = rate(1.0, nullptr);
                expect(build(input, reordered.model).forward.other_mode == canary::ModeIndex {0});
                input.forward = rate(0.5, "N2");
                expect(build(input, reordered.model).forward.other_mode == n2_mode);
                input.forward = rate(0.5, "O2");
                expect(build(input, reordered.model).forward.other_mode == o2_mode);
                input.forward = rate(1.0, "O2");
                expect(build(input, reordered.model).forward.other_mode == o2_mode);
                input.backward_model = BackwardRateModel::arrhenius;
                input.backward = rate(0.5, nullptr);
                error("backward");
                input.backward = rate(0.5, "N2");
                expect(build(input, reordered.model).backward.other_mode == n2_mode);
                input.backward = rate(1.0, nullptr);
                expect(build(input, reordered.model).backward.other_mode == canary::ModeIndex {0});
            }
        }
        // A species name: the mode of its vibration; electron: the electron
        // mode.
        expect(build(o2(rate(0.5, "N2")), split.model).forward.other_mode == canary::ModeIndex {0});
        expect(build(o2(rate(0.5, "NO")), split.model).forward.other_mode == canary::ModeIndex {1});
        expect(build(o2(rate(0.5, "electron")), split.model).forward.other_mode ==
               canary::ModeIndex {2});
        // s = 1 without a vibrating reactant is accepted.
        expect(build(ionise(rate(1.0, nullptr)), split.model).forward.other_mode ==
               canary::ModeIndex {0});
        // The backward Arrhenius rate resolves with the reactant list; an
        // equilibrium backward rate only with an other temperature.
        const auto arrhenius =
          build(o2(rate(1.0, nullptr), BackwardRateModel::arrhenius, rate(0.5, nullptr)),
                split.model);
        expect(arrhenius.backward.other_mode == canary::ModeIndex {1} &&
               !arrhenius.equilibrium_at_other_temperature);
        const auto at_electron =
          build(o2(rate(1.0, nullptr), BackwardRateModel::equilibrium, rate(1.0, "electron")),
                split.model);
        expect(at_electron.equilibrium_at_other_temperature &&
               at_electron.backward.other_mode == canary::ModeIndex {2});
        expect(!build(o2(rate(1.0, nullptr)), split.model).equilibrium_at_other_temperature);
        // One mode: every rate resolves to it.
        const auto one = canary::make_multi_temperature_mixture(canary::air11_species,
                                                                canary::air11_rotational_degrees,
                                                                canary::air11_electronic_levels,
                                                                single_mode,
                                                                si_scales);
        expect(build(ionise(rate(0.5, nullptr)), one.model).forward.other_mode ==
               canary::ModeIndex {0});
        expect(build(o2(rate(0.5, "electron")), one.model).forward.other_mode ==
               canary::ModeIndex {0});

        const auto fails = [&](const canary::ReactionInput<3>& r, std::string_view message) {
            const std::string got = error_message([&] { (void)build(r, split.model); });
            expect(got == message) << got;
        };
        fails(o2(rate(0.5, "Ar")),
              "reaction 0 (O + O2 <=> O + O + O): forward rate: other temperature \"Ar\" is not "
              "\"electron\" or a species");
        fails(o2(rate(0.5, "N")),
              "reaction 0 (O + O2 <=> O + O + O): forward rate: other temperature \"N\": the "
              "species has no vibration component");
        // The vibration of N2+ is frozen in this model.
        fails(o2(rate(1.0, "N2+")),
              "reaction 0 (O + O2 <=> O + O + O): forward rate: other temperature \"N2+\": the "
              "species has no vibration component");
        fails(ionise(rate(0.5, nullptr)),
              "reaction 0 (N + O <=> NO+ + e-): forward rate: no reactant has a vibration "
              "component to give the other temperature");
        fails(o2(rate(1.0, nullptr), BackwardRateModel::equilibrium, rate(1.0, "Ar")),
              "reaction 0 (O + O2 <=> O + O + O): backward rate: other temperature \"Ar\" is "
              "not \"electron\" or a species");

        // The two-temperature model has one mode holding every species and
        // no free electrons.
        const auto air5 = canary::make_air5_two_temperature<double>(si_scales);
        const auto build5 = [&](canary::ArrheniusInput f) {
            const canary::MechanismInput<1, 3> input {
              .temperature_lower = 300.0,
              .temperature_upper = 30000.0,
              .reactions =
                {
                  reaction("N + O <=> NO",
                           coefficients({{.name = "N", .value = 1.0}, {.name = "O", .value = 1.0}}),
                           coefficients({{.name = "NO", .value = 1.0}}),
                           f,
                           canary::BackwardRateModel::equilibrium,
                           none),
                },
            };
            return canary::make_reaction_mechanism(input,
                                                   canary::reaction_species(canary::air5_species),
                                                   si_scales,
                                                   p0,
                                                   air5.model)
              .reactions[0];
        };
        expect(build5(rate(0.5, "N")).forward.other_mode == canary::ModeIndex {0});
        expect(build5(rate(0.5, nullptr)).forward.other_mode == canary::ModeIndex {0});
        const std::string got = error_message([&] { (void)build5(rate(0.5, "electron")); });
        expect(got ==
               "reaction 0 (N + O <=> NO): forward rate: other temperature \"electron\": the model "
               "has no free electrons")
          << got;
    };

    "multi-temperature kinetics: translational rates ignore other mode logarithms"_test = [] {
        const canary::ArrheniusRate<double, si> rate {
          .A = 1.0, .n = 2.0, .activation_temperature = Temp {3.0},
          .translational_exponent = 1.0, .other_mode = canary::ModeIndex {0}};
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const std::array<double, 2> logarithms {nan, nan};
        const double value = canary::detail::rate_exponent<2>(rate, 4.0, 0.5, logarithms);
        expect(value == 6.5_d);
        const std::array<D, 2> dual_logarithms {D {nan, nan}, D {nan, nan}};
        const D dual = canary::detail::rate_exponent<2>(
          rate, D {4.0, 1.0}, D {0.5, -0.25}, dual_logarithms);
        expect(dual.real() == 6.5_d);
        expect(dual.dual() == 2.75_d);
    };

    "multi-temperature kinetics: equilibrium backward rates at T_e depend on T_e only"_test = [&] {
        const ReferenceState& s = reference_states[0];
        const auto mix = canary::mixture(gas, fractions(s.y));
        const auto rates = [&](double temp, double t_v, double t_e) {
            const std::array t_m {
              canary::ModeTemperature<double, si> {t_v},
              canary::ModeTemperature<double, si> {t_e},
            };
            return canary::rates_of_progress(mech, gas, mix, Rho {s.rho}, Temp {temp}, t_m);
        };
        const auto a = rates(12000.0, 5000.0, 8000.0);
        for (const auto& [temp, t_v] : {
               std::array {9000.0, 3000.0},
               std::array {250.0, 20000.0},
               std::array {40000.0, 8000.0},
             }) {
            const auto b = rates(temp, t_v, 8000.0);
            for (const std::size_t k : backward_at_electron) {
                expect(a.backward[k] == b.backward[k]) << temp << t_v << k;
            }
        }
        // T_e is not clamped: 40000 K is above the upper limit 30000 K.
        const auto hot = rates(12000.0, 5000.0, 40000.0);
        const auto limit = rates(12000.0, 5000.0, 30000.0);
        expect(!(hot.backward[38] == limit.backward[38]));
        const std::array<D, species_count> y_d = [&] {
            std::array<D, species_count> y {};
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = D(s.y[i]); }
            return y;
        }();
        const auto mix_d = canary::mixture(gas, fractions(y_d));
        for (const std::size_t column : std::array<std::size_t, 3> {0, 1, 2}) {
            const std::array t_m {
              canary::ModeTemperature<D, si> {D(5000.0, column == 1 ? 1.0 : 0.0)},
              canary::ModeTemperature<D, si> {D(8000.0, column == 2 ? 1.0 : 0.0)},
            };
            const auto q = canary::rates_of_progress(
              mech,
              gas,
              mix_d,
              canary::Density<D, si> {D(s.rho)},
              canary::Temperature<D, si> {D(12000.0, column == 0 ? 1.0 : 0.0)},
              t_m);
            for (const std::size_t k : backward_at_electron) {
                const double dual = q.backward[k].value().dual();
                expect(column == 2 ? dual != 0.0 : dual == 0.0) << column << k << dual;
            }
        }
    };

    "multi-temperature kinetics: exchanges vanish at equal temperatures and V-E is "
    "antisymmetric"_test = [&] {
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto at = [&](double temp, double t_v, double t_e) {
                return cell_state(gas, s.y, s.rho, temp, t_v, t_e);
            };
            const auto et =
              canary::electron_translation_source(exchange, mix, at(s.temp, s.t_v, s.temp));
            expect(et[0].value() == 0.0_d && et[1].value() == 0.0_d) << s.name;
            const auto ve =
              canary::vibration_electron_source(exchange, gas, mix, at(s.temp, s.t_v, s.t_v));
            expect(ve[0].value() == 0.0_d && ve[1].value() == 0.0_d) << s.name;
            const auto vt =
              canary::vibration_translation_source(exchange, gas, mix, at(s.temp, s.temp, s.t_e));
            expect(vt[0].value() == 0.0_d && vt[1].value() == 0.0_d) << s.name;
            const auto state = at(s.temp, s.t_v, s.t_e);
            const auto ve_state = canary::vibration_electron_source(exchange, gas, mix, state);
            expect(ve_state[0].value() == -ve_state[1].value()) << s.name;
            // The sum over the modes of the total is that of Q_VT + Q_CV +
            // Q_ET + Q_CE; the terms are summed in another order.
            const Values got = evaluate(mech, exchange, gas, s);
            double sum = 0;
            double expected = 0;
            double scale = 0;
            for (const canary::ModeIndex m : canary::mode_indices(0, mode_count)) {
                sum += got.total[m.value()];
                expected +=
                  got.vt[m.value()] + got.cv[m.value()] + got.et[m.value()] + got.ce[m.value()];
                scale += std::abs(got.vt[m.value()]) + std::abs(got.cv[m.value()]) +
                         std::abs(got.et[m.value()]) + std::abs(got.ve[m.value()]) +
                         std::abs(got.ce[m.value()]);
            }
            expect(std::abs(sum - expected) <= 1.0e-14 * scale) << s.name << sum << expected;
            // V-T has the sign of T - T_v in the vibration mode.
            if (s.temp != s.t_v) {
                expect(s.temp > s.t_v ? got.vt[0] > 0.0 : got.vt[0] < 0.0) << s.name;
            }
        }
    };

    "multi-temperature kinetics: cross-sections and the Candler fit follow their formulas"_test =
      [&] {
          const canary::EnergyExchangeInput<0, 0> no_vibration {};
          canary::ElectronExchangeInput<2, 0, 0> electrons {};
          electrons.electron_translation[0] =
            canary::air11_kim_jo_exchange_electrons.electron_translation[0];
          electrons.electron_translation[1] =
            canary::air11_kim_jo_exchange_electrons.electron_translation[5];
          const auto pairs = canary::make_multi_temperature_exchange(no_vibration,
                                                                     electrons,
                                                                     names,
                                                                     built.species,
                                                                     mech,
                                                                     gas,
                                                                     si_scales);
          const double e = canary::constants::elementary_charge<double>;
          const double eps_0 = 8.8541878128e-12;
          const double m_e = built.species.mass[electron].value();
          for (const double y_e : {9.203080774345904e-07, 1.0e-40}) {
              Fractions y = reference_states[0].y;
              y[electron] = y_e;
              const double rho = 1.0e-3;
              for (const double t_e : {8000.0, 45000.0}) {
                  const double temp = 12000.0;
                  const auto mix = canary::mixture(gas, fractions(y));
                  const auto state = cell_state(gas, y, rho, temp, 5000.0, t_e);
                  const double got =
                    canary::electron_translation_source(pairs, mix, state)[1].value();
                  // The definitions, with the normalised fractions of mix.
                  Fractions y_n {};
                  for (std::size_t i = 0; i < species_count; ++i) {
                      y_n[i] = mix.fractions[i].value();
                  }
                  const Fractions n = number_densities(built.species, y_n, rho);
                  const double t_n2 = std::min(t_e, 30000.0);
                  const double sigma_n2 = 7.5e-20 + (5.5e-24 * t_n2) - (1.0e-28 * t_n2 * t_n2);
                  const double length = e * e / (4.0 * std::numbers::pi * eps_0 * boltzmann * t_e);
                  const double n_e = std::max(n[electron], 1.0);
                  const double sigma_coulomb =
                    (2.0 * std::numbers::pi / 9.0) * length * length *
                    std::log(9.0 / (4.0 * std::numbers::pi * n_e * length * length * length));
                  const double c_e = std::sqrt(8.0 * boltzmann * t_e / (std::numbers::pi * m_e));
                  const double sum = (m_e / built.species.mass[0].value() * n[0] * sigma_n2) +
                                     (m_e / built.species.mass[5].value() * n[5] * sigma_coulomb);
                  const double expected = 3.0 * n[electron] * boltzmann * (temp - t_e) * c_e * sum;
                  // The same expressions in another order: a few ulp, and the
                  // logarithm of Lambda (about 20) loses no more.
                  expect(relative(got, expected, 1.0e-13)) << y_e << t_e << got << expected;
              }
          }

          // The blend weight is 1 at blend_lower and 0 at blend_upper, so the
          // fit is continuous at both ends.
          const auto& fit = exchange.vibration_electron[0].relaxation_time;
          const auto f = [&](double t) {
              return canary::detail::candler_exponent(fit, std::log10(t), t);
          };
          // The input fit in log10(p_e tau / (1 atm s)), and the model fit
          // in SI, in log10(p_e tau / (1 Pa s)).
          const auto quadratic = [](double a, double b, double c, double t) {
              const double x = std::log10(t);
              return (((a * x) + b) * x) + c;
          };
          const auto model_fit = [&](double a, double b, double c, double t) {
              return quadratic(a, b, c + std::log10(atmosphere), t);
          };
          for (const double lower : {6950.0, std::nextafter(6950.0, 1.0e5)}) {
              // At the bound w = 1; one ulp above, 1 - w is of order 1e-30.
              expect(relative(f(lower), model_fit(7.5, -57.0, 98.7, lower), 1.0e-14)) << lower;
          }
          for (const double upper : {7050.0, std::nextafter(7050.0, 0.0)}) {
              expect(relative(f(upper), model_fit(2.36, -17.9, 24.35, upper), 1.0e-14)) << upper;
          }
          // Between the bounds, w f_low + (1 - w) f_high.
          const double w = std::pow(std::cos(std::numbers::pi / 2.0 * 0.3), 2);
          expect(relative(f(6980.0),
                          (w * model_fit(7.5, -57.0, 98.7, 6980.0)) +
                            ((1.0 - w) * model_fit(2.36, -17.9, 24.35, 6980.0)),
                          1.0e-13));

          // V-E with p_e = 0: no free electrons give exactly zero.
          Fractions y = reference_states[0].y;
          y[electron] = 0.0;
          const auto mix = canary::mixture(gas, fractions(y));
          const auto ve =
            canary::vibration_electron_source(exchange,
                                              gas,
                                              mix,
                                              cell_state(gas, y, 1.0e-3, 12000.0, 5000.0, 8000.0));
          expect(ve[0].value() == 0.0_d && ve[1].value() == 0.0_d);
          // V-E of N2 from Candler's fit above the blend range.
          const ReferenceState& s = reference_states[0];
          const auto mix_s = canary::mixture(gas, fractions(s.y));
          const auto state = cell_state(gas, s.y, s.rho, s.temp, s.t_v, s.t_e);
          const double q =
            canary::vibration_electron_source(exchange, gas, mix_s, state)[0].value();
          const double n_e =
            s.rho * mix_s.fractions[electron].value() / built.species.mass[electron].value();
          const double p_e = n_e * boltzmann * s.t_e;
          const double inverse_tau =
            (p_e / atmosphere) * std::pow(10.0, -quadratic(2.36, -17.9, 24.35, s.t_e));
          const auto e_te = canary::species_mode_energy(gas, canary::ModeIndex {0}, Temp {s.t_e});
          const auto e_tv = canary::species_mode_energy(gas, canary::ModeIndex {0}, Temp {s.t_v});
          const double expected =
            s.rho * mix_s.fractions[0].value() * (e_te[0].value() - e_tv[0].value()) * inverse_tau;
          // 10^(-f) from exp(-f ln 10) with f of about -10: a few ulp.
          expect(relative(q, expected, 1.0e-13)) << q << expected;
      };

    "multi-temperature kinetics: net ionisation cools the electron mode"_test = [&] {
        // Weakly ionised air with a hot electron mode ionises; strongly
        // ionised air with a cold one recombines.
        for (const std::size_t index : {std::size_t {0}, std::size_t {3}}) {
            const ReferenceState& s = reference_states[index];
            const double t_e = index == 0 ? 15000.0 : 3000.0;
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto state = cell_state(gas, s.y, s.rho, s.temp, s.t_v, t_e);
            const auto q =
              canary::rates_of_progress(mech, gas, mix, state.rho, state.temp, state.mode_temp);
            const double net = q.forward[38].value() - q.backward[38].value() +
                               q.forward[39].value() - q.backward[39].value();
            const double ce = canary::electron_impact_ionisation_source(exchange, q)[1].value();
            expect(index == 0 ? net > 0.0 : net < 0.0) << index << net;
            expect(index == 0 ? ce < 0.0 : ce > 0.0) << index << ce;
        }
    };

    "multi-temperature kinetics: production rates conserve mass and charge"_test = [&] {
        const auto check =
          [&](const Fractions& y, double rho, double temp, double t_v, double t_e) {
              const auto mix = canary::mixture(gas, fractions(y));
              const auto state = cell_state(gas, y, rho, temp, t_v, t_e);
              const auto q =
                canary::rates_of_progress(mech, gas, mix, state.rho, state.temp, state.mode_temp);
              const auto omega =
                canary::production_rates(mech, gas, mix, state.rho, state.temp, state.mode_temp);
              double mass = 0;
              double mass_scale = 0;
              double charge = 0;
              double charge_scale = 0;
              for (std::size_t i = 0; i < species_count; ++i) {
                  const double m = built.species.mass[i].value();
                  const double z = built.species.charge[i].value();
                  mass += omega[i].value();
                  charge += z * omega[i].value() / m;
                  // The terms omega_i is summed from.
                  for (std::size_t k = 0; k < reaction_count; ++k) {
                      const double both = q.forward[k].value() + q.backward[k].value();
                      mass_scale += m * std::abs(coefficient_change(mech, k, i)) * both;
                      charge_scale += std::abs(z * coefficient_change(mech, k, i)) * both;
                  }
              }
              // Sums of terms of opposite signs, in double, relative to the
              // terms. The molar masses of the species data balance the
              // ionisations of reactions 24 to 26, 38 and 39 only to 2e-11
              // kg/mol (an ion's mass is that of its neutral less that of e-
              // to the 1e-10 kg/mol given), up to 7.2e-10 of the masses that
              // the reaction converts.
              expect(std::abs(mass) <= 1.0e-9 * mass_scale) << temp << mass << mass_scale;
              expect(std::abs(charge) <= 1.0e-14 * charge_scale) << temp << charge << charge_scale;
          };
        for (const ReferenceState& s : reference_states) {
            check(s.y, s.rho, s.temp, s.t_v, s.t_e);
        }
        // A fixed seed gives the same states in every run.
        std::mt19937 engine(1357);  // NOLINT(bugprone-random-generator-seed)
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (int n = 0; n < 50; ++n) {
            Fractions y {};
            double total = 0;
            for (double& v : y) {
                v = unit(engine);
                total += v;
            }
            for (double& v : y) { v /= total; }
            const auto neutral = canary::charge_neutral_fractions(gas, built.species, fractions(y));
            expect(neutral.status == canary::Status::Valid);
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = neutral.fractions[i].value(); }
            check(y,
                  std::pow(10.0, -5.0 + (3.0 * unit(engine))),
                  1000.0 + (20000.0 * unit(engine)),
                  1000.0 + (20000.0 * unit(engine)),
                  1000.0 + (20000.0 * unit(engine)));
        }
    };

    "multi-temperature kinetics: dual-number derivatives of the mode sources equal central "
    "differences"_test = [&] {
        struct Point {
            Fractions y;
            double rho;
            double temp;
            double t_v;
            double t_e;
        };
        const ReferenceState& a = reference_states[0];
        const ReferenceState& b = reference_states[3];
        // Two reference states, and T_e in the blend range of Candler's fit.
        const std::array<Point, 3> points {
          {
            {.y = a.y, .rho = a.rho, .temp = a.temp, .t_v = a.t_v, .t_e = a.t_e},
            {.y = b.y, .rho = b.rho, .temp = b.temp, .t_v = b.t_v, .t_e = b.t_e},
            {.y = a.y, .rho = 1.0e-3, .temp = 9000.0, .t_v = 6000.0, .t_e = 7000.0},
          },
        };
        constexpr std::size_t outputs = 8;
        // The sources V-T (vibration), C-V (vibration), E-T (electron), V-E
        // (both), C-E (electron) and the totals at x = (T, T_v, T_e, y_e),
        // then the size of the terms each is summed from.
        const auto sources = [&]<typename T>(const Point& p, const std::array<T, 4>& x) {
            using canary::ScalarTraits;
            const auto magnitude = [](T v) { return v < T(0.0) ? T(-v) : v; };
            std::array<T, species_count> y {};
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = T(p.y[i]); }
            y[electron] = x[3];
            const auto mix = canary::mixture(gas, fractions(y));
            const auto state = cell_state(gas, y, T(p.rho), x[0], x[1], x[2]);
            const auto q =
              canary::rates_of_progress(mech, gas, mix, state.rho, state.temp, state.mode_temp);
            const auto vt = canary::vibration_translation_source(exchange, gas, mix, state);
            const auto cv = canary::chemistry_vibration_source(exchange, mech, q, state);
            const auto et = canary::electron_translation_source(exchange, mix, state);
            const auto ve = canary::vibration_electron_source(exchange, gas, mix, state);
            const auto ce = canary::electron_impact_ionisation_source(exchange, q);
            const auto total = canary::multi_temperature_sources(mech, exchange, gas, mix, state);
            std::array<T, 2 * outputs> out {
              vt[0].value(),
              cv[0].value(),
              et[1].value(),
              ve[0].value(),
              ve[1].value(),
              ce[1].value(),
              total.modes[0].value(),
              total.modes[1].value(),
            };
            T cv_scale = T(0.0);
            for (const auto& c : exchange.chemistry_vibration) {
                cv_scale =
                  cv_scale + T(boltzmann * c.dissociation_temperature.value() *
                               (q.forward[c.reaction].value() + q.backward[c.reaction].value()));
            }
            T ce_scale = T(0.0);
            for (const auto& c : exchange.electron_impact_ionisation) {
                ce_scale = ce_scale + T(c.ionisation_energy * (q.forward[c.reaction].value() +
                                                               q.backward[c.reaction].value()));
            }
            out[outputs] = magnitude(out[0]);
            out[outputs + 1] = cv_scale;
            out[outputs + 2] = magnitude(out[2]);
            out[outputs + 3] = magnitude(out[3]);
            out[outputs + 4] = magnitude(out[4]);
            out[outputs + 5] = ce_scale;
            out[outputs + 6] = out[outputs] + cv_scale + out[outputs + 3];
            out[outputs + 7] = out[outputs + 2] + out[outputs + 4] + ce_scale;
            return out;
        };
        for (const Point& p : points) {
            const std::array<double, 4> x {p.temp, p.t_v, p.t_e, p.y[electron]};
            for (std::size_t c = 0; c < 4; ++c) {
                std::array<D, 4> x_d {};
                std::array<double, 4> plus = x;
                std::array<double, 4> minus = x;
                const double h = 1.0e-5 * x[c];
                for (std::size_t k = 0; k < 4; ++k) { x_d[k] = D(x[k], k == c ? 1.0 : 0.0); }
                plus[c] += h;
                minus[c] -= h;
                const auto ad = sources(p, x_d);
                const auto up = sources(p, plus);
                const auto down = sources(p, minus);
                for (std::size_t i = 0; i < outputs; ++i) {
                    const double fd = (up[i] - down[i]) / (2.0 * h);
                    // Central differences with h = 1e-5 x have truncation
                    // error of order 1e-10 relative, times the square of
                    // theta / T (up to about 30) in the rates. Their rounding
                    // error is about eps times the size of the terms summed
                    // in the value over 2 h; 1e-13 allows for 400 ulp of the
                    // terms.
                    const double rounding = 1.0e-13 * up[outputs + i] / h;
                    expect(std::abs(fd - ad[i].dual()) <=
                           (1.0e-6 * std::abs(ad[i].dual())) + rounding)
                      << p.temp << p.t_e << c << i << ad[i].dual() << fd;
                }
            }
        }
    };

    "multi-temperature kinetics: scaled units give the SI results after conversion"_test = [&] {
        const auto built_s = canary::make_air11_three_temperature<double>(reference_scales);
        const auto& gas_s = built_s.model;
        const auto mech_s = canary::make_air11_kim_jo(reference_scales,
                                                      canary::convert<scaled>(p0, reference_scales),
                                                      gas_s);
        const auto exchange_s =
          canary::make_air11_kim_jo_exchange(built_s.species, mech_s, gas_s, reference_scales);
        const auto to_si = [&](const auto& q) {
            return canary::convert<si>(q, reference_scales).value();
        };
        const auto in_scaled = [&](double t) {
            return canary::convert<scaled>(Temp {t}, reference_scales).value();
        };
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto mix_s = canary::mixture(gas_s, fractions(s.y));
            const auto state = cell_state(gas, s.y, s.rho, s.temp, s.t_v, s.t_e);
            const auto state_s =
              cell_state(gas_s,
                         s.y,
                         canary::convert<scaled>(Rho {s.rho}, reference_scales).value(),
                         in_scaled(s.temp),
                         in_scaled(s.t_v),
                         in_scaled(s.t_e));
            const auto q =
              canary::rates_of_progress(mech, gas, mix, state.rho, state.temp, state.mode_temp);
            const auto q_s = canary::rates_of_progress(mech_s,
                                                       gas_s,
                                                       mix_s,
                                                       state_s.rho,
                                                       state_s.temp,
                                                       state_s.mode_temp);
            // As for one temperature: the Gibbs energies of the scaled thermo
            // model differ from the SI ones by up to 6e-13 relative, times
            // m / (k_B T) in K_c.
            for (std::size_t k = 0; k < reaction_count; ++k) {
                expect(relative(to_si(q_s.forward[k]), q.forward[k].value(), 1.0e-12))
                  << s.name << k;
                expect(relative(to_si(q_s.backward[k]), q.backward[k].value(), 1.0e-10))
                  << s.name << k;
            }
            // The SI rates converted, so that only the exchange is compared.
            canary::RatesOfProgress<double, scaled, reaction_count> q_converted {};
            for (std::size_t k = 0; k < reaction_count; ++k) {
                q_converted.forward[k] = canary::convert<scaled>(q.forward[k], reference_scales);
                q_converted.backward[k] = canary::convert<scaled>(q.backward[k], reference_scales);
            }
            const auto cv = canary::chemistry_vibration_source(exchange, mech, q, state);
            const auto cv_s =
              canary::chemistry_vibration_source(exchange_s, mech_s, q_converted, state_s);
            const auto ce = canary::electron_impact_ionisation_source(exchange, q);
            const auto ce_s = canary::electron_impact_ionisation_source(exchange_s, q_converted);
            // Sums of terms of one sign, with constants converted in double.
            expect(relative(to_si(cv_s[0]), cv[0].value(), 1.0e-12)) << s.name;
            // The net rates q_f - q_b cancel by up to a factor 100 here.
            expect(relative(to_si(ce_s[1]), ce[1].value(), 1.0e-11)) << s.name;
            // The differences of energies of the exchanges cancel when the
            // temperatures are close; states where they differ by 1000 K or
            // more. The Candler fit of V-E is shifted by log10 of the scales
            // (up to about 170 times the coefficients' rounding).
            const auto vt = canary::vibration_translation_source(exchange, gas, mix, state);
            const auto vt_s =
              canary::vibration_translation_source(exchange_s, gas_s, mix_s, state_s);
            const auto et = canary::electron_translation_source(exchange, mix, state);
            const auto et_s = canary::electron_translation_source(exchange_s, mix_s, state_s);
            const auto ve = canary::vibration_electron_source(exchange, gas, mix, state);
            const auto ve_s = canary::vibration_electron_source(exchange_s, gas_s, mix_s, state_s);
            if (std::abs(s.temp - s.t_v) >= 1000.0) {
                expect(relative(to_si(vt_s[0]), vt[0].value(), 1.0e-12)) << s.name;
            }
            if (std::abs(s.temp - s.t_e) >= 1000.0) {
                expect(relative(to_si(et_s[1]), et[1].value(), 1.0e-12)) << s.name;
            }
            if (std::abs(s.t_e - s.t_v) >= 1000.0) {
                expect(relative(to_si(ve_s[0]), ve[0].value(), 1.0e-11)) << s.name;
            }
        }
    };

    "multi-temperature kinetics: neutral cross-sections are non-negative throughout the interval"_test =
      [&] {
          struct Case {
              double a;
              double b;
              double c;
              bool valid;
          };
          for (const Case value : std::array {
                 Case {-1.0, 0.0, 0.0, false},
                 Case {1.0, -1.0, 0.0, false},
                 Case {0.75, -2.0, 1.0, false},
                 Case {1.0, -2.0, 1.0, true},
                 Case {0.0, 0.0, 0.0, true},
                 Case {0.0, 2.0, -1.0, true},
                 Case {9.0, -6.0, 1.0, true},
               }) {
              auto electrons = canary::air11_kim_jo_exchange_electrons;
              auto& pair = electrons.electron_translation[0];
              pair.a = value.a;
              pair.b = value.b;
              pair.c = value.c;
              pair.maximum_temperature = 2.0;
              const auto message = error_message([&] {
                  (void)canary::make_multi_temperature_exchange(canary::air11_kim_jo_exchange,
                                                                electrons,
                                                                names,
                                                                built.species,
                                                                mech,
                                                                gas,
                                                                si_scales);
              });
              if (value.valid) {
                  expect(message == "no exception") << message;
              } else {
                  expect(message.starts_with("electron-translation pair 0 (N2): neutral "
                                             "cross-section must be finite and non-negative"))
                    << message;
              }
          }
      };

    "multi-temperature kinetics: float rounding cannot introduce a negative cross-section minimum"_test =
      [] {
          const std::array<std::string_view, 2> pair_names {"N2", "carrier"};
          const canary::SpeciesData<float, si, 2> species {
            .n = 2,
            .mass = {canary::ParticleMass<float, si> {4.65e-26F},
                     canary::ParticleMass<float, si> {9.109e-31F}},
            .charge = {canary::Charge<float, si> {0.0F}, canary::Charge<float, si> {-1.602e-19F}}};
          const double root = 1.0 + 0.75 * std::numeric_limits<float>::epsilon();
          const canary::ElectronTranslationInput input {
            .partner = "N2",
            .cross_section = canary::ElectronCrossSection::gnoffo_neutral,
            .a = root * root + 1.0e-12,
            .b = -2.0 * root,
            .c = 1.0,
            .maximum_temperature = 2.0};
          // The input minimum is positive; rounding a and b to float makes it negative.
          expect(std::fma(-root, root, input.a) > 0.0);
          expect(std::fma(-0.5 * double(float(input.b)),
                          0.5 * double(float(input.b)),
                          double(float(input.a))) < 0.0);
          expect(error_message([&] {
                     (void)canary::detail::make_electron_translation_pair(
                       input,
                       pair_names,
                       species,
                       canary::SpeciesIndex {1},
                       canary::UnitScales<float, si> {});
                 }).starts_with("neutral cross-section must be finite and non-negative"));
      };

    "multi-temperature kinetics: the builder rejects invalid input with the entry in the "
    "message"_test = [&] {
        using Input = std::remove_cvref_t<decltype(canary::air11_kim_jo_exchange)>;
        using Electrons = std::remove_cvref_t<decltype(canary::air11_kim_jo_exchange_electrons)>;
        const auto build = [&](const Input& input, const Electrons& electrons) {
            (void)canary::make_multi_temperature_exchange(input,
                                                          electrons,
                                                          names,
                                                          built.species,
                                                          mech,
                                                          gas,
                                                          si_scales);
        };
        const auto starts =
          [&](const Input& input, const Electrons& electrons, std::string_view prefix) {
              const std::string message = error_message([&] { build(input, electrons); });
              expect(message.starts_with(prefix)) << message;
          };
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        const Input base = canary::air11_kim_jo_exchange;
        const Electrons electrons = canary::air11_kim_jo_exchange_electrons;
        build(base, electrons);

        Input input = base;
        input.vibration_translation[3].species = "N";
        starts(input,
               electrons,
               "vibration-translation pair 3 (N, N): species \"N\" has no vibration component");
        input = base;
        input.vibration_translation[4].millikan_white.a = nan;
        starts(input, electrons, "vibration-translation pair 4 (N2, O): Millikan-White a must be");
        input = base;
        input.chemistry_vibration[0].species = "N";
        starts(input,
               electrons,
               "chemistry-vibration coupling 0 (N): species \"N\" has no vibration component");

        Electrons e = electrons;
        e.electron_translation[2].partner = "Ar";
        starts(base, e, "electron-translation pair 2 (Ar): unknown partner \"Ar\"");
        e = electrons;
        e.electron_translation[0].partner = "e-";
        starts(base, e, "electron-translation pair 0 (e-): the partner must be a heavy species");
        e = electrons;
        e.electron_translation[5].partner = "N";
        starts(base, e, "electron-translation pair 5 (N): a coulomb partner must be charged");
        e = electrons;
        e.electron_translation[0].partner = "N+";
        starts(base,
               e,
               "electron-translation pair 0 (N+): a gnoffo_neutral partner must be neutral");
        e = electrons;
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange): an invalid model
        e.electron_translation[0].cross_section = static_cast<canary::ElectronCrossSection>(5);
        starts(base, e, "electron-translation pair 0 (N2): unknown electron cross-section model");
        for (const double v : {nan, inf}) {
            e = electrons;
            e.electron_translation[1].a = v;
            starts(base,
                   e,
                   "electron-translation pair 1 (O2): cross-section coefficient a must be");
            e = electrons;
            e.electron_translation[1].b = v;
            starts(base,
                   e,
                   "electron-translation pair 1 (O2): cross-section coefficient b must be");
            e = electrons;
            e.electron_translation[1].c = v;
            starts(base,
                   e,
                   "electron-translation pair 1 (O2): cross-section coefficient c must be");
        }
        for (const double t_max : {nan, inf, 0.0, -1.0}) {
            e = electrons;
            e.electron_translation[3].maximum_temperature = t_max;
            starts(base,
                   e,
                   "electron-translation pair 3 (O): cross-section maximum temperature must be");
        }
        // A coulomb pair ignores the coefficients.
        e = electrons;
        e.electron_translation[6].a = nan;
        e.electron_translation[6].maximum_temperature = nan;
        build(base, e);

        e = electrons;
        e.vibration_electron[0].species = "Ar";
        starts(base, e, "vibration-electron exchange 0 (Ar): unknown species \"Ar\"");
        e = electrons;
        e.vibration_electron[0].species = "O";
        starts(base, e, "vibration-electron exchange 0 (O): species \"O\" has no vibration");
        for (const double v : {nan, inf}) {
            e = electrons;
            e.vibration_electron[0].relaxation_time.b_high = v;
            starts(base, e, "vibration-electron exchange 0 (N2): relaxation-time parameters");
        }
        for (const double upper : {6950.0, 6000.0}) {
            e = electrons;
            e.vibration_electron[0].relaxation_time.blend_upper = upper;
            starts(base,
                   e,
                   "vibration-electron exchange 0 (N2): relaxation-time blend bounds must satisfy");
        }

        for (const int reaction : {-1, 40}) {
            e = electrons;
            e.electron_impact_ionisation[0].reaction = reaction;
            starts(base,
                   e,
                   "electron-impact ionisation 0 (O): reaction index must be between 0 and 39");
        }
        e = electrons;
        e.electron_impact_ionisation[0].species = "Ar";
        starts(base, e, "electron-impact ionisation 0 (Ar): unknown species \"Ar\"");
        e = electrons;
        e.electron_impact_ionisation[0].species = "e-";
        starts(base,
               e,
               "electron-impact ionisation 0 (e-): the ionised species must be a heavy species");
        // Reaction 38 ionises O, not N.
        e = electrons;
        e.electron_impact_ionisation[0].species = "N";
        starts(base,
               e,
               "electron-impact ionisation 0 (N): reaction 38 must consume one \"N\" (net change -1)");
        for (const double energy : {nan, inf, 0.0, -1.0}) {
            e = electrons;
            e.electron_impact_ionisation[1].ionisation_energy = energy;
            starts(base, e, "electron-impact ionisation 1 (N): ionisation energy must be");
        }

        auto species = built.species;
        species.n = 10;
        expect(error_message([&] {
                   (void)canary::make_air11_kim_jo_exchange(species, mech, gas, si_scales);
               }).starts_with("species data"));

        // Five-species air has no free electrons.
        const auto air5 = canary::make_air5_four_temperature<double>(si_scales);
        const auto air5_mech = canary::make_air5_park(si_scales, p0, air5.model);
        const std::array<std::string_view, 5> air5_names {"N2", "O2", "N", "O", "NO"};
        const canary::EnergyExchangeInput<0, 0> nothing {};
        const auto without_electrons = [&](const auto& extra) {
            return error_message([&] {
                (void)canary::make_multi_temperature_exchange(nothing,
                                                              extra,
                                                              air5_names,
                                                              air5.species,
                                                              air5_mech,
                                                              air5.model,
                                                              si_scales);
            });
        };
        canary::ElectronExchangeInput<1, 0, 0> pair {};
        pair.electron_translation[0] = electrons.electron_translation[0];
        expect(without_electrons(pair) ==
               "electron-translation pair 0 (N2): the model has no free electrons");
        canary::ElectronExchangeInput<0, 1, 0> relaxation {};
        relaxation.vibration_electron[0] = electrons.vibration_electron[0];
        expect(without_electrons(relaxation) ==
               "vibration-electron exchange 0 (N2): the model has no free electrons");
        canary::ElectronExchangeInput<0, 0, 1> ionisation {};
        ionisation.electron_impact_ionisation[0] = {
          .reaction = 0,
          .species = "N",
          .ionisation_energy = 1.4e6,
        };
        expect(without_electrons(ionisation) ==
               "electron-impact ionisation 0 (N): the model has no free electrons");
        // Without electron exchanges the model needs no free electrons.
        const canary::ElectronExchangeInput<0, 0, 0> no_electrons {};
        const auto vibration_only = canary::make_multi_temperature_exchange(nothing,
                                                                            no_electrons,
                                                                            air5_names,
                                                                            air5.species,
                                                                            air5_mech,
                                                                            air5.model,
                                                                            si_scales);
        expect(!vibration_only.electron && !vibration_only.electron_mode);
    };
    "electron-impact ionisation validates participants and charge at construction"_test = [&] {
        const auto source = *canary::detail::species_index(names, "O");
        const auto ion = *canary::detail::species_index(names, "O+");
        const auto carrier = *canary::electron_species_index(gas);
        const auto spectator = *canary::detail::species_index(names, "N2");
        const canary::EnergyExchangeInput<0, 0> heavy {};
        canary::ElectronExchangeInput<0, 0, 1> input {};
        input.electron_impact_ionisation[0] =
          canary::air11_kim_jo_exchange_electrons.electron_impact_ionisation[0];
        const auto build = [&](const auto& species, const auto& mechanism, const auto& labels) {
            return canary::make_multi_temperature_exchange(
              heavy, input, labels, species, mechanism, gas, si_scales);
        };
        const auto rejects = [&](const auto& species, const auto& mechanism, auto reason) {
            const auto message = error_message([&] { (void)build(species, mechanism, names); });
            expect(message.starts_with("electron-impact ionisation 0 (O): ")) << message;
            expect(message.find(reason) != std::string::npos) << message;
        };
        const auto reference = build(built.species, mech, names);
        auto renamed = names;
        renamed[ion.value()] = "oxygen ion";
        renamed[carrier.value()] = "carrier";
        const auto arbitrary_names = build(built.species, mech, renamed);
        expect(arbitrary_names.electron == reference.electron);
        expect(arbitrary_names.electron_impact_ionisation[0].species == source);
        expect(arbitrary_names.electron_impact_ionisation[0].ionisation_energy ==
               reference.electron_impact_ionisation[0].ionisation_energy);

        auto changed = mech;
        changed.reactions[38].reactants[spectator.value()] = 2;
        changed.reactions[38].products[spectator.value()] = 2;
        expect(error_message([&] { (void)build(built.species, changed, names); }) == "no exception");
        changed.reactions[38].products[spectator.value()] = 3;
        rejects(built.species, changed, "no other net species changes");
        changed.reactions[38].products[spectator.value()] = 1;
        rejects(built.species, changed, "no other net species changes");
        for (const int count : {0, 1, 3}) {
            changed = mech;
            changed.reactions[38].products[carrier.value()] = count;
            rejects(built.species, changed, "net electron production +1");
        }
        changed = mech;
        changed.reactions[38].reactants[carrier.value()] = 0;
        changed.reactions[38].products[carrier.value()] = 1;
        rejects(built.species, changed, "incident electron");
        changed = mech;
        changed.reactions[38].products[source.value()] = 1;
        rejects(built.species, changed, "must consume one");
        for (const int count : {0, 2}) {
            changed = mech;
            changed.reactions[38].products[ion.value()] = count;
            rejects(built.species, changed, "must produce one ion");
        }
        const double q = built.species.charge[ion.value()].value();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        for (const double charge : {q, -q, nan, inf}) {
            auto species = built.species;
            species.charge[source.value()] = canary::Charge<double, si> {charge};
            rejects(species, mech, "must be neutral");
        }
        for (const double charge : {0.0, q, nan, inf, -inf}) {
            auto species = built.species;
            species.charge[carrier.value()] = canary::Charge<double, si> {charge};
            rejects(species, mech, "electron charge must be finite and negative");
        }
        for (const double charge : {0.0, -q, 2.0 * q, nan, inf}) {
            auto species = built.species;
            species.charge[ion.value()] = canary::Charge<double, si> {charge};
            rejects(species, mech, "must have charge opposite to the electron charge");
        }
    };
    "electron-impact ionisation charge checks support float and double in SI and scaled units"_test =
      [] {
          const auto check = []<typename P, typename U>(const canary::UnitScales<P, U>& scales) {
              const auto built = canary::make_air11_three_temperature<P>(scales);
              const auto mechanism = canary::make_air11_kim_jo(
                scales,
                canary::convert<U>(canary::Pressure<P, si> {P(atmosphere)}, scales),
                built.model);
              const auto exchange = canary::make_air11_kim_jo_exchange(
                built.species, mechanism, built.model, scales);
              expect(exchange.electron == canary::electron_species_index(built.model));
              for (const auto& coupling : exchange.electron_impact_ionisation) {
                  expect(std::isfinite(coupling.ionisation_energy) && coupling.ionisation_energy > P(0));
              }
          };
          check(canary::UnitScales<float, si> {});
          check(canary::UnitScales<double, si> {});
          check(canary::UnitScales<float, scaled> {{
            .mass = 4.65e-26F, .density = 1e-3F, .temperature = 1000.0F, .length = 1.0F}});
          check(reference_scales);
      };
    "multi-temperature exchange rejects duplicate physical electron processes"_test = [&] {
        const auto build = [&](const auto& electrons) {
            return canary::make_multi_temperature_exchange(canary::air11_kim_jo_exchange,
                                                           electrons,
                                                           names,
                                                           built.species,
                                                           mech,
                                                           gas,
                                                           si_scales);
        };
        auto electrons = canary::air11_kim_jo_exchange_electrons;
        electrons.electron_translation[1] = electrons.electron_translation[0];
        expect(error_message([&] {
                   (void)build(electrons);
               }).starts_with("electron-translation pair 1 (N2): duplicates entry 0"));
        canary::ElectronExchangeInput<0, 2, 0> vibration {};
        vibration.vibration_electron.fill(
          canary::air11_kim_jo_exchange_electrons.vibration_electron[0]);
        expect(error_message([&] {
                   (void)build(vibration);
               }).starts_with("vibration-electron exchange 1 (N2): duplicates entry 0"));
        electrons = canary::air11_kim_jo_exchange_electrons;
        electrons.electron_impact_ionisation[1] = electrons.electron_impact_ionisation[0];
        expect(error_message([&] { (void)build(electrons); }).find("duplicates entry 0") !=
               std::string::npos);
        for (std::size_t i : {std::size_t(0), electron}) {
            auto species = built.species;
            species.mass[i] = canary::ParticleMass<double, si> {0.0};
            expect(error_message([&] {
                       (void)canary::make_air11_kim_jo_exchange(species, mech, gas, si_scales);
                   }).find("particle mass must be finite and positive") != std::string::npos);
        }
    };
}
