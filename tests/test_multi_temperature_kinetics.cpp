#include "dual.hpp"
#include "nasa9_blend_offsets.hpp"
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
#include <canary/data/air5_park_exchange.hpp>
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
#include <vector>

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
    canary::OptionalSpeciesIndex returned_electron;
};

canary::OptionalModeIndex species_vibration_mode(const SpeciesQueryModel& gas,
                                                 canary::SpeciesIndex index)
{ return canary::species_vibration_mode(gas.thermo, index); }

canary::OptionalModeIndex electron_energy_mode(const SpeciesQueryModel& gas)
{ return canary::electron_energy_mode(gas.thermo); }

canary::OptionalSpeciesIndex electron_species_index(const SpeciesQueryModel& gas)
{ return gas.returned_electron; }

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

static_assert(!canary::HasReactionTemperatureModes<SpeciesQueryModel>);
static_assert(canary::HasReactionTemperatureModes<Gas>);
struct TemperatureQueryModel {
    using parameter_type = double;
    using units = si;
    static constexpr std::size_t mode_count = 2;
    canary::OptionalModeIndex selected;
    friend canary::OptionalModeIndex species_vibration_mode(const TemperatureQueryModel& g,
                                                            canary::SpeciesIndex)
    { return g.selected; }
    friend canary::OptionalModeIndex species_electronic_mode(const TemperatureQueryModel& g,
                                                             canary::SpeciesIndex)
    { return g.selected; }
    friend canary::OptionalModeIndex electron_energy_mode(const TemperatureQueryModel& g)
    { return g.selected; }
};

template <typename G>
concept AcceptsScalarReactionTemperature =
  requires(const G& gas,
           const Mechanism& mechanism,
           const canary::NormalisedComposition<double, si, 11>& mix,
           Rho rho,
           Temp t) {
      canary::rates_of_progress(mechanism, gas, mix, rho, t, t);
      canary::production_rates(mechanism, gas, mix, rho, t, t);
  };
static_assert(!AcceptsScalarReactionTemperature<Gas>);

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
// electron-impact ionisation source. Its vibration-electron relaxation time
// is two quadratics in log10(T_e); for each state both are set to the
// constant log10(p_e tau / (1 atm s)) of the fit of air11_kim_jo_exchange at
// the state's T_e, so that its relaxation time at the state is that of the
// fit.
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
     .forward = {0.3311303189768377, 0.052627452969181254, 0.022448973269704985, 0.7025898268797712, 1.8043698160862143, 0.002649094430315865, 0.002105134210015285, 0.00448987674212909, 0.035130867313617975, 0.09022158439289515, 3.757021967031462, 1.1129291152244, 0.14030075107560014, 150.97993731443114, 39.65268737543988, 0.030056764352446936, 0.04451792784182268, 0.02806066324679883, 7.549292548637031, 1.9827023532283474, 2.281594200427389, 96.86190041881451, 13.839049306984272, 485.04950711670284, 132.0638911465566, 13.863990185844353, 42.5310028348855, 0.29336840211948706, 505.0594083052646, 26.676769158815784, 0.3577808907733399, 0.3716723115506606, 0.2994478545082911, 0.0008755228249258909, 15.54304492867512, 3.185563637673142, 44.658250254821446, 7.113714033138841, 1.4569416353480358e-05, 1.0065880537400348e-05},
     .backward = {0.00024092765228365137, 3.233845620010255e-05, 1.379441903154147e-05, 0.0005085202957490465, 0.0011087470697178016, 1.9274589646330224e-06, 1.2935604253122927e-06, 2.7589342477628097e-06, 2.5427010686419434e-05, 5.5439254430611795e-05, 7.034232353912862e-05, 7.572516736183032e-06, 2.6268360716401334e-06, 0.0016517966352282494, 0.0004338203919644147, 5.6274960891158e-07, 3.0290586258444556e-07, 5.253768197656986e-07, 8.259306668158507e-05, 2.1691763382449616e-05, 0.002916028667270711, 0.059519588112821036, 408.70612822019496, 82.53188397798903, 4806.771921048249, 8856.440429922985, 4900.148994199093, 31.059970408170955, 75.88039514271601, 160.92347260393558, 6.311094298776658, 0.3368476267356278, 3.024299066658695, 0.030883887029004877, 9.23703834804622, 55.90978505118614, 4.007480858105974, 22.51804272893099, 0.5108842557802363, 0.6613783687259123},
     .source = {-27.236319247012094, -24.962422780482374, 211.8704810281965, 366.01847376423103, 24.225885244661498, -8.414993808725754, -0.6073673595083084, -123.34500269987626, -280.0505497386352, -137.48810328826903, -0.010080745275893923},
     .total = {3053795132.878916, -2040379837.7646344},
     .vt = {914865860.7972299, 0.0},
     .cv = {-12993436.55042893, 0.0},
     .et = {0.0, 111542870.86748034},
     .ve = {2151922708.6321144, -2151922708.6321144},
   },
   {
     .name = "strongly ionised, T_e above T",
     .y = {0.10100951469793618, 0.010100951469793618, 0.4040380587917447, 0.3030285440938085, 0.010100951469793618, 0.08080761175834894, 0.0606057088187617, 0.010100951469793618, 0.005050475734896809, 0.015151427204690425, 5.8044904318440994e-06},
     .rho = 0.0001, .temp = 10000.0, .t_v = 9000.0, .t_e = 14000.0,
     .forward = {0.0015436797865539034, 0.00022687414659030279, 0.00024194082009703965, 0.032674095910462116, 0.058339066266077856, 0.00015437100177759523, 0.00011343901813825023, 0.0003629178653544314, 0.006535075141185032, 0.01166821334208387, 0.0059401130670267204, 0.0011010873358985613, 0.0005545631168057373, 1.8612427210022535, 0.36662127295767394, 0.0005940229397420298, 0.000550553106848065, 0.0008318598840605525, 0.37226312461911926, 0.07332676887740838, 0.07597563028424267, 13.036584635310339, 0.07657889474983046, 0.9394322571411249, 1.9794729291111564, 0.05076968986777312, 0.4286845401078314, 0.004245006643104379, 6.1417289925808, 0.8236250103810452, 0.004971439702997369, 0.0009621862755506518, 0.008628112706562255, 6.590530295421265e-05, 0.06410497790492585, 0.0771423898198302, 0.3501637348847107, 0.3012488975273424, 0.1983917438380295, 0.43808762302903503},
     .backward = {2.640476979964928e-07, 3.8025854426307823e-08, 4.055114495441056e-08, 5.585398998099e-06, 9.778076852479154e-06, 2.6405286907190576e-08, 1.9013253183837084e-08, 6.082782954372929e-08, 1.117123553964628e-06, 1.955682428471771e-06, 3.066548742884365e-08, 5.032364819835478e-09, 2.8628997621789197e-09, 9.006930278200952e-06, 1.7741545510286254e-06, 3.0666087977718997e-09, 2.5162255490772146e-09, 4.294428158079952e-09, 1.8014566132374844e-06, 3.548430773440072e-07, 9.725645628045439e-06, 0.0021850319968530997, 1.7194071326865983, 1.0895604578132572, 1260.6689412925487, 1340.1850791009358, 894.3528681530415, 0.4852479774358687, 4.583322114742326, 2.2714652875516053, 0.07300373293963283, 0.0019802748306447697, 0.03638305194972803, 0.0004423407642706464, 0.11833982048681882, 3.1974360726776134, 0.17090475991796653, 0.9868361769248365, 0.32753568988842613, 0.5269828380156827},
     .source = {-0.4127252429615826, -0.20150460538791815, 43.025418822089705, 63.188748780506266, 0.0429246737410809, -0.05104090367064907, -0.0042957236383597466, -25.01192211724084, -42.93628486038126, -37.637402581908546, -0.0019161709497631809},
     .total = {61477309.51722118, -85472019.28362702},
     .vt = {394932.73415702383, 0.0},
     .cv = {-629840.2088934951, 0.0},
     .et = {0.0, -23759802.291669354},
     .ve = {61712216.99195767, -61712216.99195767},
   },
   {
     .name = "equal temperatures",
     .y = {0.4999995398459613, 0.049999953984596134, 0.19999981593838453, 0.19999981593838453, 0.01999998159383845, 0.009999990796919225, 0.009999990796919225, 0.003999996318767691, 0.0019999981593838453, 0.003999996318767691, 9.203080774345904e-07},
     .rho = 0.01, .temp = 7000.0, .t_v = 7000.0, .t_e = 7000.0,
     .forward = {7.889144766709588, 1.303823119453573, 0.5561639164672647, 16.75948102778312, 44.70250695269393, 0.06311439413347045, 0.05215381892370655, 0.11123481699395397, 0.838006873579311, 2.235201989997017, 170.3436430983337, 63.773972308542156, 6.361246028735305, 7751.3888796605115, 2035.789691252314, 1.362775832686811, 2.5510026277330624, 1.2722724666398533, 387.5846244984217, 101.79297492013166, 81.08238260215231, 2399.713040154466, 434.8174756767775, 8456.607464613888, 1153.608333910327, 2.6686038443610265, 34.091619120509314, 0.2263003286065131, 18660.374091866433, 450.84465803486006, 2.6377043760877337, 3.296325761356547, 10.766466117761034, 0.07314527909117138, 178.9879879766242, 15.09907219039365, 946.7429829013997, 85.978685529195, 3.503848696933863, 4.2536308623646715},
     .backward = {0.5312153652619788, 0.08779289709072192, 0.037449283384732995, 1.1284992352204282, 3.0100421859676088, 0.004249806148110983, 0.003511776090898008, 0.007489993616118897, 0.05642717184237337, 0.15050727000992756, 0.05765297463681445, 0.021584364060306126, 0.002152970015684728, 2.6234652397399887, 0.6890150363163513, 0.0004612328296407379, 0.000863389364698752, 0.000430601875809138, 0.13117839984746307, 0.03445193313077409, 2.748821938124777, 161.5846174516668, 43554.32705916809, 16796.42166225658, 637445.63524915, 673755.8295727227, 735642.2493692847, 2911.9603050855117, 11266.765706503027, 12696.778185758802, 252.580141278417, 56.05000376509645, 317.1692081479675, 4.948720076173864, 816.4568369791563, 6898.964705850181, 546.4617599826104, 3357.5699566135618, 8820.179049524282, 10236.24377282351},
     .source = {-48.009907730826356, -2200.47642145432, 29232.077159504537, 33157.48764168009, 1029.6243767885517, -462.66811943265486, -148.44201981219757, -20481.27732131944, -21396.251740798583, -18680.930951901937, -1.1326540289057063},
     .total = {-2577241493.055196, 0.0},
     .vt = {0.0, 0.0},
     .cv = {-2577241493.055196, 0.0},
     .et = {0.0, 0.0},
     .ve = {0.0, 0.0},
   },
   {
     .name = "expansion, frozen modes",
     .y = {0.10100951469793618, 0.010100951469793618, 0.4040380587917447, 0.3030285440938085, 0.010100951469793618, 0.08080761175834894, 0.0606057088187617, 0.010100951469793618, 0.005050475734896809, 0.015151427204690425, 5.8044904318440994e-06},
     .rho = 0.0001, .temp = 3000.0, .t_v = 6000.0, .t_e = 5000.0,
     .forward = {1.6119429237806475e-09, 3.232054402411065e-10, 3.4466945858299213e-10, 3.445004886386648e-08, 8.31099703477131e-08, 1.6119744918590038e-10, 1.6160549074858384e-10, 5.170136404085771e-10, 6.8902796441500495e-09, 1.6622564037080822e-08, 3.178521125121944e-06, 3.7879241777514213e-06, 2.967436077547268e-07, 0.002674184052296462, 0.0005267516967631719, 3.1785833728618953e-07, 1.8939945602626306e-06, 4.451235498033094e-07, 0.0005348577592171515, 0.00010535395181172432, 4.221261147160491e-06, 1.857194898414229e-05, 0.0008187210757294595, 0.00040222184689005854, 0.00034759396528702354, 1.3372122359842202e-11, 1.0178449205757974e-08, 3.494199797631361e-11, 0.19523588761241895, 0.0008797001360215288, 1.017399397211982e-06, 7.226187532637879e-08, 0.00014419097305196522, 1.2204965412864353e-05, 1.945262032123211e-05, 6.467833448559583e-07, 0.001110678010913019, 7.612363316936983e-05, 0.00013770330574934254, 0.00017633476605538178},
     .backward = {1.1429905806240397e-06, 2.6198139040863747e-07, 2.7937953000916177e-07, 2.4529517944916796e-05, 6.736664324793535e-05, 1.1430129648013435e-07, 1.309929409987683e-07, 4.190769569764301e-07, 4.906095745888711e-06, 1.347379064950896e-05, 2.2944264099129995e-08, 6.093982306806669e-08, 2.142053942080338e-09, 2.9538091646649228e-05, 5.818313021744715e-06, 2.2944713436582476e-09, 3.0470433931124946e-08, 3.2131396587894473e-09, 5.907849721903216e-06, 1.1637025062949504e-06, 8.154301028596971e-06, 0.015053908170092325, 1.1151976264789052, 3.9601633233724916, 2647.7476924155444, 491.8640203934629, 1844.6109958379438, 0.4101435460550079, 5.3070553884692595, 2.138368083983717, 0.010457079400294212, 0.004774754006389753, 0.04644968934765906, 0.0014156579418250774, 0.0689608470267417, 3.1231955058862617, 0.22781651622650215, 1.8110828453228263, 165.26129034779856, 304.96261612254784},
     .source = {0.26099539880613537, -0.11706486432372398, 92.90002912533828, 60.85148154715887, -0.08713362143257451, -4.235582152831216, -2.640957593775767, -51.87322624857706, -15.77090099039122, -79.28464829092164, -0.0029921994320760955},
     .total = {-4853285.962910276, -37399717.59029706},
     .vt = {-50414.723649730324, 0.0},
     .cv = {-1626.6582417519821, 0.0},
     .et = {0.0, -42200962.17131586},
     .ve = {-4801244.581018793, 4801244.581018793},
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
{
    return mech.reactions()[canary::ReactionIndex {k}].products()[canary::SpeciesIndex {i}] -
           mech.reactions()[canary::ReactionIndex {k}].reactants()[canary::SpeciesIndex {i}];
}

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
    const auto r_records = mech.reactions();
    const auto& r = r_records[canary::ReactionIndex {k}];
    double out = 1.0;
    if (r.third_body) {
        out = 0.0;
        for (std::size_t i = 0; i < species_count; ++i) {
            out += r.efficiencies()[canary::SpeciesIndex {i}] * n[i];
        }
    }
    for (std::size_t i = 0; i < species_count; ++i) {
        const int nu = forward ? r.reactants()[canary::SpeciesIndex {i}]
                               : r.products()[canary::SpeciesIndex {i}];
        for (int c = 0; c < nu; ++c) { out *= n[i]; }
    }
    return out;
}

// ln(k_f) of reaction k at the rate-controlling temperature T_c, from the
// formula A T_c^n exp(-theta / T_c).
double forward_log_rate(const Mechanism& mech, std::size_t k, double t_c)
{
    const auto f_records = mech.reactions();
    const auto& f = f_records[canary::ReactionIndex {k}].forward;
    return std::log(f.A) + (f.n * std::log(t_c)) - (f.activation_temperature.value() / t_c);
}

// k_b of an equilibrium backward rate at the temperature t_eq from the
// formula k_f(t_eq) / K_c(t_eq), with the Gibbs energies of gas at t_eq.
double equilibrium_backward_rate(const Mechanism& mech, const Gas& gas, std::size_t k, double t_eq)
{
    const auto props = canary::species_properties(gas, Temp {t_eq});
    double ln_k_c = mech.reactions()[canary::ReactionIndex {k}].coefficient_change *
                    std::log(atmosphere / (boltzmann * t_eq));
    for (std::size_t i = 0; i < species_count; ++i) {
        ln_k_c -= coefficient_change(mech, k, i) * props.gibbs_energy[i].value() *
                  mech.mass()[canary::SpeciesIndex {i}].value() / (boltzmann * t_eq);
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
    for (const auto k : canary::reaction_indices(0, reaction_count)) {
        out.forward[k.value()] = q.forward()[k].value();
        out.backward[k.value()] = q.backward()[k].value();
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
    if (backward_model == canary::BackwardRateModel::equilibrium) {
        r.equilibrium_temperature = canary::TranslationTemperatureInput {};
    }
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

namespace
{
// Each field has its own accessor offset in addition to its layout mapping.
template <typename T> struct StateOffsetAccessor {
    using offset_policy = StateOffsetAccessor;
    using element_type = T;
    using reference = T&;
    using data_handle_type = T*;
    std::size_t offset_value;
    constexpr reference access(data_handle_type p, std::size_t i) const noexcept
    { return p[offset_value + i]; }
    constexpr data_handle_type offset(data_handle_type p, std::size_t i) const noexcept
    { return p + i; }
};

template <typename T, typename E, typename V> struct CustomModeState {
    using value_type = T;
    using units = si;
    canary::Density<T, si> rho;
    canary::SpecificEnergy<T, si> e;
    canary::Pressure<T, si> p;
    canary::Temperature<T, si> temp;
    E mode_e;
    V mode_temp;
};
using StateEnergyView =
  canary::IndexedView<canary::ModeIndex,
                      Kokkos::mdspan<const canary::SpecificEnergy<double, si>, Extents>>;
using StateTemperatureView = canary::IndexedView<canary::ModeIndex, RawModes>;
using BorrowedState =
  canary::MultiTemperatureStateView<double, si, StateEnergyView, StateTemperatureView>;
static_assert(canary::MultiTemperatureStateOf<BorrowedState, double, si, 2>);
static_assert(
  canary::MultiTemperatureStateOf<canary::MultiTemperatureState<double, si, 2>, double, si, 2>);
static_assert(
  !canary::MultiTemperatureStateOf<canary::MultiTemperatureState<double, si, 1>, double, si, 2>);
static_assert(!canary::MultiTemperatureStateOf<BorrowedState, float, si, 2>);
static_assert(!canary::MultiTemperatureStateOf<BorrowedState, double, scaled, 2>);
static_assert(!canary::MultiTemperatureStateOf<canary::ThermoState<double, si>, double, si, 2>);
using WrongDomain = canary::IndexedView<canary::SpeciesIndex, RawModes>;
static_assert(
  !canary::
    MultiTemperatureStateOf<CustomModeState<double, StateEnergyView, WrongDomain>, double, si, 2>);
static_assert(!canary::MultiTemperatureStateOf<
              CustomModeState<double, StateTemperatureView, StateTemperatureView>,
              double,
              si,
              2>);
using FixedTemperatures = canary::IndexedView<
  canary::ModeIndex,
  Kokkos::mdspan<const canary::ModeTemperature<double, si>, Kokkos::extents<std::size_t, 2>>>;
static_assert(
  canary::MultiTemperatureStateOf<CustomModeState<double, StateEnergyView, FixedTemperatures>,
                                  double,
                                  si,
                                  2>);
static_assert(
  !canary::MultiTemperatureStateOf<CustomModeState<double, StateEnergyView, FixedTemperatures>,
                                   double,
                                   si,
                                   3>);
static_assert(std::is_trivially_copyable_v<BorrowedState>);

// Air11 model whose species_mode_energy records each call as (mode, temperature).
struct CountedModeGas {
    using parameter_type = double;
    using units = si;
    static constexpr std::size_t mode_count = Gas::mode_count;
    Gas model;
    std::vector<std::pair<std::size_t, double>>* calls;
};
[[maybe_unused]] canary::OptionalModeIndex species_vibration_mode(const CountedModeGas& gas,
                                                                  canary::SpeciesIndex i)
{ return canary::species_vibration_mode(gas.model, i); }
[[maybe_unused]] canary::OptionalModeIndex electron_energy_mode(const CountedModeGas& gas)
{ return canary::electron_energy_mode(gas.model); }
template <canary::Scalar T>
auto species_properties(const CountedModeGas& gas, canary::Temperature<T, si> t)
{ return canary::species_properties(gas.model, t); }
template <canary::Scalar T>
auto species_mode_energy(const CountedModeGas& gas,
                         canary::ModeIndex m,
                         canary::Temperature<T, si> t)
{
    gas.calls->emplace_back(m.value(), canary::ScalarTraits<T>::value(t.value()));
    return canary::species_mode_energy(gas.model, m, t);
}
static_assert(canary::ModeRelaxingModel<CountedModeGas, double, double, species_count>);
static_assert(canary::ReactingModel<CountedModeGas, double, double, species_count>);

// Electron-translation pairs of the air11 Kim-Jo exchange with cross-section K, in list order.
template <canary::ElectronCrossSection K> constexpr auto electron_pairs()
{
    constexpr auto& all = canary::air11_kim_jo_exchange_electrons.electron_translation;
    constexpr std::size_t count = static_cast<std::size_t>(
      std::ranges::count_if(all, [](const auto& pair) { return pair.cross_section == K; }));
    canary::ElectronExchangeInput<count, 0, 0> out {};
    std::size_t k = 0;
    for (const auto& pair : all) {
        if (pair.cross_section == K) { out.electron_translation[k++] = pair; }
    }
    return out;
}

// The read-only contract of EnergyExchange also holds for the multi-temperature exchange.
template <typename X>
concept BorrowsSpeciesModes = requires { std::declval<X>().vibration_mode(); };
template <typename X>
concept BorrowsElectronPairs = requires { std::declval<X>().electron_translation(); };
template <typename X>
concept ReadsElectron = requires {
    std::declval<X>().electron();
    std::declval<X>().electron_mode();
    std::declval<X>().coulomb_cross_section();
    std::declval<X>().coulomb_logarithm_factor();
    std::declval<X>().minimum_electron_density();
};
template <typename X, typename Index>
concept SpeciesModesIndexedBy = requires(const X& x, Index i) { x.vibration_mode()[i]; };
static_assert(!std::is_default_constructible_v<Exchange>);
static_assert(std::is_trivially_copyable_v<Exchange>);
static_assert(BorrowsSpeciesModes<const Exchange&> && BorrowsElectronPairs<const Exchange&>);
static_assert(!BorrowsSpeciesModes<Exchange&&> && !BorrowsSpeciesModes<const Exchange&&>);
static_assert(!BorrowsElectronPairs<Exchange&&> && !BorrowsElectronPairs<const Exchange&&>);
static_assert(ReadsElectron<Exchange&&>);
static_assert(SpeciesModesIndexedBy<Exchange, canary::SpeciesIndex>);
static_assert(!SpeciesModesIndexedBy<Exchange, std::size_t>);
static_assert(!SpeciesModesIndexedBy<Exchange, canary::ModeIndex>);
static_assert(
  std::is_const_v<decltype(std::declval<const Exchange&>().vibration_mode())::element_type>);
static_assert(noexcept(std::declval<const Exchange&>().vibration_mode()) &&
              noexcept(std::declval<const Exchange&>().electron_impact_ionisation()) &&
              noexcept(std::declval<const Exchange&>().electron()) &&
              noexcept(std::declval<const Exchange&>().minimum_electron_density()));
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
    const auto& exchange_names = canary::detail::air11_kim_jo_exchange_species_names;
    const canary::EnergyExchangeInput<0, 0> no_vibration {};
    const canary::ElectronExchangeInput<0, 0, 0> no_electrons {};
    const auto build_exchange = [&](const auto& heavy, const auto& electrons) {
        return canary::make_multi_temperature_exchange(heavy,
                                                       electrons,
                                                       exchange_names,
                                                       built.species,
                                                       mech,
                                                       gas,
                                                       si_scales);
    };

    "built exchanges prepare their mode and Coulomb requirements from their records"_test = [&] {
        const auto check = [&](const auto& x, bool translation_0, bool electron_0, bool coulomb) {
            const auto modes = x.vibration_mode();
            for (const canary::ModeIndex m : canary::mode_indices(0, mode_count)) {
                bool translation = false;
                bool electron_mode = false;
                for (const auto& pair : x.vibration_translation()) {
                    translation =
                      translation || modes[pair.species] == canary::OptionalModeIndex {m};
                }
                for (const auto& c : x.vibration_electron()) {
                    electron_mode =
                      electron_mode || modes[c.species] == canary::OptionalModeIndex {m};
                }
                expect(canary::detail::uses_translation_mode(x, m) == translation) << m.value();
                expect(canary::detail::uses_electron_mode(x, m) == electron_mode) << m.value();
            }
            const canary::ModeIndex vibration {0};
            const canary::ModeIndex free_electron {1};
            expect(canary::detail::uses_translation_mode(x, vibration) == translation_0);
            expect(canary::detail::uses_electron_mode(x, vibration) == electron_0);
            expect(!canary::detail::uses_translation_mode(x, free_electron));
            expect(!canary::detail::uses_electron_mode(x, free_electron));
            expect(canary::detail::uses_coulomb(x) == coulomb);
        };
        check(exchange, true, true, true);
        check(build_exchange(no_vibration, no_electrons), false, false, false);
        check(build_exchange(canary::air11_kim_jo_exchange, no_electrons), true, false, false);
        check(build_exchange(no_vibration,
                             electron_pairs<canary::ElectronCrossSection::gnoffo_neutral>()),
              false,
              false,
              false);
        check(build_exchange(no_vibration, electron_pairs<canary::ElectronCrossSection::coulomb>()),
              false,
              false,
              true);
        auto copy = exchange;
        copy = canary::make_air11_kim_jo_exchange(built.species, mech, gas, si_scales);
        check(copy, true, true, true);
        expect(&copy.vibration_translation() != &exchange.vibration_translation());
        expect(copy.mass().base().data_handle() != exchange.mass().base().data_handle());
        for (const auto i : canary::species_indices(0, species_count)) {
            expect(copy.mass()[i] == exchange.mass()[i] &&
                   copy.vibration_mode()[i] == exchange.vibration_mode()[i])
              << i.value();
        }
    };

    "exchange sources evaluate mode energies only for configured processes"_test = [&] {
        std::vector<std::pair<std::size_t, double>> calls;
        const CountedModeGas counted {gas, &calls};
        Fractions y {};
        y[0] = 0.7;
        y[1] = 0.2;
        y[5] = 0.1 - 1e-6;
        y[electron] = 1e-6;
        const auto mix =
          canary::mixture(gas,
                          canary::as_kind<canary::kind::mass_fraction>(
                            Kokkos::mdspan<const double, Extents>(y.data(), y.size())));
        canary::MultiTemperatureState<double, si, mode_count> state {};
        state.rho = Rho {0.01};
        state.temp = Temp {8000};
        state.mode_temp = {canary::ModeTemperature<double, si> {5000},
                           canary::ModeTemperature<double, si> {11000}};
        using Calls = std::vector<std::pair<std::size_t, double>>;
        const auto zero = [](const auto& modes) {
            return std::ranges::all_of(modes, [](const auto& q) { return q.value() == 0.0; });
        };

        const auto empty = build_exchange(no_vibration, no_electrons);
        calls.clear();
        expect(zero(canary::vibration_translation_source(empty, counted, mix, state)));
        expect(zero(canary::vibration_electron_source(empty, counted, mix, state)));
        expect(zero(canary::electron_translation_source(empty, mix, state)));
        const auto empty_total =
          canary::multi_temperature_sources(mech, empty, counted, mix, state);
        expect(calls.empty());
        expect(zero(empty_total.modes));
        const auto production =
          canary::production_rates(mech, gas, mix, state.rho, state.temp, state.mode_temp);
        for (const auto i : canary::species_indices(0, species_count)) {
            expect(empty_total.species[i.value()] == production[i.value()]) << i.value();
        }

        calls.clear();
        (void)canary::vibration_translation_source(exchange, counted, mix, state);
        expect(calls == Calls {{0, 5000.0}, {0, 8000.0}});
        calls.clear();
        (void)canary::vibration_electron_source(exchange, counted, mix, state);
        expect(calls == Calls {{0, 5000.0}, {0, 11000.0}});
        calls.clear();
        (void)canary::multi_temperature_sources(mech, exchange, counted, mix, state);
        expect(calls == Calls {{0, 5000.0}, {0, 8000.0}, {0, 11000.0}});
        calls.clear();
        (void)canary::multi_temperature_sources(
          mech,
          build_exchange(canary::air11_kim_jo_exchange, no_electrons),
          counted,
          mix,
          state);
        expect(calls == Calls {{0, 5000.0}, {0, 8000.0}});
    };

    "stateful accessor states match contiguous inputs across model sizes"_test = [] {
        const auto run = []<typename T, typename Layout, bool Const, bool FiveSpecies>() {
            using P = canary::ScalarTraits<T>::value_type;
            using Element = std::conditional_t<Const, const T, T>;
            const canary::UnitScales<P, si> scales {};
            const auto built = [&] {
                if constexpr (FiveSpecies) {
                    return canary::make_air5_four_temperature<P>(scales);
                } else {
                    return canary::make_air11_three_temperature<P>(scales);
                }
            }();
            const auto& model = built.model;
            constexpr auto n = std::tuple_size_v<decltype(built.species.mass)>;
            constexpr auto m_count = std::remove_cvref_t<decltype(model)>::mode_count;
            const auto mechanism = [&] {
                if constexpr (FiveSpecies) {
                    return canary::make_air5_park(scales,
                                                  canary::Pressure<P, si> {P(atmosphere)},
                                                  model);
                } else {
                    return canary::make_air11_kim_jo(scales,
                                                     canary::Pressure<P, si> {P(atmosphere)},
                                                     model);
                }
            }();
            const auto exchanges = [&] {
                if constexpr (FiveSpecies) {
                    return canary::make_multi_temperature_exchange(
                      canary::air5_park_exchange,
                      canary::ElectronExchangeInput<0, 0, 0> {},
                      canary::detail::air5_species_names,
                      built.species,
                      mechanism,
                      model,
                      scales);
                } else {
                    return canary::make_air11_kim_jo_exchange(built.species,
                                                              mechanism,
                                                              model,
                                                              scales);
                }
            }();
            const auto seed = [](P value, P derivative) {
                if constexpr (requires { T(value, derivative); }) {
                    return T(value, derivative);
                } else {
                    return T(value);
                }
            };
            const auto equal = [](const T& a, const T& b) {
                expect(canary::ScalarTraits<T>::value(a) == canary::ScalarTraits<T>::value(b));
                if constexpr (requires { a.dual(); }) { expect(a.dual() == b.dual()); }
            };
            std::array<T, n> y {};
            for (const auto i : canary::species_indices(0, n)) {
                const auto value = [&] {
                    if constexpr (FiveSpecies) {
                        return std::array<P, 5> {P(0.65),
                                                 P(0.2),
                                                 P(0.05),
                                                 P(0.05),
                                                 P(0.05)}[i.value()];
                    } else {
                        return P(reference_states[0].y[i.value()]);
                    }
                }();
                const P derivative = i.value() == 0 ? P(1) : (i.value() == 1 ? P(-1) : P(0));
                y[i.value()] = seed(value, derivative);
            }
            std::array<T, m_count> temperatures {};
            for (const auto m : canary::mode_indices(0, m_count)) {
                temperatures[m.value()] =
                  seed(P(5000 + 1000 * m.value()), m.value() == 0 ? P(1) : P(0));
            }
            const auto reference_mix =
              canary::mixture(model,
                              canary::as_kind<canary::kind::mass_fraction>(
                                Kokkos::mdspan<const T, Extents>(y.data(), n)));
            const auto contiguous_modes = canary::as_kind<canary::kind::mode_temperature, si>(
              Kokkos::mdspan<const T, Extents>(temperatures.data(), m_count));
            const canary::Density<T, si> rho {T(P(0.001))};
            const canary::Temperature<T, si> temp {T(P(10000))};
            const auto reference =
              canary::thermo_state(model, reference_mix, rho, temp, contiguous_modes);
            expect(reference.status == canary::Status::Valid);

            std::array<T, 3 * n + 8> y_storage;
            std::array<T, 3 * m_count + 8> t_storage;
            std::array<T, 3 * m_count + 8> e_storage;
            y_storage.fill(T(P(-7)));
            t_storage.fill(T(P(-11)));
            e_storage.fill(T(P(-13)));
            using Grid = Kokkos::mdspan<T, Kokkos::dextents<std::size_t, 2>, Layout>;
            using OtherLayout = std::conditional_t<std::same_as<Layout, Kokkos::layout_left>,
                                                   Kokkos::layout_right,
                                                   Kokkos::layout_left>;
            const Grid y_grid(y_storage.data() + 2, 3, n);
            const Grid t_grid(t_storage.data() + 3, 3, m_count);
            const Kokkos::mdspan<T, Kokkos::dextents<std::size_t, 2>, OtherLayout> e_grid(
              e_storage.data() + 4,
              3,
              m_count);
            for (const auto i : canary::species_indices(0, n)) {
                y_storage[2 + y_grid.mapping()(1, i.value())] = y[i.value()];
            }
            for (const auto m : canary::mode_indices(0, m_count)) {
                t_storage[3 + t_grid.mapping()(1, m.value())] = temperatures[m.value()];
                e_storage[4 + e_grid.mapping()(1, m.value())] =
                  reference.state.mode_e[m.value()].value();
            }
            const auto slice = [](auto& storage, const auto& grid, std::size_t offset) {
                using Raw = Kokkos::
                  mdspan<Element, Extents, Kokkos::layout_stride, StateOffsetAccessor<Element>>;
                return Raw(storage.data() + 1 + grid.mapping()(1, 0),
                           typename Raw::mapping_type(Extents(grid.extent(1)),
                                                      std::array<std::size_t, 1> {grid.stride(1)}),
                           StateOffsetAccessor<Element> {offset});
            };
            const auto y_before = y_storage;
            const auto t_before = t_storage;
            const auto e_before = e_storage;
            const auto mix = canary::mixture(
              model,
              canary::as_kind<canary::kind::mass_fraction>(slice(y_storage, y_grid, 1)));
            const auto mode_temp = canary::as_indexed<canary::ModeIndex>(
              canary::as_kind<canary::kind::mode_temperature, si>(slice(t_storage, t_grid, 2)));
            const auto mode_e = canary::as_indexed<canary::ModeIndex>(
              canary::as_kind<canary::kind::specific_energy, si>(slice(e_storage, e_grid, 3)));
            const auto closed = canary::thermo_state(model, mix, rho, temp, mode_temp);
            expect(closed.status == canary::Status::Valid);
            equal(closed.state.e.value(), reference.state.e.value());
            equal(closed.state.p.value(), reference.state.p.value());
            const canary::MultiTemperatureStateView<T, si, decltype(mode_e), decltype(mode_temp)>
              state {closed.state.rho,
                     closed.state.e,
                     closed.state.p,
                     closed.state.temp,
                     mode_e,
                     mode_temp};
            const auto compare_modes = [&](const auto& a, const auto& b) {
                for (const auto m : canary::mode_indices(0, m_count)) {
                    equal(a[m.value()].value(), b[m.value()].value());
                }
            };
            compare_modes(closed.state.mode_e, reference.state.mode_e);
            compare_modes(closed.state.mode_temp, reference.state.mode_temp);
            for (const auto m : canary::mode_indices(0, m_count)) {
                equal(state.mode_e[m].value(), reference.state.mode_e[m.value()].value());
            }
            const auto rates =
              canary::rates_of_progress(mechanism, model, mix, rho, temp, mode_temp);
            const auto reference_rates = canary::rates_of_progress(mechanism,
                                                                   model,
                                                                   reference_mix,
                                                                   rho,
                                                                   temp,
                                                                   reference.state.mode_temp);
            for (const auto r : canary::reaction_indices(0, rates.reaction_count)) {
                equal(rates.forward()[r].value(), reference_rates.forward()[r].value());
                equal(rates.backward()[r].value(), reference_rates.backward()[r].value());
            }
            compare_modes(canary::vibration_translation_source(exchanges, model, mix, state),
                          canary::vibration_translation_source(exchanges,
                                                               model,
                                                               reference_mix,
                                                               reference.state));
            compare_modes(canary::chemistry_vibration_source(exchanges, mechanism, rates, state),
                          canary::chemistry_vibration_source(exchanges,
                                                             mechanism,
                                                             reference_rates,
                                                             reference.state));
            compare_modes(
              canary::electron_translation_source(exchanges, mix, state),
              canary::electron_translation_source(exchanges, reference_mix, reference.state));
            compare_modes(
              canary::vibration_electron_source(exchanges, model, mix, state),
              canary::vibration_electron_source(exchanges, model, reference_mix, reference.state));
            const auto sources =
              canary::multi_temperature_sources(mechanism, exchanges, model, mix, state);
            const auto expected = canary::multi_temperature_sources(mechanism,
                                                                    exchanges,
                                                                    model,
                                                                    reference_mix,
                                                                    reference.state);
            compare_modes(sources.modes, expected.modes);
            for (const auto i : canary::species_indices(0, n)) {
                equal(sources.species[i.value()].value(), expected.species[i.value()].value());
            }
            const auto unchanged = [&](const auto& after, const auto& before) {
                for (std::size_t i = 0; i < after.size(); ++i) { equal(after[i], before[i]); }
            };
            unchanged(y_storage, y_before);
            unchanged(t_storage, t_before);
            unchanged(e_storage, e_before);
        };
        run.template operator()<double, Kokkos::layout_left, true, true>();
        run.template operator()<double, Kokkos::layout_right, false, true>();
        run.template operator()<D, Kokkos::layout_left, false, true>();
        run.template operator()<D, Kokkos::layout_right, true, true>();
        run.template operator()<double, Kokkos::layout_left, false, false>();
        run.template operator()<double, Kokkos::layout_right, true, false>();
        run.template operator()<D, Kokkos::layout_left, true, false>();
        run.template operator()<D, Kokkos::layout_right, false, false>();
    };

    "borrowed states preserve all multi-temperature sources"_test = [&] {
        const auto run = []<typename T, typename Layout, bool Const>() {
            using P = canary::ScalarTraits<T>::value_type;
            using Element = std::conditional_t<Const, const T, T>;
            const canary::UnitScales<P, si> scales {};
            const auto built_gas = canary::make_air11_three_temperature<P>(scales);
            const auto& model = built_gas.model;
            const auto mechanism =
              canary::make_air11_kim_jo(scales, canary::Pressure<P, si> {P(atmosphere)}, model);
            const auto exchanges =
              canary::make_air11_kim_jo_exchange(built_gas.species, mechanism, model, scales);
            const auto seed = [](P value) {
                if constexpr (requires { T(value, P(1)); }) {
                    return T(value, P(1));
                } else {
                    return T(value);
                }
            };
            const auto equal = [](const auto& a, const auto& b) {
                expect(canary::ScalarTraits<T>::value(a) == canary::ScalarTraits<T>::value(b));
                if constexpr (requires { a.dual(); }) { expect(a.dual() == b.dual()); }
            };
            // Select the middle cell of each multidimensional field. Scalar
            // quantities are loaded once; mode descriptors borrow the selected cell.
            std::array<T, 35> composition {};
            using Grid = Kokkos::mdspan<T, Kokkos::dextents<std::size_t, 2>, Layout>;
            const Grid species_grid(composition.data() + 1, 3, species_count);
            for (const auto i : canary::species_indices(0, species_count)) {
                composition[1 + species_grid.mapping()(1, i.value())] =
                  T(P(reference_states[0].y[i.value()]));
            }
            using Raw = Kokkos::mdspan<Element, Extents, Kokkos::layout_stride>;
            const auto slice = [](auto& storage, const auto& grid) {
                return Raw(storage.data() + 1 + grid.mapping()(1, 0),
                           typename Raw::mapping_type(Extents(grid.extent(1)),
                                                      std::array<std::size_t, 1> {grid.stride(1)}));
            };
            const auto mix = canary::mixture(
              model,
              canary::as_kind<canary::kind::mass_fraction>(slice(composition, species_grid)));
            std::array<T, 8> temperatures {};
            std::array<T, 8> energies {};
            const Grid temp_grid(temperatures.data() + 1, 3, 2);
            using OtherLayout = std::conditional_t<std::same_as<Layout, Kokkos::layout_left>,
                                                   Kokkos::layout_right,
                                                   Kokkos::layout_left>;
            const Kokkos::mdspan<T, Kokkos::dextents<std::size_t, 2>, OtherLayout> energy_grid(
              energies.data() + 1,
              3,
              2);
            temperatures[1 + temp_grid.mapping()(1, 0)] = seed(P(5000));
            temperatures[1 + temp_grid.mapping()(1, 1)] = T(P(8000));
            const auto mode_temp = canary::as_indexed<canary::ModeIndex>(
              canary::as_kind<canary::kind::mode_temperature, si>(slice(temperatures, temp_grid)));
            const auto closed = canary::thermo_state(model,
                                                     mix,
                                                     canary::Density<T, si> {T(P(0.001))},
                                                     canary::Temperature<T, si> {T(P(10000))},
                                                     mode_temp);
            expect(closed.status == canary::Status::Valid);
            const auto& owned = closed.state;
            for (const auto m : canary::mode_indices(0, 2)) {
                energies[1 + energy_grid.mapping()(1, m.value())] = owned.mode_e[m.value()].value();
            }
            const auto mode_e = canary::as_indexed<canary::ModeIndex>(
              canary::as_kind<canary::kind::specific_energy, si>(slice(energies, energy_grid)));
            const canary::MultiTemperatureStateView<T, si, decltype(mode_e), decltype(mode_temp)>
              view {owned.rho, owned.e, owned.p, owned.temp, mode_e, mode_temp};
            const auto rates =
              canary::rates_of_progress(mechanism, model, mix, view.rho, view.temp, view.mode_temp);
            const auto compare_modes = [&](const auto& a, const auto& b) {
                for (const auto m : canary::mode_indices(0, 2)) {
                    equal(a[m.value()].value(), b[m.value()].value());
                }
            };
            compare_modes(canary::vibration_translation_source(exchanges, model, mix, view),
                          canary::vibration_translation_source(exchanges, model, mix, owned));
            compare_modes(canary::chemistry_vibration_source(exchanges, mechanism, rates, view),
                          canary::chemistry_vibration_source(exchanges, mechanism, rates, owned));
            compare_modes(canary::electron_translation_source(exchanges, mix, view),
                          canary::electron_translation_source(exchanges, mix, owned));
            compare_modes(canary::vibration_electron_source(exchanges, model, mix, view),
                          canary::vibration_electron_source(exchanges, model, mix, owned));
            const auto expected =
              canary::multi_temperature_sources(mechanism, exchanges, model, mix, owned);
            const auto actual =
              canary::multi_temperature_sources(mechanism, exchanges, model, mix, view);
            compare_modes(actual.modes, expected.modes);
            for (const auto i : canary::species_indices(0, species_count)) {
                equal(actual.species[i.value()].value(), expected.species[i.value()].value());
            }
            const CustomModeState<T, decltype(mode_e), decltype(mode_temp)> custom {view.rho,
                                                                                    view.e,
                                                                                    view.p,
                                                                                    view.temp,
                                                                                    view.mode_e,
                                                                                    view.mode_temp};
            const auto custom_sources =
              canary::multi_temperature_sources(mechanism, exchanges, model, mix, custom);
            compare_modes(custom_sources.modes, expected.modes);
            const auto copy = view;
            temperatures[1 + temp_grid.mapping()(1, 0)] = T(P(6000));
            equal(copy.mode_temp[canary::ModeIndex {0}].value(), T(P(6000)));
            equal(copy.temp.value(), owned.temp.value());
        };
        run.template operator()<double, Kokkos::layout_left, true>();
        run.template operator()<double, Kokkos::layout_right, false>();
        run.template operator()<float, Kokkos::layout_left, false>();
        run.template operator()<float, Kokkos::layout_right, true>();
        run.template operator()<D, Kokkos::layout_left, false>();
        run.template operator()<D, Kokkos::layout_right, true>();
        run.template operator()<canary::test::Dual<float>, Kokkos::layout_left, true>();
        run.template operator()<canary::test::Dual<float>, Kokkos::layout_right, false>();
    };

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
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                equal(typed.forward()[k].value(), expected.forward()[k].value());
                equal(typed.backward()[k].value(), expected.backward()[k].value());
                equal(untyped.forward()[k].value(), expected.forward()[k].value());
                equal(untyped.backward()[k].value(), expected.backward()[k].value());
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
        const SpeciesQueryModel custom {gas, canary::electron_species_index(gas)};
        const auto result = canary::make_multi_temperature_exchange(input,
                                                                    electrons,
                                                                    renamed,
                                                                    built.species,
                                                                    mech,
                                                                    custom,
                                                                    si_scales);
        expect(result.electron() == canary::OptionalSpeciesIndex {canary::SpeciesIndex {10}});
        expect(result.electron_translation()[0].partner == canary::SpeciesIndex {2});
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
            invalid.returned_electron = canary::OptionalSpeciesIndex {canary::SpeciesIndex {index}};
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
        for (const auto k : canary::reaction_indices(0, reaction_count)) {
            const auto r_records = mech.reactions();
            const auto& r = r_records[k];
            const bool forward_electron =
              r.forward.translational_exponent != 1.0 &&
              r.forward.other_mode == canary::OptionalModeIndex {canary::ModeIndex {1}};
            expect(forward_electron == contains(forward_at_electron, k.value())) << k.value();
            const bool backward_electron = r.equilibrium_at_other_temperature;
            expect(backward_electron == contains(backward_at_electron, k.value())) << k.value();
            if (backward_electron) {
                expect(r.backward.other_mode == canary::OptionalModeIndex {canary::ModeIndex {1}})
                  << k.value();
            }
            // Dissociations of N2, O2 and NO: their vibration mode.
            if (k.value() <= 20) {
                expect(r.forward.other_mode == canary::OptionalModeIndex {canary::ModeIndex {0}})
                  << k.value();
            }
        }
        expect(exchange.electron() == canary::OptionalSpeciesIndex {canary::SpeciesIndex {10}} &&
               exchange.electron_mode() == canary::OptionalModeIndex {canary::ModeIndex {1}});
        for (const std::size_t i : std::array<std::size_t, 6> {0, 1, 4, 7, 8, 9}) {
            expect(exchange.vibration_mode()[canary::SpeciesIndex {i}] ==
                   canary::OptionalModeIndex {canary::ModeIndex {0}})
              << i;
        }
        for (const std::size_t i : std::array<std::size_t, 5> {2, 3, 5, 6, 10}) {
            expect(!exchange.vibration_mode()[canary::SpeciesIndex {i}]) << i;
        }
        const auto& ve = exchange.vibration_electron()[0];
        expect(ve.species == canary::SpeciesIndex {0});
        // In SI the fit is that of the input with ln(101325) added to the
        // constant, for p_e in Pa.
        expect(ve.relaxation_time.reference_temperature.value() == 10000.0_d);
        expect(ve.relaxation_time.exponent == 4.475475);
        expect(ve.relaxation_time.inverse_coefficients ==
               std::array {3.192608, -0.08318516, 0.002015418});
        expect(relative(ve.relaxation_time.constant, -24.91199 + std::log(atmosphere), 1.0e-15));
        const auto& ce = exchange.electron_impact_ionisation()[1];
        expect(ce.reaction == canary::ReactionIndex {39} && ce.species == canary::SpeciesIndex {2});
        expect(relative(ce.ionisation_energy, 1402330.360173364 / avogadro, 1.0e-15));
        // (2 pi / 9) C^2 and 9 / (4 pi C^3) with C = e^2 / (4 pi eps_0 k_B),
        // the exact e and k_B and CODATA 2018 eps_0; a few ulp.
        const double e = canary::constants::elementary_charge<double>;
        const double c = e * e / (4.0 * std::numbers::pi * 8.8541878128e-12 * boltzmann);
        expect(relative(exchange.coulomb_cross_section(),
                        2.0 * std::numbers::pi / 9.0 * c * c,
                        1.0e-15));
        expect(relative(exchange.coulomb_logarithm_factor(),
                        9.0 / (4.0 * std::numbers::pi * c * c * c),
                        1.0e-15));
        expect(exchange.minimum_electron_density() == 1.0_d);
        expect(exchange.electron_translation()[0].cross_section ==
               canary::ElectronCrossSection::gnoffo_neutral);
        expect(exchange.electron_translation()[5].cross_section ==
               canary::ElectronCrossSection::coulomb);
        // kim_2021 of N2-N2 with omega = 0.55.
        expect(relative(exchange.vibration_translation()[0].correction_exponent, -1.05, 1.0e-15));
    };

    "multi-temperature kinetics: rates and mode sources equal the reference"_test = [&] {
        const double r_u =
          canary::constants::universal_gas_constant<double> / reference_gas_constant;
        for (const ReferenceState& s : reference_states) {
            const Values got = evaluate(mech, exchange, gas, s);
            const Fractions n = number_densities(built.species, s.y, s.rho);
            const double t_l = std::clamp(s.temp, 300.0, 30000.0);
            // The reference joins NASA9 segments by blending their
            // coefficients; canary's Gibbs energies differ by the blend
            // offsets at T_l, which scale the equilibrium backward rates, and
            // its vibrational energy differences e_v(T) - e_v(T_v) and
            // e_v(T_e) - e_v(T_v) exceed the reference's by
            // R_i (dh_i(T) - dh_i(T_v)) and R_i (dh_i(T_e) - dh_i(T_v)), which
            // the same relaxation times turn into source differences (see
            // nasa9_blend_offsets.hpp).
            const auto off = canary::test::blend_offsets(canary::air11_species, t_l);
            std::array<canary::PowerDensity<double, si>, mode_count> dq_vt {};
            std::array<canary::PowerDensity<double, si>, mode_count> dq_ve {};
            {
                const auto mix = canary::mixture(gas, fractions(s.y));
                const auto state = cell_state(gas, s.y, s.rho, s.temp, s.t_v, s.t_e);
                const auto d = canary::detail::partial_densities(exchange.mass(), mix, state.rho);
                const auto off_t = canary::test::blend_offsets(canary::air11_species, s.temp);
                const auto off_v = canary::test::blend_offsets(canary::air11_species, s.t_v);
                const auto off_e = canary::test::blend_offsets(canary::air11_species, s.t_e);
                std::array<double, species_count> dvt {};
                std::array<double, species_count> dve {};
                for (std::size_t i = 0; i < species_count; ++i) {
                    const double r_i = gas.thermally_perfect.gas_constant[i].value();
                    dvt[i] = r_i * static_cast<double>(off_t[i].h_over_r - off_v[i].h_over_r);
                    dve[i] = r_i * static_cast<double>(off_e[i].h_over_r - off_v[i].h_over_r);
                }
                dq_vt = canary::detail::vibration_translation_terms(exchange, d, dvt, state);
                dq_ve = canary::detail::vibration_electron_terms(exchange, d, dve, state);
            }
            // The reference rates in particles per volume and time, with
            // the reactions that differ in the reference replaced by the
            // formulas at canary's temperatures.
            PerReaction forward {};
            PerReaction backward {};
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                const int dnu = mech.reactions()[k].coefficient_change;
                forward[k.value()] = avogadro * s.forward[k.value()];
                double offset = 0.0;
                std::array<int, species_count> dnu_i {};
                for (std::size_t i = 0; i < species_count; ++i) {
                    dnu_i[i] = coefficient_change(mech, k.value(), i);
                    offset += dnu_i[i] * enthalpy_offset[i];
                }
                backward[k.value()] = avogadro * s.backward[k.value()] * std::pow(r_u, dnu) *
                                      std::exp(offset / t_l) *
                                      canary::test::blend_rate_factor(off, dnu_i, t_l);
                // The reference evaluates these forward rates at T_v; the
                // formula at T_c = T^s T_e^(1 - s).
                if (contains(forward_at_electron, k.value())) {
                    const double s_t = mech.reactions()[k].forward.translational_exponent;
                    const double t_c = std::pow(t_l, s_t) * std::pow(s.t_e, 1.0 - s_t);
                    forward[k.value()] = std::exp(forward_log_rate(mech, k.value(), t_c)) *
                                         concentration_product(mech, k.value(), n, true);
                }
                // The reference evaluates these equilibrium backward rates at
                // T; the formula at T_e.
                if (contains(backward_at_electron, k.value())) {
                    backward[k.value()] = equilibrium_backward_rate(mech, gas, k.value(), s.t_e) *
                                          concentration_product(mech, k.value(), n, false);
                }
                // Forward rates contain no constant that differs between the
                // implementations, and the formulas compute the same
                // expression in another order. Backward rates contain R_u
                // through the Gibbs energies, corrected above, and differ by
                // rounding amplified by the exponent of K_c (up to about
                // 100).
                expect(relative(got.forward[k.value()], forward[k.value()], 1.0e-12))
                  << s.name << k.value() << got.forward[k.value()] << forward[k.value()];
                expect(relative(got.backward[k.value()], backward[k.value()], 1.0e-10))
                  << s.name << k.value() << got.backward[k.value()] << backward[k.value()];
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
            for (const auto& c : exchange.chemistry_vibration()) {
                cv_scale += boltzmann * c.dissociation_temperature.value() *
                            (got.forward[c.reaction.value()] + got.backward[c.reaction.value()]) /
                            2.0;
            }
            for (const canary::ModeIndex m : canary::mode_indices(0, mode_count)) {
                const double vt = s.vt[m.value()] + dq_vt[m.value()].value();
                const double ve = s.ve[m.value()] + dq_ve[m.value()].value();
                expect(std::abs(got.vt[m.value()] - vt) <= 5.0e-5 * std::abs(s.vt[m.value()]))
                  << s.name << m.value() << got.vt[m.value()] << vt;
                expect(std::abs(got.cv[m.value()] - s.cv[m.value()]) <= 5.0e-5 * cv_scale)
                  << s.name << m.value() << got.cv[m.value()] << s.cv[m.value()];
                expect(std::abs(got.et[m.value()] - s.et[m.value()]) <=
                       5.0e-5 * std::abs(s.et[m.value()]))
                  << s.name << m.value() << got.et[m.value()] << s.et[m.value()];
                expect(std::abs(got.ve[m.value()] - ve) <= 5.0e-5 * std::abs(s.ve[m.value()]))
                  << s.name << m.value() << got.ve[m.value()] << ve;
                const double scale = std::abs(s.vt[m.value()]) + cv_scale +
                                     std::abs(s.et[m.value()]) + std::abs(s.ve[m.value()]) +
                                     std::abs(got.ce[m.value()]);
                const double total = s.total[m.value()] + got.ce[m.value()] +
                                     dq_vt[m.value()].value() + dq_ve[m.value()].value();
                expect(std::abs(got.total[m.value()] - total) <= 5.0e-5 * scale)
                  << s.name << m.value() << got.total[m.value()] << total;
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
        for (const auto k : canary::reaction_indices(0, reaction_count)) {
            const auto f = with_gas.reactions()[k].forward.other_mode;
            const auto b = with_gas.reactions()[k].backward.other_mode;
            expect(!f || *f == canary::ModeIndex {0});
            expect(!b || *b == canary::ModeIndex {0});
        }
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(one.model, fractions(s.y));
            const std::array t_m {canary::ModeTemperature<double, si> {s.t_v}};
            const auto modes =
              canary::rates_of_progress(with_gas, one.model, mix, Rho {s.rho}, Temp {s.temp}, t_m);
            const auto two = canary::rates_of_progress(with_gas,
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
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                expect(modes.forward()[k] == two.forward()[k] &&
                       modes.backward()[k] == two.backward()[k])
                  << s.name << k.value();
                expect(three.forward()[k] == two.forward()[k] &&
                       three.backward()[k] == two.backward()[k])
                  << s.name << k.value();
            }
            const auto omega_modes =
              canary::production_rates(with_gas, one.model, mix, Rho {s.rho}, Temp {s.temp}, t_m);
            const auto omega_two = canary::production_rates(with_gas,
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

    "multi-temperature kinetics: explicit physical temperature requirements"_test = [&] {
        using namespace canary;
        const auto split = make_multi_temperature_mixture(air11_species,
                                                          air11_rotational_degrees,
                                                          air11_electronic_levels,
                                                          split_modes,
                                                          si_scales);
        const auto rate = [](double s, std::optional<PhysicalTemperatureInput> other = {}) {
            return ArrheniusInput {.A = 1.0e15,
                                   .n = -1.0,
                                   .activation_temperature = 50000.0,
                                   .translational_exponent = s,
                                   .other_temperature = other};
        };
        const auto base =
          reaction("N2 + O2 <=> N + N + O2",
                   coefficients({{.name = "N2", .value = 1.0}, {.name = "O2", .value = 1.0}}),
                   coefficients({{.name = "N", .value = 2.0}, {.name = "O2", .value = 1.0}}),
                   rate(0.5),
                   BackwardRateModel::equilibrium,
                   {});
        const auto build = [&](const ReactionInput<3>& r, const auto& model) {
            const MechanismInput<1, 3> input {.temperature_lower = 300.0,
                                              .temperature_upper = 30000.0,
                                              .reactions = {r}};
            const auto mechanism =
              make_reaction_mechanism(input, reaction_species(air11_species), si_scales, p0, model);
            return mechanism.reactions()[ReactionIndex {0}];
        };
        const auto rejects = [&](const auto& input, const auto& model, std::string_view reason) {
            const auto message = error_message([&] { (void)build(input, model); });
            expect(message.starts_with("reaction 0 (N2 + O2 <=> N + N + O2): ")) << message;
            expect(message.find(reason) != std::string::npos) << message;
        };
        for (const bool reorder_modes : {false, true}) {
            auto modes = split_modes;
            if (reorder_modes) { std::swap(modes[0], modes[1]); }
            const auto reordered = make_multi_temperature_mixture(air11_species,
                                                                  air11_rotational_degrees,
                                                                  air11_electronic_levels,
                                                                  modes,
                                                                  si_scales);
            for (const bool reorder_reactants : {false, true}) {
                auto input = base;
                if (reorder_reactants) {
                    std::swap(input.reactants.entries[0], input.reactants.entries[1]);
                }
                rejects(input,
                        reordered.model,
                        "forward rate: other temperature requirement is missing");
                rejects(input, gas, "forward rate: other temperature requirement is missing");
                input.forward = rate(1.0);
                expect(!build(input, reordered.model).forward.other_mode);
                for (const double s : {0.5, 1.0}) {
                    input.forward = rate(s, VibrationTemperatureInput {"N2"});
                    expect(build(input, reordered.model).forward.other_mode ==
                           OptionalModeIndex {ModeIndex {reorder_modes ? 1U : 0U}});
                    input.forward = rate(s, VibrationTemperatureInput {"O2"});
                    expect(build(input, reordered.model).forward.other_mode ==
                           OptionalModeIndex {ModeIndex {reorder_modes ? 0U : 1U}});
                }
                input.backward_model = BackwardRateModel::arrhenius;
                input.equilibrium_temperature.reset();
                input.backward = rate(0.5);
                rejects(input,
                        reordered.model,
                        "backward rate: other temperature requirement is missing");
                input.backward = rate(0.5, VibrationTemperatureInput {"N2"});
                expect(build(input, reordered.model).backward.other_mode ==
                       OptionalModeIndex {ModeIndex {reorder_modes ? 1U : 0U}});
                input.backward = rate(1.0);
                expect(!build(input, reordered.model).backward.other_mode);
            }
        }
        const std::array<PhysicalTemperatureInput, 5> requirements {
          TranslationTemperatureInput {},
          RotationTemperatureInput {},
          VibrationTemperatureInput {"NO"},
          BoundElectronicTemperatureInput {"N"},
          FreeElectronTemperatureInput {}};
        const std::array<OptionalModeIndex, 5> selected {{{},
                                                          {},
                                                          OptionalModeIndex {ModeIndex {1}},
                                                          OptionalModeIndex {ModeIndex {2}},
                                                          OptionalModeIndex {ModeIndex {2}}}};
        for (std::size_t i = 0; i < requirements.size(); ++i) {
            auto input = base;
            input.forward.other_temperature = requirements[i];
            input.equilibrium_temperature = requirements[i];
            const auto result = build(input, split.model);
            expect(result.forward.other_mode == selected[i]);
            expect(result.backward.other_mode == selected[i]);
            expect(result.equilibrium_at_other_temperature == selected[i].has_value());
            input.backward_model = BackwardRateModel::arrhenius;
            input.equilibrium_temperature.reset();
            input.backward = rate(0.25, requirements[i]);
            expect(build(input, split.model).backward.other_mode == selected[i]);
        }
        auto invalid = base;
        invalid.forward = rate(1.0);
        invalid.equilibrium_temperature.reset();
        rejects(invalid, split.model, "equilibrium temperature requirement is missing");
        invalid.equilibrium_temperature = TranslationTemperatureInput {};
        invalid.backward.other_temperature = FreeElectronTemperatureInput {};
        rejects(invalid, split.model, "an equilibrium backward rate uses equilibrium_temperature");
        invalid.backward = rate(1.0);
        invalid.backward_model = BackwardRateModel::arrhenius;
        rejects(invalid,
                split.model,
                "equilibrium temperature is only valid for an equilibrium backward rate");
        for (const double s : {0.5, 1.0}) {
            for (const auto& [source, reason] : std::array {
                   std::pair {PhysicalTemperatureInput {VibrationTemperatureInput {nullptr}},
                              "temperature requirement needs a species"},
                   std::pair {PhysicalTemperatureInput {VibrationTemperatureInput {""}},
                              "temperature requirement needs a species"},
                   std::pair {PhysicalTemperatureInput {VibrationTemperatureInput {"Ar"}},
                              "unknown temperature species"},
                   std::pair {PhysicalTemperatureInput {VibrationTemperatureInput {"N"}},
                              "vibration temperature requires a molecular heavy species"},
                   std::pair {PhysicalTemperatureInput {VibrationTemperatureInput {"N2+"}},
                              "temperature component is not assigned to a mode"},
                   std::pair {PhysicalTemperatureInput {BoundElectronicTemperatureInput {"e-"}},
                              "bound electronic temperature requires a heavy species"}}) {
                invalid = base;
                invalid.forward = rate(s, source);
                rejects(invalid, split.model, reason);
            }
        }
        const auto one = make_multi_temperature_mixture(air11_species,
                                                        air11_rotational_degrees,
                                                        air11_electronic_levels,
                                                        single_mode,
                                                        si_scales);
        rejects(base, one.model, "other temperature requirement is missing");
        auto input = base;
        input.forward = rate(0.5, FreeElectronTemperatureInput {});
        expect(build(input, one.model).forward.other_mode == OptionalModeIndex {ModeIndex {0}});
        const auto air5 = make_air5_two_temperature<double>(si_scales);
        const auto build5 = [&](PhysicalTemperatureInput requirement) {
            auto r = base;
            r.forward = rate(0.5, requirement);
            const MechanismInput<1, 3> mechanism {.temperature_lower = 300.0,
                                                  .temperature_upper = 30000.0,
                                                  .reactions = {r}};
            return make_reaction_mechanism(mechanism,
                                           reaction_species(air5_species),
                                           si_scales,
                                           p0,
                                           air5.model);
        };
        expect(error_message([&] {
                   (void)build5(VibrationTemperatureInput {"N"});
               }).find("vibration temperature requires a molecular heavy species") !=
               std::string::npos);
        const auto atomic = build5(BoundElectronicTemperatureInput {"N"});
        expect(atomic.reactions()[ReactionIndex {0}].forward.other_mode ==
               OptionalModeIndex {ModeIndex {0}});
        expect(error_message([&] {
                   (void)build5(FreeElectronTemperatureInput {});
               }).find("free-electron temperature requires a free-electron species") !=
               std::string::npos);
        static_assert(AcceptsScalarReactionTemperature<std::remove_cvref_t<decltype(one.model)>>);
        for (const auto& source : requirements) {
            if (std::holds_alternative<TranslationTemperatureInput>(source) ||
                std::holds_alternative<RotationTemperatureInput>(source)) {
                continue;
            }
            input = base;
            input.forward.other_temperature = source;
            rejects(input,
                    TemperatureQueryModel {OptionalModeIndex {ModeIndex {2}}},
                    "temperature mode is outside the model mode range");
            rejects(input,
                    TemperatureQueryModel {},
                    "temperature component is not assigned to a mode");
        }
        // Requirement validation uses the original input exponent before storage narrowing.
        const UnitScales<float, si> float_scales {};
        const auto float_gas = make_air11_three_temperature(float_scales);
        auto narrow = base;
        narrow.forward.translational_exponent = std::nextafter(1.0, 0.0);
        expect(static_cast<float>(narrow.forward.translational_exponent) == 1.0_f);
        const auto build_float = [&](const auto& r) {
            const MechanismInput<1, 3> input {.temperature_lower = 300.0,
                                              .temperature_upper = 30000.0,
                                              .reactions = {r}};
            return make_reaction_mechanism(input,
                                           reaction_species(air11_species),
                                           float_scales,
                                           Pressure<float, si> {101325.0F},
                                           float_gas.model);
        };
        expect(error_message([&] {
                   (void)build_float(narrow);
               }).find("other temperature requirement is missing") != std::string::npos);
        narrow.forward.other_temperature = VibrationTemperatureInput {"N2"};
        const auto narrowed = build_float(narrow);
        expect(narrowed.reactions()[ReactionIndex {0}].forward.other_mode ==
               OptionalModeIndex {ModeIndex {0}});
        auto metadata = air11_species;
        metadata[4].metadata.role = SpeciesRole::unspecified;
        auto referenced = base;
        referenced.forward.other_temperature = VibrationTemperatureInput {"NO"};
        const MechanismInput<1, 3> unused_species {.temperature_lower = 300.0,
                                                   .temperature_upper = 30000.0,
                                                   .reactions = {referenced}};
        expect(error_message([&] {
                   (void)make_reaction_mechanism(unused_species,
                                                 reaction_species(metadata),
                                                 si_scales,
                                                 p0,
                                                 split.model);
               }).find("species role must be specified") != std::string::npos);
        expect(
          error_message([&] {
              (void)make_reaction_mechanism(unused_species, reaction_species(metadata), si_scales);
          }).find("species role must be specified") != std::string::npos);
        const MechanismInput<1, 3> collapse {.temperature_lower = 300.0,
                                             .temperature_upper = 30000.0,
                                             .reactions = {base}};
        const auto single =
          make_reaction_mechanism(collapse, reaction_species(air11_species), si_scales);
        expect(!single.reactions()[ReactionIndex {0}].forward.other_mode);
    };

    "physical temperature identities preserve independent derivatives at equal values"_test = [&] {
        using namespace canary;
        const std::array<EnergyModeInput<1>, 3> components {
          {{.count = 1, .components = {{{.species = "N2", .component = ModeComponent::vibration}}}},
           {.count = 1, .components = {{{.species = "N", .component = ModeComponent::electronic}}}},
           {.count = 1,
            .components = {{{.species = "e-", .component = ModeComponent::electron}}}}}};
        for (const bool permute : {false, true}) {
            auto assigned = components;
            if (permute) { std::swap(assigned[0], assigned[2]); }
            const auto split = make_multi_temperature_mixture(air11_species,
                                                              air11_rotational_degrees,
                                                              air11_electronic_levels,
                                                              assigned,
                                                              si_scales);
            auto r =
              reaction("N2 <=> N + N",
                       coefficients({{"N2", 1.0}}),
                       coefficients({{"N", 2.0}}),
                       ArrheniusInput {.A = 1e9,
                                       .n = 1.0,
                                       .activation_temperature = 0.0,
                                       .translational_exponent = 0.0,
                                       .other_temperature = VibrationTemperatureInput {"N2"}},
                       BackwardRateModel::equilibrium,
                       {});
            r.equilibrium_temperature = BoundElectronicTemperatureInput {"N"};
            for (const bool explicit_reverse : {false, true}) {
                if (explicit_reverse) {
                    r.backward_model = BackwardRateModel::arrhenius;
                    r.equilibrium_temperature.reset();
                    r.backward = {.A = 1e9,
                                  .n = 1.0,
                                  .activation_temperature = 0.0,
                                  .translational_exponent = 0.0,
                                  .other_temperature = FreeElectronTemperatureInput {}};
                }
                const MechanismInput<1, 3> input {.temperature_lower = 300.0,
                                                  .temperature_upper = 30000.0,
                                                  .reactions = {r}};
                const auto mechanism = make_reaction_mechanism(input,
                                                               reaction_species(air11_species),
                                                               si_scales,
                                                               p0,
                                                               split.model);
                std::array<D, 11> y {};
                y[0] = D(0.5);
                y[2] = D(0.5);
                const auto mix = mixture(split.model, fractions(y));
                for (std::size_t column = 0; column < 3; ++column) {
                    std::array<ModeTemperature<D, si>, 3> temperatures {};
                    for (std::size_t m = 0; m < 3; ++m) {
                        temperatures[m] =
                          ModeTemperature<D, si> {D(8000.0, m == column ? 1.0 : 0.0)};
                    }
                    const auto q = rates_of_progress(mechanism,
                                                     split.model,
                                                     mix,
                                                     Density<D, si> {D(0.01)},
                                                     Temperature<D, si> {D(8000.0)},
                                                     temperatures);
                    const auto f = q.forward()[ReactionIndex {0}].value().dual();
                    const auto b = q.backward()[ReactionIndex {0}].value().dual();
                    const std::size_t forward_mode = permute ? 2 : 0;
                    const std::size_t backward_mode = explicit_reverse ? (permute ? 0 : 2) : 1;
                    expect(column == forward_mode ? f != 0.0 : f == 0.0);
                    expect(column == backward_mode ? b != 0.0 : b == 0.0);
                }
            }
        }
    };

    "multi-temperature kinetics: translational rates ignore other mode logarithms"_test = [] {
        const canary::ArrheniusRate<double, si> rate {
          .A = 1.0,
          .log_A = 0.0,
          .n = 2.0,
          .activation_temperature = Temp {3.0},
          .translational_exponent = 1.0,
          .other_mode = canary::OptionalModeIndex {canary::ModeIndex {0}}};
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
                expect(a.backward()[canary::ReactionIndex {k}] ==
                       b.backward()[canary::ReactionIndex {k}])
                  << temp << t_v << k;
            }
        }
        // T_e is not clamped: 40000 K is above the upper limit 30000 K.
        const auto hot = rates(12000.0, 5000.0, 40000.0);
        const auto limit = rates(12000.0, 5000.0, 30000.0);
        expect(!(hot.backward()[canary::ReactionIndex {38}] ==
                 limit.backward()[canary::ReactionIndex {38}]));
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
                const double dual = q.backward()[canary::ReactionIndex {k}].value().dual();
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

    "multi-temperature kinetics: cross-sections and the relaxation-time fit follow their formulas"_test =
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

          const auto& fit = exchange.vibration_electron()[0].relaxation_time;
          const auto f = [&](double t) { return canary::detail::power_exponential_log(fit, t); };
          // The input fit, ln(p_e tau / (1 atm s)), and the model fit in SI,
          // ln(p_e tau / (1 Pa s)).
          const auto input_fit = [](double t) {
              const double u = 10000.0 / t;
              return -24.91199 + (4.475475 * std::log(t / 10000.0)) + (3.192608 * u) -
                     (0.08318516 * u * u) + (0.002015418 * u * u * u);
          };
          for (const double t :
               {300.0, 1160.4518121550082, 7000.0, 12000.0, 58022.590607750404, 2.0e5}) {
              expect(relative(f(t), input_fit(t) + std::log(atmosphere), 1.0e-14)) << t;
          }
          // Values of Eq. (19) of J.-H. Lee, AIAA Paper 85-1035, from the rate
          // coefficients of its Table 2 at 0.1, 0.6, 1.0, 2.0 and 5.0 eV, as
          // ln(p_e tau / (1 atm s)); the fit follows them to 0.8 % of tau.
          const std::array<std::array<double, 2>, 5> lee {
            {
              {1160.4518121550082, -11.927111},
              {6962.710872930048, -22.110617},
              {11604.518121550082, -21.552285},
              {23209.036243100163, -19.784929},
              {58022.590607750404, -16.493319},
            },
          };
          for (const auto& [t, value] : lee) {
              expect(std::abs(input_fit(t) - value) <= 0.008) << t << input_fit(t) << value;
          }

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
          // V-E of N2 from the relaxation-time fit.
          const ReferenceState& s = reference_states[0];
          const auto mix_s = canary::mixture(gas, fractions(s.y));
          const auto state = cell_state(gas, s.y, s.rho, s.temp, s.t_v, s.t_e);
          const double q =
            canary::vibration_electron_source(exchange, gas, mix_s, state)[0].value();
          const double n_e =
            s.rho * mix_s.fractions[electron].value() / built.species.mass[electron].value();
          const double p_e = n_e * boltzmann * s.t_e;
          const double inverse_tau = (p_e / atmosphere) * std::exp(-input_fit(s.t_e));
          const auto e_te = canary::species_mode_energy(gas, canary::ModeIndex {0}, Temp {s.t_e});
          const auto e_tv = canary::species_mode_energy(gas, canary::ModeIndex {0}, Temp {s.t_v});
          const double expected =
            s.rho * mix_s.fractions[0].value() * (e_te[0].value() - e_tv[0].value()) * inverse_tau;
          // exp(-f) with f of about -20: a few ulp.
          expect(relative(q, expected, 1.0e-13)) << q << expected;
      };

    "multi-temperature kinetics: V-E without electrons or molecules is zero when the fit overflows"_test =
      [&] {
          // The constant -1000 makes exp(-f) overflow.
          auto electrons = canary::air11_kim_jo_exchange_electrons;
          expect(std::string_view {electrons.vibration_electron[0].species} == "N2");
          electrons.vibration_electron[0].relaxation_time.constant = -1000.0;
          const auto overflowing = build_exchange(canary::air11_kim_jo_exchange, electrons);
          const ReferenceState& s = reference_states[0];
          const double t_e = 12000.0;
          Fractions y = s.y;
          y[electron] = 0.0;
          const auto no_electrons =
            canary::vibration_electron_source(overflowing,
                                              gas,
                                              canary::mixture(gas, fractions(y)),
                                              cell_state(gas, y, s.rho, s.temp, s.t_v, t_e));
          expect(no_electrons[0].value() == 0.0_d && no_electrons[1].value() == 0.0_d)
            << no_electrons[0].value();
          // With electrons and no N2 (the only V-E species) both sources are
          // exactly zero; the product 0 * inf would give NaN.
          y = s.y;
          y[0] = 0.0;
          const auto mix = canary::mixture(gas, fractions(y));
          const auto state = cell_state(gas, y, s.rho, s.temp, s.t_v, t_e);
          const auto got = canary::vibration_electron_source(overflowing, gas, mix, state);
          const auto expected = canary::vibration_electron_source(exchange, gas, mix, state);
          expect(got[0].value() == 0.0_d && got[1].value() == 0.0_d) << got[1].value();
          expect(expected[0].value() == 0.0_d && expected[1].value() == 0.0_d);
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
            const double net = q.forward()[canary::ReactionIndex {38}].value() -
                               q.backward()[canary::ReactionIndex {38}].value() +
                               q.forward()[canary::ReactionIndex {39}].value() -
                               q.backward()[canary::ReactionIndex {39}].value();
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
                  for (const auto k : canary::reaction_indices(0, reaction_count)) {
                      const double both = q.forward()[k].value() + q.backward()[k].value();
                      mass_scale += m * std::abs(coefficient_change(mech, k.value(), i)) * both;
                      charge_scale += std::abs(z * coefficient_change(mech, k.value(), i)) * both;
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
        // Two reference states, and a state with T_e near the minimum of the
        // vibration-electron relaxation time.
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
            for (const auto& c : exchange.chemistry_vibration()) {
                cv_scale = cv_scale +
                           T(boltzmann * c.dissociation_temperature.value() *
                             (q.forward()[c.reaction].value() + q.backward()[c.reaction].value()));
            }
            T ce_scale = T(0.0);
            for (const auto& c : exchange.electron_impact_ionisation()) {
                ce_scale = ce_scale + T(c.ionisation_energy * (q.forward()[c.reaction].value() +
                                                               q.backward()[c.reaction].value()));
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
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                expect(relative(to_si(q_s.forward()[k]), q.forward()[k].value(), 1.0e-12))
                  << s.name << k.value();
                expect(relative(to_si(q_s.backward()[k]), q.backward()[k].value(), 1.0e-10))
                  << s.name << k.value();
            }
            // The SI rates converted, so that only the exchange is compared.
            canary::RatesOfProgress<double, scaled, reaction_count> q_converted {};
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                q_converted.forward()[k] =
                  canary::convert<scaled>(q.forward()[k], reference_scales);
                q_converted.backward()[k] =
                  canary::convert<scaled>(q.backward()[k], reference_scales);
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
            // more. The constant of the V-E relaxation-time fit is shifted
            // by the logarithm of the scales.
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

    "multi-temperature kinetics: the coulomb cross-section needs a singly charged partner and "
    "changes sign at the weak-coupling limit"_test = [&] {
        const std::size_t ion = (*canary::detail::species_index(names, "N+")).value();
        canary::ElectronExchangeInput<1, 0, 0> electrons {};
        electrons.electron_translation[0] =
          canary::air11_kim_jo_exchange_electrons.electron_translation[5];
        const auto build = [&](const auto& species) {
            return canary::make_multi_temperature_exchange(no_vibration,
                                                           electrons,
                                                           names,
                                                           species,
                                                           mech,
                                                           gas,
                                                           si_scales);
        };
        const double e = built.species.charge[ion].value();
        expect(fatal(e == -built.species.charge[electron].value()));

        auto doubly_charged = built.species;
        doubly_charged.charge[ion] = canary::Charge<double, si> {2.0 * e};
        const std::string message = error_message([&] { (void)build(doubly_charged); });
        expect(message == "electron-translation pair 0 (N+): a coulomb partner must carry one "
                          "elementary charge, got 2 electron charges")
          << message;

        // A wrong electron charge is reported as such, not as a wrong partner.
        const double inf = std::numeric_limits<double>::infinity();
        for (const double charge : {0.0, e, std::numeric_limits<double>::quiet_NaN(), inf, -inf}) {
            auto species = built.species;
            species.charge[electron] = canary::Charge<double, si> {charge};
            const std::string electron_message = error_message([&] { (void)build(species); });
            expect(electron_message == "electron-translation pair 0 (N+): the electron charge "
                                       "must be finite and negative")
              << charge << electron_message;
        }

        // The cross-section does not depend on the sign of the charge.
        auto negatively_charged = built.species;
        negatively_charged.charge[ion] = canary::Charge<double, si> {-e};
        const auto pair = build(built.species);
        const auto pair_negative = build(negatively_charged);
        {
            const Fractions& y = reference_states[0].y;
            const auto mix = canary::mixture(gas, fractions(y));
            const auto state = cell_state(gas, y, 1.0e-3, 12000.0, 5000.0, 8000.0);
            const auto positive = canary::electron_translation_source(pair, mix, state);
            const auto negative = canary::electron_translation_source(pair_negative, mix, state);
            expect(positive[1].value() != 0.0_d);
            expect(negative[0].value() == positive[0].value());
            expect(negative[1].value() == positive[1].value());
        }

        // N+ and electrons in the ratio of a neutral plasma, with T > T_e, so
        // that Q_ET has the sign of ln(Lambda). The logarithm is zero at the
        // electron number density coulomb_logarithm_factor T_e^3.
        const double m_e = built.species.mass[electron].value();
        const double m_ion = built.species.mass[ion].value();
        Fractions y {};
        y[ion] = m_ion / (m_ion + m_e);
        y[electron] = m_e / (m_ion + m_e);
        const double t_e = 1.0e4;
        const double limit = pair.coulomb_logarithm_factor() * t_e * t_e * t_e;
        const auto source = [&]<typename T>(double n_e) {
            std::array<T, species_count> y_t {};
            for (std::size_t i = 0; i < species_count; ++i) { y_t[i] = T(y[i]); }
            const auto mix = canary::mixture(gas, fractions(y_t));
            const double rho = n_e * m_e / y[electron];
            const auto state = cell_state(gas, y_t, T(rho), T(12000.0), T(5000.0), T(t_e));
            return canary::ScalarTraits<T>::value(
              canary::electron_translation_source(pair, mix, state)[1].value());
        };
        expect(source.operator()<double>(0.999 * limit) > 0.0);
        expect(source.operator()<double>(1.001 * limit) < 0.0);
        expect(source.operator()<D>(0.999 * limit) > 0.0);
        expect(source.operator()<D>(1.001 * limit) < 0.0);
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
            e.vibration_electron[0].relaxation_time.exponent = v;
            starts(base, e, "vibration-electron exchange 0 (N2): relaxation-time parameters");
            e = electrons;
            e.vibration_electron[0].relaxation_time.inverse_coefficients[2] = v;
            starts(base, e, "vibration-electron exchange 0 (N2): relaxation-time parameters");
        }
        for (const double reference : {0.0, -10000.0}) {
            e = electrons;
            e.vibration_electron[0].relaxation_time.reference_temperature = reference;
            starts(
              base,
              e,
              "vibration-electron exchange 0 (N2): relaxation-time reference temperature must be");
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
        expect(!vibration_only.electron() && !vibration_only.electron_mode());
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
        expect(arbitrary_names.electron() == reference.electron());
        expect(arbitrary_names.electron_impact_ionisation()[0].species == source);
        expect(arbitrary_names.electron_impact_ionisation()[0].ionisation_energy ==
               reference.electron_impact_ionisation()[0].ionisation_energy);

        auto spectator_input = canary::air11_kim_jo;
        auto& spectator_reaction = spectator_input.reactions[38];
        spectator_reaction.reactants.entries[spectator_reaction.reactants.count++] = {"N2", 2};
        spectator_reaction.products.entries[spectator_reaction.products.count++] = {"N2", 2};
        const auto with_spectator =
          canary::make_reaction_mechanism(spectator_input,
                                          canary::reaction_species(canary::air11_species),
                                          si_scales,
                                          p0,
                                          gas);
        expect(error_message([&] { (void)build(built.species, with_spectator, names); }) ==
               "no exception");

        const auto rejects_record = [&](const auto& record, auto reason) {
            const auto message = error_message([&] {
                canary::detail::check_ionisation_participants(record,
                                                              names,
                                                              built.species,
                                                              source,
                                                              carrier,
                                                              canary::ReactionIndex {38});
            });
            expect(message.find(reason) != std::string::npos) << message;
        };
        auto changed = with_spectator.reactions()[canary::ReactionIndex {38}];
        changed.products()[spectator] = 3;
        rejects_record(changed, "no other net species changes");
        changed.products()[spectator] = 1;
        rejects_record(changed, "no other net species changes");
        for (const int count : {0, 1, 3}) {
            changed = mech.reactions()[canary::ReactionIndex {38}];
            changed.products()[carrier] = count;
            rejects_record(changed, "net electron production +1");
        }
        changed = mech.reactions()[canary::ReactionIndex {38}];
        changed.reactants()[carrier] = 0;
        changed.products()[carrier] = 1;
        rejects_record(changed, "incident electron");
        changed = mech.reactions()[canary::ReactionIndex {38}];
        changed.products()[source] = 1;
        rejects_record(changed, "must consume one");
        for (const int count : {0, 2}) {
            changed = mech.reactions()[canary::ReactionIndex {38}];
            changed.products()[ion] = count;
            rejects_record(changed, "must produce one ion");
        }
        // A reaction may conserve charge and nuclei while consuming two neutrals.
        // This exchange model requires one neutral consumed per reaction.
        auto unsupported_input = canary::air11_kim_jo;
        auto& unsupported_reaction = unsupported_input.reactions[38];
        for (int i = 0; i < unsupported_reaction.reactants.count; ++i) {
            auto& coefficient = unsupported_reaction.reactants.entries[i];
            if (std::string_view(coefficient.name) == "O") { coefficient.value = 2; }
        }
        for (int i = 0; i < unsupported_reaction.products.count; ++i) {
            auto& coefficient = unsupported_reaction.products.entries[i];
            coefficient.value = std::string_view(coefficient.name) == "e-" ? 3 : 2;
        }
        const auto unsupported =
          canary::make_reaction_mechanism(unsupported_input,
                                          canary::reaction_species(canary::air11_species),
                                          si_scales,
                                          p0,
                                          gas);
        rejects(built.species, unsupported, "must consume one");
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
              expect(exchange.electron() == canary::electron_species_index(built.model));
              for (const auto& coupling : exchange.electron_impact_ionisation()) {
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
