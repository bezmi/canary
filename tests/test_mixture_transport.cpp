#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_transport.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <canary/transport/constant.hpp>
#include <canary/transport/mixture.hpp>
#include <canary/transport/prandtl.hpp>
#include <canary/transport/transport.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <numbers>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;
using Air = canary::Air5<double, si>;
using AirTwoTemperature = canary::Air5TwoTemperature<double, si>;
using Wilke = canary::Air5WilkeTransport<double, si>;
using GuptaYos = canary::Air5GuptaYosTransport<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;

constexpr std::size_t species_count = 5;
using Fractions = std::array<double, species_count>;
using Row = std::array<double, species_count>;
using Matrix = std::array<Row, species_count>;

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();
constexpr double boltzmann = canary::constants::boltzmann<double>;
constexpr double avogadro = canary::constants::avogadro<double>;

const canary::UnitScales<double, si> si_scales {};
const canary::UnitScales<double, scaled> reference_scales {{
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
}};

constexpr Fractions cold_air {0.767, 0.233, 0.0, 0.0, 0.0};
constexpr Fractions dissociated {0.3, 0.1, 0.3, 0.25, 0.05};

template <typename T> auto fractions(const std::array<T, species_count>& y)
{
    return canary::as_kind<canary::kind::mass_fraction>(
      Kokkos::mdspan<const T, Extents>(y.data(), species_count));
}

template <canary::UnitSystem U, typename T> auto mode_temperatures(const std::array<T, 1>& t_ve)
{
    return canary::as_kind<canary::kind::mode_temperature, U>(
      Kokkos::mdspan<const T, Extents>(t_ve.data(), 1));
}

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

// A state of a reference implementation of the two models for five-species
// air with the data of air5_transport.hpp, evaluated from (p, T, y) and,
// for the two-temperature model, T_ve: viscosity, conductivity, the mode
// conductivity k_ve and the binary diffusion coefficients, all in SI. Its
// Wilke model mixes the conductivities by their ratio and leaves D_ii out
// (nan); its Gupta-Yos model gives k = k_tr + k_rot and k_ve separately.
struct ReferenceState {
    Fractions y;
    double p;
    double temp;
    double t_ve;
    double mu;
    double k;
    double k_ve;
    Matrix d;
};

// clang-format off
const std::array<ReferenceState, 8> wilke_reference {
 {
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .p = 100000.0, .temp = 300.0, .t_ve = quiet_nan,
    .mu = 1.8503076010484828e-05, .k = 0.026257139263519774, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 2.1394794606812487e-05, 2.984027646982789e-05, 3.3236308153695796e-05, 2.1026025083461284e-05},
      {2.1394794606812487e-05, quiet_nan, 3.0282826607893495e-05, 3.3790319277682844e-05, 2.1012585249546513e-05},
      {2.984027646982789e-05, 3.0282826607893495e-05, quiet_nan, 4.469646780053729e-05, 2.9508146846331078e-05},
      {3.3236308153695796e-05, 3.3790319277682844e-05, 4.469646780053729e-05, quiet_nan, 3.283267671464772e-05},
      {2.1026025083461284e-05, 2.1012585249546513e-05, 2.9508146846331078e-05, 3.283267671464772e-05, quiet_nan},
    }, },
  },
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .p = 100000.0, .temp = 1000.0, .t_ve = quiet_nan,
    .mu = 4.3179551043161135e-05, .k = 0.06599255107214451, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 0.0001671363555095312, 0.00022880274387195116, 0.0002560448785797809, 0.0001634588303880179},
      {0.0001671363555095312, quiet_nan, 0.0002331192609297815, 0.00026142842852664084, 0.00016415053207962386},
      {0.00022880274387195116, 0.0002331192609297815, quiet_nan, 0.00034018680511855096, 0.00022625611300363242},
      {0.0002560448785797809, 0.00026142842852664084, 0.00034018680511855096, quiet_nan, 0.00025293539474890135},
      {0.0001634588303880179, 0.00016415053207962386, 0.00022625611300363242, 0.00025293539474890135, quiet_nan},
    }, },
  },
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .p = 10000.0, .temp = 5000.0, .t_ve = quiet_nan,
    .mu = 0.00013156347801856192, .k = 0.24765541922452486, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 0.02408286769304058, 0.03291635983551894, 0.036847264168200215, 0.023541379908250126},
      {0.02408286769304058, quiet_nan, 0.033546246740007495, 0.03763469818217849, 0.023652636996684774},
      {0.03291635983551894, 0.033546246740007495, quiet_nan, 0.04892133061153857, 0.03254999264685999},
      {0.036847264168200215, 0.03763469818217849, 0.04892133061153857, quiet_nan, 0.03639978022406241},
      {0.023541379908250126, 0.023652636996684774, 0.03254999264685999, 0.03639978022406241, quiet_nan},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 10000.0, .temp = 500.0, .t_ve = quiet_nan,
    .mu = 3.8587134761617286e-05, .k = 0.06858418270175409, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 0.0005194425580161682, 0.0007165948525965409, 0.0008003455574505192, 0.0005090225018523826},
      {0.0005194425580161682, quiet_nan, 0.0007289091422104047, 0.0008157241279781566, 0.0005101629266906712},
      {0.0007165948525965409, 0.0007289091422104047, quiet_nan, 0.001068700778145924, 0.0007086189754684153},
      {0.0008003455574505192, 0.0008157241279781566, 0.001068700778145924, quiet_nan, 0.0007906259271114447},
      {0.0005090225018523826, 0.0005101629266906712, 0.0007086189754684153, 0.0007906259271114447, quiet_nan},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 10000.0, .temp = 3000.0, .t_ve = quiet_nan,
    .mu = 0.00010224671618862985, .k = 0.20102314766304671, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 0.010334854575552513, 0.01412563831516619, 0.015812535764406357, 0.010102483349023009},
      {0.010334854575552513, quiet_nan, 0.01439594608619361, 0.016150452537816776, 0.01015022657620255},
      {0.01412563831516619, 0.01439594608619361, quiet_nan, 0.020993968720565918, 0.013968416483122774},
      {0.015812535764406357, 0.016150452537816776, 0.020993968720565918, quiet_nan, 0.015620503709098886},
      {0.010102483349023009, 0.01015022657620255, 0.013968416483122774, 0.015620503709098886, quiet_nan},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 1000.0, .temp = 8000.0, .t_ve = quiet_nan,
    .mu = 0.00021423048012021451, .k = 0.44557288425290853, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 0.5245070174216977, 0.7168939322889665, 0.802506116582829, 0.5127138146017373},
      {0.5245070174216977, quiet_nan, 0.7306124024394297, 0.8196558460658439, 0.5151369115761195},
      {0.7168939322889665, 0.7306124024394297, quiet_nan, 1.0654703390702844, 0.7089147263302336},
      {0.802506116582829, 0.8196558460658439, 1.0654703390702844, quiet_nan, 0.7927602477822596},
      {0.5127138146017373, 0.5151369115761195, 0.7089147263302336, 0.7927602477822596, quiet_nan},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 1000.0, .temp = 12000.0, .t_ve = quiet_nan,
    .mu = 0.00029498584296530986, .k = 0.6573831429864334, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 1.0265403792859058, 1.403070968959045, 1.5706270954118937, 1.00345927933169},
      {1.0265403792859058, quiet_nan, 1.4299201112654802, 1.6041917365384981, 1.008201650366793},
      {1.403070968959045, 1.4299201112654802, quiet_nan, 2.08528826051479, 1.38745444365202},
      {1.5706270954118937, 1.6041917365384981, 2.08528826051479, quiet_nan, 1.5515529409721949},
      {1.00345927933169, 1.008201650366793, 1.38745444365202, 1.5515529409721949, quiet_nan},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 1000.0, .temp = 20000.0, .t_ve = quiet_nan,
    .mu = 0.000353152884524178, .k = 0.818422154403889, .k_ve = quiet_nan,
    .d = {{
      {quiet_nan, 2.392102315186896, 3.2695151412876657, 3.6599638816386544, 2.338317433703766},
      {2.392102315186896, quiet_nan, 3.3320805276747882, 3.738178102304012, 2.349368374281901},
      {3.2695151412876657, 3.3320805276747882, quiet_nan, 4.859256368735783, 3.233124632842109},
      {3.6599638816386544, 3.738178102304012, 4.859256368735783, quiet_nan, 3.615516210688608},
      {2.338317433703766, 2.349368374281901, 3.233124632842109, 3.615516210688608, quiet_nan},
    }, },
  },
 },
};

const std::array<ReferenceState, 7> gupta_yos_reference {
 {
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .p = 100000.0, .temp = 300.0, .t_ve = 300.0,
    .mu = 2.0176093931151995e-05, .k = 0.0285383813414462, .k_ve = 8.10348503206131e-05,
    .d = {{
      {1.8232938950801226e-05, 2.726261956309996e-05, 2.9124272933611707e-05, 2.4896700108572597e-05, 2.2357065616328285e-05},
      {2.726261956309996e-05, 2.7370063847004554e-05, 2.962588540670919e-05, 3.3273242672395306e-05, 2.732996062692988e-05},
      {2.9124272933611707e-05, 2.962588540670919e-05, 1.1542788704668242e-05, 3.261179505665383e-05, 2.9317067478712306e-05},
      {2.4896700108572597e-05, 3.3273242672395306e-05, 3.261179505665383e-05, 1.1449046865747248e-05, 2.8625262225667332e-05},
      {2.2357065616328285e-05, 2.732996062692988e-05, 2.9317067478712306e-05, 2.8625262225667332e-05, 2.471528880164355e-05},
    }, },
  },
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .p = 10000.0, .temp = 5000.0, .t_ve = 1000.0,
    .mu = 0.00012462142036472716, .k = 0.18672352304353157, .k_ve = 0.024922109914359534,
    .d = {{
      {0.02708084399999634, 0.023786034317222975, 0.04164667315421273, 0.03725994289124046, 0.025343606878248288},
      {0.023786034317222975, 0.023678050747403246, 0.03952761746905348, 0.038790272690040664, 0.023798099839032742},
      {0.04164667315421273, 0.03952761746905348, 0.010527227953202175, 0.05960863115845299, 0.04045103521969657},
      {0.03725994289124046, 0.038790272690040664, 0.05960863115845299, 0.010483594216460463, 0.03819256031251923},
      {0.025343606878248288, 0.023798099839032742, 0.04045103521969657, 0.03819256031251923, 0.02459595829276787},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 10000.0, .temp = 10000.0, .t_ve = 3000.0,
    .mu = 0.0002456283717571788, .k = 0.47270096323775423, .k_ve = 0.03397081269206418,
    .d = {{
      {0.09540210364716518, 0.08008572286968697, 0.14827622698810666, 0.1323731098963816, 0.08724624792842603},
      {0.08008572286968697, 0.0784992676550771, 0.1378458424519865, 0.1323546683238446, 0.07956380962862236},
      {0.14827622698810666, 0.1378458424519865, 0.032261789277865195, 0.21257928043985977, 0.14244528814698146},
      {0.1323731098963816, 0.1323546683238446, 0.21257928043985977, 0.03216757169689863, 0.1331900577058379},
      {0.08724624792842603, 0.07956380962862236, 0.14244528814698146, 0.1331900577058379, 0.08346685937122587},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 10000.0, .temp = 3000.0, .t_ve = 6000.0,
    .mu = 9.436775412795857e-05, .k = 0.18000715359661293, .k_ve = 0.01347691793555784,
    .d = {{
      {0.010779791435087737, 0.010004365737671968, 0.016532251488428833, 0.014764300424858766, 0.010374810375132314},
      {0.010004365737671968, 0.010039034858603142, 0.015917799293305477, 0.015920053337237155, 0.010044852823488862},
      {0.016532251488428833, 0.015917799293305477, 0.004621194427122991, 0.023284450013386065, 0.0161792468959619},
      {0.014764300424858766, 0.015920053337237155, 0.023284450013386065, 0.004598157502703757, 0.015380170839492368},
      {0.010374810375132314, 0.010044852823488862, 0.0161792468959619, 0.015380170839492368, 0.010221512195963709},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 1000.0, .temp = 15000.0, .t_ve = 8000.0,
    .mu = 0.0003424738548873458, .k = 0.6596937882778369, .k_ve = 0.04917843852428201,
    .d = {{
      {2.0028026966645225, 1.66326854150959, 3.1435925567748493, 2.796046268020408, 1.821503554139132},
      {1.66326854150959, 1.6116910077417923, 2.885324187933214, 2.7409808085775205, 1.643654981594326},
      {3.1435925567748493, 2.885324187933214, 0.6220677030440672, 4.462883291058378, 2.999416765964573},
      {2.796046268020408, 2.7409808085775205, 4.462883291058378, 0.6207246514355459, 2.787871496557686},
      {1.821503554139132, 1.643654981594326, 2.999416765964573, 2.787871496557686, 1.7336699139921656},
    }, },
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .p = 1000.0, .temp = 20000.0, .t_ve = 20000.0,
    .mu = 0.00043491818661630727, .k = 0.8378067309972396, .k_ve = 0.05097652738838327,
    .d = {{
      {3.3972511016402422, 2.8196380424585126, 5.378469033043512, 4.76604407375331, 3.088652118960538},
      {2.8196380424585126, 2.7070434512829697, 4.890538313690424, 4.615162968603995, 2.7743886936445388},
      {5.378469033043512, 4.890538313690424, 0.9918363407526711, 7.5462223708603835, 5.10607141947871},
      {4.76604407375331, 4.615162968603995, 7.5462223708603835, 0.9902508587222588, 4.725358912728324},
      {3.088652118960538, 2.7743886936445388, 5.10607141947871, 4.725358912728324, 2.9332880044221676},
    }, },
  },
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .p = 10000.0, .temp = 2000.0, .t_ve = 2000.0,
    .mu = 6.492851769821179e-05, .k = 0.09603883959277462, .k_ve = 0.022323863937722753,
    .d = {{
      {0.005210488959510057, 0.005118458964865422, 0.00799785463314579, 0.007117820522631503, 0.005161682013767614},
      {0.005118458964865422, 0.005158402618878566, 0.007784452705201628, 0.007917495532225406, 0.005148432502784589},
      {0.00799785463314579, 0.007784452705201628, 0.002406996893110492, 0.01102175291551528, 0.007871466832617144},
      {0.007117820522631503, 0.007917495532225406, 0.01102175291551528, 0.0023934827489379296, 0.007521530476157129},
      {0.005161682013767614, 0.005148432502784589, 0.007871466832617144, 0.007521530476157129, 0.005160537161230234},
    }, },
  },
 },
};
// clang-format on

// D_ij of canary over D_ij of the reference Wilke model, which computes D_ij
// from the constant 0.00266 with p in atm.
constexpr double reference_diffusion_ratio = 0.9882309643094598;

// canary results at one state, in SI: viscosity, conductivity, the mode
// conductivity (quiet_nan without modes), D_ij and D_i.
struct Result {
    double mu;
    double k;
    double k_ve;
    Matrix d;
    Row d_mix;
};

// The state of gas at (p, T) with the mass fractions y and, for a model with
// modes, T_ve; the inputs are converted from SI to the units of the gas.
template <typename G, typename T>
auto state_at(const G& gas,
              const canary::mixture_t<G, T>& mix,
              T p,
              T temp,
              T t_ve,
              const canary::UnitScales<double, typename G::units>& scales)
{
    using U = G::units;
    const auto p_q = canary::convert<U>(canary::Pressure<T, si> {p}, scales);
    const auto temp_q = canary::convert<U>(canary::Temperature<T, si> {temp}, scales);
    if constexpr (std::same_as<canary::state_t<G, T>, canary::ThermoState<T, U>>) {
        return canary::thermo_state(gas, mix, p_q, temp_q).state;
    } else {
        const std::array<T, 1> modes {
          canary::convert<U>(canary::Temperature<T, si> {t_ve}, scales).value(),
        };
        return canary::thermo_state(gas, mix, p_q, temp_q, mode_temperatures<U>(modes)).state;
    }
}

template <typename M, typename G>
Result evaluate(const M& model,
                const G& gas,
                const Fractions& y,
                double p,
                double temp,
                double t_ve,
                const canary::UnitScales<double, typename G::units>& scales)
{
    const auto mix = canary::mixture(gas, fractions(y));
    const auto state = state_at(gas, mix, p, temp, t_ve, scales);
    const auto coefficients = canary::transport(model, gas, mix, state);
    const auto d = canary::binary_diffusion_coefficients(model, gas, mix, state);
    const auto d_mix = canary::mixture_diffusion_coefficients(model, gas, mix, state);
    Result out {};
    out.mu = canary::convert<si>(coefficients.viscosity, scales).value();
    out.k = canary::convert<si>(coefficients.conductivity, scales).value();
    out.k_ve = quiet_nan;
    if constexpr (requires { coefficients.mode_conductivity; }) {
        out.k_ve = canary::convert<si>(coefficients.mode_conductivity[0], scales).value();
    }
    for (std::size_t i = 0; i < species_count; ++i) {
        for (std::size_t j = 0; j < species_count; ++j) {
            out.d[i][j] = canary::convert<si>(d[i][j], scales).value();
        }
        out.d_mix[i] = canary::convert<si>(d_mix[i], scales).value();
    }
    return out;
}

// A CEA fit in SI at temperature t: the segment k of the input, or by
// default the first segment that contains t, with t clamped to the ends.
template <std::size_t S>
double cea_value(const canary::TransportCurveInput<S>& curve, double t, std::size_t k = S)
{
    const auto n = static_cast<std::size_t>(curve.segments);
    const double tc = std::clamp(t, curve.breaks[0], curve.breaks[n]);
    if (k == S) {
        k = 0;
        while (k + 1 < n && tc > curve.breaks[k + 1]) { ++k; }
    }
    const std::array<double, 4>& c = curve.coefficients[k];
    return std::exp((c[0] * std::log(tc)) + (c[1] / tc) + (c[2] / (tc * tc)) + c[3]);
}

// Mole fractions of the mass fractions y, from the molar masses.
Row mole_fractions(const Fractions& y)
{
    Row x {};
    double sum = 0;
    for (std::size_t i = 0; i < species_count; ++i) {
        x[i] = y[i] / canary::air5_species[i].molar_mass;
        sum += x[i];
    }
    for (double& v : x) { v /= sum; }
    return x;
}

// Wilke's rule for the species values v with phi_ij from the species values
// w (w = v for the viscosity and the conductivity-ratio rule).
double wilke_rule(const Row& x, const Row& v, const Row& w)
{
    double total = 0;
    for (std::size_t i = 0; i < species_count; ++i) {
        const double m_i = canary::air5_species[i].molar_mass;
        double sum = 0;
        for (std::size_t j = 0; j < species_count; ++j) {
            const double m_j = canary::air5_species[j].molar_mass;
            const double a = 1.0 + (std::sqrt(w[i] / w[j]) * std::pow(m_j / m_i, 0.25));
            sum += x[j] * a * a / std::sqrt(8.0 * (1.0 + (m_i / m_j)));
        }
        total += x[i] * v[i] / sum;
    }
    return total;
}

// Delta(2)_ij in SI at T from the fit of the pair.
double delta_2(std::size_t i, std::size_t j, double t)
{
    const std::string_view a = canary::air5_species[i].name;
    const std::string_view b = canary::air5_species[j].name;
    for (const auto& pair : canary::air5_collision_integrals) {
        if ((pair.first == a && pair.second == b) || (pair.first == b && pair.second == a)) {
            const double m_i = canary::air5_species[i].molar_mass / avogadro;
            const double m_j = canary::air5_species[j].molar_mass / avogadro;
            const double mu = m_i * m_j / (m_i + m_j);
            const double l = std::log(t);
            const auto& f = pair.pi_omega_22;
            const double pi_omega = std::exp((f.A * l * l * l) + (f.B * l * l) + (f.C * l) + f.D);
            return (16.0 / 5.0) * std::sqrt(2.0 * mu / (std::numbers::pi * boltzmann * t)) *
                   pi_omega;
        }
    }
    return quiet_nan;
}

// A copy of air5_species whose molecules have the constant heat capacity
// cp / R = cp_over_r at every temperature.
std::array<canary::Nasa9Species<4>, species_count> constant_cp_species(double cp_over_r)
{
    auto species = canary::air5_species;
    for (const std::size_t i : {std::size_t {0}, std::size_t {1}, std::size_t {4}}) {
        canary::Nasa9Data<4> fit {};
        fit.segments = 1;
        fit.breaks[0] = 200.0;
        fit.breaks[1] = 20000.0;
        fit.a[0][2] = cp_over_r;
        species[i].fit = fit;
    }
    return species;
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

using SutherlandPrandtl =
  canary::Transport<canary::ConstantViscosity<double, si>, canary::ConstantPrandtl<double>>;
using CaloricallyPerfectAir = canary::CaloricallyPerfectMixture<double, si, species_count>;

}  // namespace

static_assert(canary::HasNormalisedComposition<Air, double, 5>);
static_assert(canary::HasNormalisedComposition<AirTwoTemperature, D, 5>);
static_assert(!canary::HasNormalisedComposition<Air, double, 4>);
static_assert(!canary::HasNormalisedComposition<CaloricallyPerfectAir, double, 5>);
static_assert(canary::HasSpeciesIsobaricSpecificHeat<Air, double, 5>);
static_assert(canary::HasSpeciesIsobaricSpecificHeat<AirTwoTemperature, D, 5>);
static_assert(canary::HasSpeciesIsobaricSpecificHeat<AirTwoTemperature, float, 5>);
static_assert(!canary::HasSpeciesIsobaricSpecificHeat<Air, double, 4>);
static_assert(!canary::HasSpeciesIsobaricSpecificHeat<CaloricallyPerfectAir, double, 5>);
static_assert(canary::HasVibrationalElectronicMode<AirTwoTemperature>);
static_assert(!canary::HasVibrationalElectronicMode<Air>);
static_assert(
  !canary::HasVibrationalElectronicMode<canary::MultiTemperatureMixture<double, si, 11, 1>>);
static_assert(canary::TransportModel<Wilke, Air, double>);
static_assert(canary::TransportModel<Wilke, Air, D>);
static_assert(canary::TransportModel<Wilke, Air, float>);
static_assert(canary::TransportModel<Wilke, AirTwoTemperature, D>);
static_assert(canary::TransportModel<GuptaYos, Air, D>);
static_assert(canary::TransportModel<GuptaYos, AirTwoTemperature, double>);
static_assert(canary::TransportModel<GuptaYos, AirTwoTemperature, float>);
static_assert(canary::TransportModel<canary::Air5GuptaYosTransport<float, scaled>,
                                     canary::Air5TwoTemperature<float, scaled>,
                                     double>);
static_assert(!canary::TransportModel<GuptaYos, CaloricallyPerfectAir, double>);
static_assert(!canary::TransportModel<GuptaYos, canary::Air5<double, scaled>, double>);
// The single mode of a multi-temperature model may hold any energy components,
// so it has no Gupta-Yos mode conductivity; Wilke gives a zero mode
// conductivity for any mode.
static_assert(canary::TransportModel<Wilke,
                                     canary::MultiTemperatureMixture<double, si, species_count, 1>,
                                     double>);
static_assert(!canary::TransportModel<GuptaYos,
                                      canary::MultiTemperatureMixture<double, si, species_count, 1>,
                                      double>);
// Constant-Prandtl conductivity requires a single-temperature cp.
static_assert(!canary::TransportModel<SutherlandPrandtl, AirTwoTemperature, double>);
static_assert(std::is_trivially_copyable_v<Wilke>);
static_assert(std::is_trivially_copyable_v<GuptaYos>);
static_assert(
  std::is_trivially_copyable_v<canary::MultiTemperatureTransportCoefficients<D, si, 1>>);
static_assert(std::same_as<decltype(canary::transport(
                             std::declval<const Wilke&>(),
                             std::declval<const AirTwoTemperature&>(),
                             std::declval<const canary::mixture_t<AirTwoTemperature, double>&>(),
                             std::declval<const canary::state_t<AirTwoTemperature, double>&>())),
                           canary::MultiTemperatureTransportCoefficients<double, si, 1>>);
static_assert(canary::air5_species_transport.size() == 5);
static_assert(canary::air5_collision_integrals.size() == 15);

void test_mixture_transport()
{
    using namespace boost::ut;
    const auto air = canary::make_air5<double>(si_scales);
    const auto air_two_temperature = canary::make_air5_two_temperature<double>(si_scales);
    const Wilke wilke = canary::make_air5_wilke_transport(air.species, si_scales);
    const Wilke wilke_viscosity_ratio =
      canary::make_air5_wilke_transport(air.species,
                                        si_scales,
                                        canary::ConductivityMixing::viscosity_ratio);
    const GuptaYos gupta_yos = canary::make_air5_gupta_yos_transport(air.species, si_scales);

    "single-species transport retains self diffusion and zero mixture diffusion"_test = [] {
        const auto check = []<std::floating_point P, canary::UnitSystem U>(
                             const canary::UnitScales<P, U>& scales) {
            const std::array species {canary::air5_species[0]};
            const std::array input {canary::air5_species_transport[0]};
            const std::array pairs {canary::air5_collision_integrals[0]};
            const std::array<std::string_view, 1> names {"N2"};
            const auto built = canary::make_thermally_perfect_mixture(species, scales);
            const auto wilke = canary::make_wilke_transport(input, names, built.species, scales);
            const auto gupta = canary::make_gupta_yos_transport(
              pairs, names, built.species, std::array<int, 1> {2}, scales);
            const auto state_check = [&]<typename T>() {
                using V = canary::ScalarTraits<T>::value_type;
                std::array<T, 1> weights {T(V(1))};
                const auto composition = canary::as_kind<canary::kind::mass_fraction>(
                  Kokkos::mdspan<const T, Kokkos::extents<std::size_t, 1>>(weights.data()));
                const auto mix = canary::mixture(built.model, composition);
                const auto p = canary::convert<U>(canary::Pressure<T, si> {T(V(10000))}, scales);
                const auto temp = canary::convert<U>(canary::Temperature<T, si> {T(V(3000))}, scales);
                const auto state = canary::thermo_state(built.model, mix, p, temp);
                expect(fatal(state.status == canary::Status::Valid));
                const auto transport_check = [&](const auto& model) {
                    const auto all = canary::binary_diffusion_coefficients(
                      model, built.model, mix, state.state);
                    const auto mixture = canary::mixture_diffusion_coefficients(
                      model, built.model, mix, state.state);
                    expect(canary::ScalarTraits<T>::value(all[0][0].value()) > V(0));
                    expect(canary::math::isfinite(all[0][0].value()));
                    expect(canary::ScalarTraits<T>::value(mixture[0].value()) == V(0));
                    if constexpr (requires { mixture[0].value().dual(); }) {
                        expect(mixture[0].value().dual() == V(0));
                    }
                    const auto c = canary::transport(model, built.model, mix, state.state);
                    expect(canary::math::isfinite(c.viscosity.value()));
                    expect(canary::math::isfinite(c.conductivity.value()));
                };
                transport_check(wilke);
                transport_check(gupta);
            };
            state_check.template operator()<P>();
            state_check.template operator()<canary::test::Dual<P>>();
        };
        check(canary::UnitScales<float, si> {});
        check(canary::UnitScales<double, si> {});
        check(canary::UnitScales<float, scaled> {
          {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}});
        check(reference_scales);
    };

    "mixture transport: Wilke values equal the reference"_test = [&] {
        for (const ReferenceState& s : wilke_reference) {
            const Result got = evaluate(wilke, air.model, s.y, s.p, s.temp, s.t_ve, si_scales);
            // The reference uses the same fits and rules; the evaluations
            // differ in the order of operations only, by a few ulp of the
            // exponents (up to about 50) of the fits.
            expect(relative(got.mu, s.mu, 1.0e-12)) << s.temp << got.mu << s.mu;
            expect(relative(got.k, s.k, 1.0e-12)) << s.temp << got.k << s.k;
            for (std::size_t i = 0; i < species_count; ++i) {
                for (std::size_t j = 0; j < species_count; ++j) {
                    if (i == j) {
                        expect(std::isfinite(got.d[i][j]) && got.d[i][j] > 0.0);
                        continue;
                    }
                    // The constant of the reference is rounded; the ratio is
                    // stated to 1e-12.
                    expect(relative(got.d[i][j], reference_diffusion_ratio * s.d[i][j], 1.0e-12))
                      << s.temp << i << j << got.d[i][j] / s.d[i][j];
                }
            }
        }
    };

    "mixture transport: Gupta-Yos values equal the reference"_test = [&] {
        for (const ReferenceState& s : gupta_yos_reference) {
            const Result got =
              evaluate(gupta_yos, air_two_temperature.model, s.y, s.p, s.temp, s.t_ve, si_scales);
            // The reference uses rounded CGS constants (k_B = 1.38066e-16
            // erg/K, R_u = 8.31451 J/(mol K), 1.546e-20), which change the
            // results by up to about 1e-4.
            expect(relative(got.mu, s.mu, 3.0e-4)) << s.temp << got.mu << s.mu;
            expect(relative(got.k, s.k, 3.0e-4)) << s.temp << got.k << s.k;
            expect(relative(got.k_ve, s.k_ve, 3.0e-4)) << s.temp << s.t_ve << got.k_ve << s.k_ve;
            for (std::size_t i = 0; i < species_count; ++i) {
                for (std::size_t j = 0; j < species_count; ++j) {
                    expect(relative(got.d[i][j], s.d[i][j], 3.0e-4))
                      << s.temp << i << j << got.d[i][j] << s.d[i][j];
                }
            }
        }
    };

    "mixture transport: Wilke values of a pure species are the species values"_test = [&] {
        for (std::size_t i = 0; i < species_count; ++i) {
            Fractions y {};
            y[i] = 1.0;
            const auto& data = canary::air5_species_transport[i];
            for (const double temp : {150.0, 300.0, 750.0, 2500.0, 7000.0, 12000.0, 20000.0}) {
                const double mu = cea_value(data.viscosity, temp);
                const double k = cea_value(data.conductivity, temp);
                for (const Wilke* model : {&wilke, &wilke_viscosity_ratio}) {
                    const Result got =
                      evaluate(*model, air.model, y, 1.0e5, temp, quiet_nan, si_scales);
                    // The fits evaluated in another order of operations:
                    // a few ulp of exponents up to about 50.
                    expect(relative(got.mu, mu, 1.0e-13)) << i << temp << got.mu << mu;
                    expect(relative(got.k, k, 1.0e-13)) << i << temp << got.k << k;
                }
            }
        }
    };

    "mixture transport: the fits take the lower segment at a break and are constant beyond the ends"_test =
      [&] {
          const auto pure = [&](std::size_t i, double temp) {
              Fractions y {};
              y[i] = 1.0;
              return evaluate(wilke, air.model, y, 1.0e5, temp, quiet_nan, si_scales);
          };
          for (std::size_t i = 0; i < species_count; ++i) {
              const auto& data = canary::air5_species_transport[i];
              const auto n = static_cast<std::size_t>(data.viscosity.segments);
              for (std::size_t k = 0; k + 1 < n; ++k) {
                  const double t_break = data.viscosity.breaks[k + 1];
                  const Result got = pure(i, t_break);
                  // As for the pure-species values.
                  expect(relative(got.mu, cea_value(data.viscosity, t_break, k), 1.0e-13))
                    << i << t_break;
                  expect(relative(got.k, cea_value(data.conductivity, t_break, k), 1.0e-13))
                    << i << t_break;
              }
              // N and O start at 1000 K, the molecules at 200 K; every fit
              // ends at 15000 K. The value beyond an end is the value at the
              // end, computed by the same function.
              const double low = data.viscosity.breaks[0];
              const double high = data.viscosity.breaks[n];
              expect(low == (i == 2 || i == 3 ? 1000.0 : 200.0) && high == 15000.0);
              const Result at_low = pure(i, low);
              const Result at_high = pure(i, high);
              for (const double below : {0.5 * low, 0.99 * low}) {
                  const Result got = pure(i, below);
                  expect(got.mu == at_low.mu && got.k == at_low.k) << i << below;
              }
              for (const double above : {1.01 * high, 2.0 * high}) {
                  const Result got = pure(i, above);
                  expect(got.mu == at_high.mu && got.k == at_high.k) << i << above;
              }
          }
          // Beyond the ends the derivative with respect to T is zero.
          const std::array<D, species_count> y {1.0, 0.0, 0.0, 0.0, 0.0};
          const auto mix = canary::mixture(air.model, fractions(y));
          for (const double temp : {100.0, 30000.0}) {
              const auto c = canary::thermo_state(air.model,
                                                  mix,
                                                  canary::Pressure<D, si> {D(1.0e5)},
                                                  canary::Temperature<D, si> {D(temp, 1.0)});
              const auto coefficients = canary::transport(wilke, air.model, mix, c.state);
              expect(coefficients.viscosity.value().dual() == 0.0 &&
                     coefficients.conductivity.value().dual() == 0.0)
                << temp;
          }
      };

    "mixture transport: the viscosity-ratio rule equals a direct evaluation of the formula"_test =
      [&] {
          for (const Fractions& y : {cold_air, dissociated}) {
              const Row x = mole_fractions(y);
              for (const double temp : {300.0, 2500.0, 9000.0}) {
                  Row mu {};
                  Row k {};
                  for (std::size_t i = 0; i < species_count; ++i) {
                      mu[i] = cea_value(canary::air5_species_transport[i].viscosity, temp);
                      k[i] = cea_value(canary::air5_species_transport[i].conductivity, temp);
                  }
                  const Result ratio = evaluate(wilke_viscosity_ratio,
                                                air.model,
                                                y,
                                                1.0e5,
                                                temp,
                                                quiet_nan,
                                                si_scales);
                  const Result own =
                    evaluate(wilke, air.model, y, 1.0e5, temp, quiet_nan, si_scales);
                  // Sums of positive terms of a few operations each: tens
                  // of ulp.
                  expect(relative(ratio.mu, wilke_rule(x, mu, mu), 1.0e-13)) << temp;
                  expect(relative(ratio.k, wilke_rule(x, k, mu), 1.0e-13)) << temp;
                  expect(relative(own.k, wilke_rule(x, k, k), 1.0e-13)) << temp;
                  expect(ratio.mu == own.mu);
              }
          }
      };

    "mixture transport: the Gupta-Yos viscosity of a pure species is m_i / Delta(2)_ii"_test = [&] {
        for (std::size_t i = 0; i < species_count; ++i) {
            Fractions y {};
            y[i] = 1.0;
            for (const double temp : {300.0, 2500.0, 9000.0, 25000.0}) {
                const Result got = evaluate(gupta_yos, air.model, y, 1.0e5, temp, temp, si_scales);
                const double m = canary::air5_species[i].molar_mass / avogadro;
                // The constants of Delta are folded into the exponent of
                // the fit, of magnitude about 50: a few 1e-15.
                expect(relative(got.mu, m / delta_2(i, i, temp), 1.0e-12)) << i << temp << got.mu;
            }
        }
    };

    "mixture transport: binary, trace and pure mixtures have the limiting diffusion coefficients"_test =
      [&] {
          const auto check = [&](const auto& model, const auto& gas) {
              // N2 and O2 only: D_1 = D_2 = D_12; for the absent species,
              // D_i = 1 / sum_j x_j / D_ij.
              const Fractions binary {0.6, 0.4, 0.0, 0.0, 0.0};
              const Row x = mole_fractions(binary);
              const Result r = evaluate(model, gas, binary, 1.0e4, 4000.0, 3000.0, si_scales);
              // D_i is a ratio of two sums of two terms: a few ulp.
              expect(relative(r.d_mix[0], r.d[0][1], 8.0 * eps)) << r.d_mix[0] << r.d[0][1];
              expect(relative(r.d_mix[1], r.d[0][1], 8.0 * eps)) << r.d_mix[1] << r.d[0][1];
              for (std::size_t i = 2; i < species_count; ++i) {
                  const double expected = 1.0 / ((x[0] / r.d[i][0]) + (x[1] / r.d[i][1]));
                  expect(relative(r.d_mix[i], expected, 8.0 * eps)) << i;
              }
              // A trace of N: D_N differs from 1 / sum_j x_j / D_Nj by
              // the factor 1 - x_N.
              Fractions trace = dissociated;
              trace[2] = 1.0e-9;
              const Row xt = mole_fractions(trace);
              const Result t = evaluate(model, gas, trace, 1.0e4, 4000.0, 3000.0, si_scales);
              double sum = 0;
              for (std::size_t j = 0; j < species_count; ++j) {
                  if (j != 2) { sum += xt[j] / t.d[2][j]; }
              }
              expect(relative(t.d_mix[2], 1.0 / sum, 1.0e-8)) << t.d_mix[2] * sum;
              // Pure O2: D_O2 = 0, and every other species diffuses at its
              // binary coefficient with O2.
              const Fractions pure {0.0, 1.0, 0.0, 0.0, 0.0};
              const Result p = evaluate(model, gas, pure, 1.0e4, 4000.0, 3000.0, si_scales);
              expect(p.d_mix[1] == 0.0);
              for (const std::size_t i :
                   {std::size_t {0}, std::size_t {2}, std::size_t {3}, std::size_t {4}}) {
                  expect(relative(p.d_mix[i], p.d[i][1], 4.0 * eps)) << i;
              }
              // D_ij is symmetric.
              for (std::size_t i = 0; i < species_count; ++i) {
                  for (std::size_t j = 0; j < species_count; ++j) {
                      expect(r.d[i][j] == r.d[j][i]) << i << j;
                  }
              }
          };
          check(wilke, air.model);
          check(gupta_yos, air.model);
          check(gupta_yos, air_two_temperature.model);
      };

    "mixture transport: Gupta-Yos k on the single-temperature model is k + k_ve at T_ve = T"_test =
      [&] {
          for (const Fractions& y : {cold_air, dissociated}) {
              for (const double temp : {300.0, 2000.0, 8000.0, 20000.0}) {
                  const Result one =
                    evaluate(gupta_yos, air.model, y, 1.0e4, temp, temp, si_scales);
                  const Result two =
                    evaluate(gupta_yos, air_two_temperature.model, y, 1.0e4, temp, temp, si_scales);
                  // The same terms; the closures of the two models give the
                  // same T and p.
                  expect(relative(one.k, two.k + two.k_ve, 4.0 * eps)) << temp;
                  expect(relative(one.mu, two.mu, 4.0 * eps)) << temp;
                  expect(std::isnan(one.k_ve));
                  expect(two.k_ve > 0.0);
              }
          }
      };

    "mixture transport: k_ve is zero when cp / R <= 7/2 and limited at cp / R >= 9/2"_test = [&] {
        const auto k_ve_of = [&](double cp_over_r, double t_ve) {
            const auto gas = canary::make_two_temperature_mixture(constant_cp_species(cp_over_r),
                                                                  canary::air5_rotational_degrees,
                                                                  si_scales);
            const auto model = canary::make_air5_gupta_yos_transport(gas.species, si_scales);
            return evaluate(model, gas.model, dissociated, 1.0e4, 5000.0, t_ve, si_scales);
        };
        for (const double t_ve : {300.0, 6000.0}) {
            const Result below = k_ve_of(3.4, t_ve);
            expect(below.k_ve == 0.0) << t_ve << below.k_ve;
            // At cp / R = 7/2 the ratio of the thermo model's R_i to
            // k_B / m_i is 1 within rounding, so c_i is at most a few ulp.
            const Result at = k_ve_of(3.5, t_ve);
            expect(at.k_ve <= 1.0e-14 * at.k) << t_ve << at.k_ve;
            // c_i = 1 above 9/2, and 1/2 at 4 within a few ulp of the ratio
            // of the gas constants.
            const Result above = k_ve_of(5.0, t_ve);
            const Result far = k_ve_of(6.0, t_ve);
            const Result middle = k_ve_of(4.0, t_ve);
            expect(above.k_ve == far.k_ve && above.k_ve > 0.0) << t_ve;
            expect(relative(middle.k_ve, 0.5 * above.k_ve, 32.0 * eps)) << t_ve;
            // The other coefficients do not depend on cp.
            expect(below.mu == above.mu && below.k == above.k) << t_ve;
        }
    };

    "mixture transport: an absent species with a non-finite collision fit contributes nothing"_test =
      [&] {
          // The (N, N2) pair: N is absent from cold air. ln Delta = +800
          // overflows and -800 underflows the exponential at every T.
          const auto reference = evaluate(gupta_yos,
                                          air_two_temperature.model,
                                          cold_air,
                                          1.0e4,
                                          4000.0,
                                          3000.0,
                                          si_scales);
          for (const double d : {800.0, -800.0}) {
              auto pairs = canary::air5_collision_integrals;
              expect(std::string_view {pairs[3].first} == "N" &&
                     std::string_view {pairs[3].second} == "N2");
              pairs[3].pi_omega_11.D = d;
              pairs[3].pi_omega_22.D = d;
              const auto model =
                canary::make_gupta_yos_transport(pairs,
                                                 canary::detail::air5_species_names,
                                                 air.species,
                                                 canary::air5_rotational_degrees,
                                                 si_scales);
              const auto got = evaluate(model,
                                        air_two_temperature.model,
                                        cold_air,
                                        1.0e4,
                                        4000.0,
                                        3000.0,
                                        si_scales);
              expect(got.mu == reference.mu && got.k == reference.k && got.k_ve == reference.k_ve)
                << d << got.mu << got.k << got.k_ve;
              // D_i of the present species are unchanged; D_N depends on its
              // own fit and is only required to be finite.
              for (std::size_t i = 0; i < species_count; ++i) {
                  if (i == 2) {
                      expect(std::isfinite(got.d_mix[i])) << d << got.d_mix[i];
                  } else {
                      expect(got.d_mix[i] == reference.d_mix[i]) << d << i << got.d_mix[i];
                  }
              }
          }
      };

    "mixture transport: scaled units give the SI results after conversion"_test = [&] {
        const auto air_s = canary::make_air5<double>(reference_scales);
        const auto air_two_temperature_s =
          canary::make_air5_two_temperature<double>(reference_scales);
        const auto wilke_s = canary::make_air5_wilke_transport(air_s.species, reference_scales);
        const auto wilke_ratio_s =
          canary::make_air5_wilke_transport(air_s.species,
                                            reference_scales,
                                            canary::ConductivityMixing::viscosity_ratio);
        const auto gupta_yos_s =
          canary::make_air5_gupta_yos_transport(air_s.species, reference_scales);
        const auto compare = [](const Result& a, const Result& b, double temp) {
            // The constants are folded in double for each unit system; the
            // folded exponents (up to about 50) differ by a few ulp, and the
            // closures solve nothing at (p, T).
            const double bound = 1.0e-12;
            expect(relative(a.mu, b.mu, bound)) << temp << a.mu << b.mu;
            expect(relative(a.k, b.k, bound)) << temp << a.k << b.k;
            expect((std::isnan(a.k_ve) && std::isnan(b.k_ve)) || relative(a.k_ve, b.k_ve, bound))
              << temp << a.k_ve << b.k_ve;
            for (std::size_t i = 0; i < species_count; ++i) {
                for (std::size_t j = 0; j < species_count; ++j) {
                    expect(relative(a.d[i][j], b.d[i][j], bound)) << temp << i << j;
                }
                expect(relative(a.d_mix[i], b.d_mix[i], bound)) << temp << i;
            }
        };
        for (const ReferenceState& s : wilke_reference) {
            compare(evaluate(wilke_s, air_s.model, s.y, s.p, s.temp, quiet_nan, reference_scales),
                    evaluate(wilke, air.model, s.y, s.p, s.temp, quiet_nan, si_scales),
                    s.temp);
            compare(
              evaluate(wilke_ratio_s, air_s.model, s.y, s.p, s.temp, quiet_nan, reference_scales),
              evaluate(wilke_viscosity_ratio, air.model, s.y, s.p, s.temp, quiet_nan, si_scales),
              s.temp);
        }
        for (const ReferenceState& s : gupta_yos_reference) {
            compare(
              evaluate(gupta_yos_s,
                       air_two_temperature_s.model,
                       s.y,
                       s.p,
                       s.temp,
                       s.t_ve,
                       reference_scales),
              evaluate(gupta_yos, air_two_temperature.model, s.y, s.p, s.temp, s.t_ve, si_scales),
              s.temp);
            compare(evaluate(gupta_yos_s, air_s.model, s.y, s.p, s.temp, s.temp, reference_scales),
                    evaluate(gupta_yos, air.model, s.y, s.p, s.temp, s.temp, si_scales),
                    s.temp);
        }
    };

    "mixture transport: dual-number derivatives equal central differences"_test = [&] {
        // Outputs mu, k, k_ve (0 without modes), D_ij for j >= i and D_i, at
        // (p, T, T_ve) with the mass fraction of O2 set to y_o2.
        constexpr std::size_t outputs = 3 + 15 + species_count;
        const auto values =
          [&]<typename T>(const auto& model, const auto& gas, T temp, T t_ve, T y_o2) {
              std::array<T, species_count> y {};
              for (std::size_t i = 0; i < species_count; ++i) { y[i] = T(dissociated[i]); }
              y[1] = y_o2;
              const auto mix = canary::mixture(gas, fractions(y));
              const auto state = state_at(gas, mix, T(1.0e4), temp, t_ve, si_scales);
              const auto c = canary::transport(model, gas, mix, state);
              const auto d = canary::binary_diffusion_coefficients(model, gas, mix, state);
              const auto d_mix = canary::mixture_diffusion_coefficients(model, gas, mix, state);
              std::array<T, outputs> out {};
              out[0] = c.viscosity.value();
              out[1] = c.conductivity.value();
              out[2] = T(0.0);
              if constexpr (requires { c.mode_conductivity; }) {
                  out[2] = c.mode_conductivity[0].value();
              }
              std::size_t n = 3;
              for (std::size_t i = 0; i < species_count; ++i) {
                  for (std::size_t j = i; j < species_count; ++j) { out[n++] = d[i][j].value(); }
              }
              for (std::size_t i = 0; i < species_count; ++i) { out[n++] = d_mix[i].value(); }
              return out;
          };
        // T between the breaks of the fits; T_ve where 7/2 < cp_i / R_i <
        // 9/2 for N2, O2 and NO.
        const double temp = 3000.0;
        const double t_ve = 1500.0;
        const double y_o2 = dissociated[1];
        const auto check = [&](const auto& model, const auto& gas, const char* name) {
            // Inputs T, T_ve and y_O2, and the step of each.
            const std::array<double, 3> inputs {temp, t_ve, y_o2};
            const std::array<double, 3> steps {1.0e-5 * temp, 1.0e-5 * t_ve, 1.0e-6};
            for (std::size_t c = 0; c < 3; ++c) {
                std::array<D, 3> seeded {D(temp), D(t_ve), D(y_o2)};
                seeded[c] = D(inputs[c], 1.0);
                const auto ad = values(model, gas, seeded[0], seeded[1], seeded[2]);
                const double h = steps[c];
                const auto shifted = [&](double sign) {
                    std::array<double, 3> x = inputs;
                    x[c] += sign * h;
                    return values(model, gas, x[0], x[1], x[2]);
                };
                const auto plus = shifted(1.0);
                const auto minus = shifted(-1.0);
                for (std::size_t i = 0; i < outputs; ++i) {
                    const double fd = (plus[i] - minus[i]) / (2.0 * h);
                    // Central differences with a step of 1e-5 of the
                    // temperature or 1e-6 of a mass fraction have a truncation
                    // error of order 1e-10 relative; their rounding error is
                    // a few eps of the value over 2 h, bounded by 1e-13 of
                    // the value over h.
                    const double rounding = 1.0e-13 * std::abs(ad[i].real()) / h;
                    expect(std::abs(fd - ad[i].dual()) <=
                           (1.0e-6 * std::abs(ad[i].dual())) + rounding)
                      << name << c << i << ad[i].dual() << fd;
                }
            }
        };
        check(wilke, air.model, "Wilke");
        check(wilke_viscosity_ratio, air.model, "Wilke, viscosity ratio");
        check(gupta_yos, air.model, "Gupta-Yos");
        check(gupta_yos, air_two_temperature.model, "Gupta-Yos, two temperatures");
    };

    "mixture transport: float states agree with double states within float precision"_test = [&] {
        const auto mix_d = canary::mixture(air_two_temperature.model, fractions(dissociated));
        std::array<float, species_count> y_f {};
        for (std::size_t i = 0; i < species_count; ++i) {
            y_f[i] = static_cast<float>(dissociated[i]);
        }
        const auto mix_f = canary::mixture(air_two_temperature.model, fractions(y_f));
        const auto state_d =
          state_at(air_two_temperature.model, mix_d, 1.0e4, 6000.0, 2000.0, si_scales);
        const auto state_f =
          state_at(air_two_temperature.model, mix_f, 1.0e4F, 6000.0F, 2000.0F, si_scales);
        const auto d = canary::transport(gupta_yos, air_two_temperature.model, mix_d, state_d);
        const auto f = canary::transport(gupta_yos, air_two_temperature.model, mix_f, state_f);
        // Rounding of float inputs and operations, amplified by the
        // exponents of the fits (up to about 50).
        const double bound = 1.0e-5;
        expect(relative(f.viscosity.value(), d.viscosity.value(), bound));
        expect(relative(f.conductivity.value(), d.conductivity.value(), bound));
        expect(relative(f.mode_conductivity[0].value(), d.mode_conductivity[0].value(), bound));
        const auto w_d =
          canary::mixture_diffusion_coefficients(wilke, air_two_temperature.model, mix_d, state_d);
        const auto w_f =
          canary::mixture_diffusion_coefficients(wilke, air_two_temperature.model, mix_f, state_f);
        for (std::size_t i = 0; i < species_count; ++i) {
            expect(relative(w_f[i].value(), w_d[i].value(), bound)) << i;
        }
    };

    "transport curve construction retains long double scale precision"_test = [] {
        const canary::TransportCurveInput<1> input {.segments = 1,
                                                    .breaks = {1.0, 2.0},
                                                    .coefficients = {{{1.0, 1.0, 1.0, 0.0}}}};
        const long double tau = std::nextafter(1.0L, 2.0L);
        const auto curve = canary::make_transport_curve(input, tau, 1.0L);
        expect(curve.breaks[0] == 1.0L / tau && curve.breaks[1] == 2.0L / tau);
        expect(curve.coefficients[0][1] == 1.0L / tau);
        expect(curve.coefficients[0][2] == 1.0L / (tau * tau));
        expect(curve.coefficients[0][3] == std::log(tau));
    };

    "mixture transport: make_transport_curve rejects scales that are not finite and positive"_test =
      [] {
          const canary::TransportCurveInput<3> curve = canary::air5_species_transport[0].viscosity;
          for (const double bad : {0.0, -1.0, quiet_nan, std::numeric_limits<double>::infinity()}) {
              expect(error_message([&] {
                         (void)canary::make_transport_curve(curve, bad, 1.0);
                     }).starts_with("temperature scale must be finite and positive"));
              expect(error_message([&] {
                         (void)canary::make_transport_curve(curve, 1.0, bad);
                     }).starts_with("value scale must be finite and positive"));
          }
      };

    // One test per error case of the Wilke builder.
    using WilkeInput = std::array<canary::SpeciesTransportInput<3>, species_count>;
    struct WilkeCase {
        const char* name;
        void (*change)(WilkeInput&);
        const char* message;
    };
    const std::array wilke_cases {
      WilkeCase {
        .name = "an entry without a name",
        .change = [](WilkeInput& in) { in[2].name = nullptr; },
        .message = "transport data entry 2 has no name",
      },
      WilkeCase {
        .name = "an unknown species",
        .change = [](WilkeInput& in) { in[2].name = "Ar"; },
        .message = "unknown species \"Ar\" in the transport data",
      },
      WilkeCase {
        .name = "a species listed twice",
        .change = [](WilkeInput& in) { in[2].name = "N2"; },
        .message = "species \"N2\" is listed twice in the transport data",
      },
      WilkeCase {
        .name = "a curve with no segments",
        .change = [](WilkeInput& in) { in[0].viscosity.segments = 0; },
        .message = "species N2: viscosity: segment count must be between 1 and 3, got 0",
      },
      WilkeCase {
        .name = "a curve with more segments than S",
        .change = [](WilkeInput& in) { in[0].conductivity.segments = 4; },
        .message = "species N2: conductivity: segment count must be between 1 and 3, got 4",
      },
      WilkeCase {
        .name = "breaks not increasing",
        .change = [](WilkeInput& in) { in[1].conductivity.breaks[2] = 900.0; },
        .message = "species O2: conductivity: segment boundaries must be finite, positive and "
                   "increasing, got 900 at index 2",
      },
      WilkeCase {
        .name = "a break that is not finite",
        .change = [](WilkeInput& in) { in[1].viscosity.breaks[0] = quiet_nan; },
        .message = "species O2: viscosity: segment boundaries must be finite, positive and "
                   "increasing, got nan at index 0",
      },
      WilkeCase {
        .name = "a coefficient that is not finite",
        .change =
          [](WilkeInput& in) {
              in[4].viscosity.coefficients[1][2] = std::numeric_limits<double>::infinity();
          },
        .message = "species NO: viscosity: coefficient 2 of segment 1 is not finite",
      },
      // ln(mu) = 2e5 / T - 5e7 / T^2 + 540 is 290 and 690 at the ends and has
      // its maximum 740 at T = 500 K, above ln(DBL_MAX) = 709.8.
      WilkeCase {
        .name = "a fit that overflows inside a segment",
        .change =
          [](WilkeInput& in) {
              in[0].viscosity = {.segments = 1,
                                 .breaks = {200.0, 1000.0},
                                 .coefficients = {{{0.0, 2.0e5, -5.0e7, 540.0}}}};
          },
        .message = "species N2: viscosity: the fit at temperature 500 in segment 0 must be finite "
                   "and positive, got inf",
      },
      // ln(k) = -2e5 / T + 5e7 / T^2 - 555 is -305 and -705 at the ends and
      // has its minimum -755 at T = 500 K, below the logarithm -745.1 of the
      // smallest positive double.
      WilkeCase {
        .name = "a fit that underflows inside a segment",
        .change =
          [](WilkeInput& in) {
              in[1].conductivity = {.segments = 1,
                                    .breaks = {200.0, 1000.0},
                                    .coefficients = {{{0.0, -2.0e5, 5.0e7, -555.0}}}};
          },
        .message = "species O2: conductivity: the fit at temperature 500 in segment 0 must be "
                   "finite and positive, got 0",
      },
      WilkeCase {
        .name = "a sigma that is not positive",
        .change = [](WilkeInput& in) { in[3].lennard_jones.sigma = 0.0; },
        .message = "species O: Lennard-Jones sigma must be finite and positive, got 0",
      },
      WilkeCase {
        .name = "a sigma that is not finite",
        .change = [](WilkeInput& in) { in[3].lennard_jones.sigma = quiet_nan; },
        .message = "species O: Lennard-Jones sigma must be finite and positive, got nan",
      },
      WilkeCase {
        .name = "an epsilon that is not positive",
        .change = [](WilkeInput& in) { in[3].lennard_jones.epsilon = -80.0; },
        .message = "species O: Lennard-Jones epsilon must be finite and positive, got -80",
      },
      WilkeCase {
        .name = "an epsilon that is not finite",
        .change =
          [](WilkeInput& in) {
              in[3].lennard_jones.epsilon = std::numeric_limits<double>::infinity();
          },
        .message = "species O: Lennard-Jones epsilon must be finite and positive, got inf",
      },
    };
    for (const WilkeCase& c : wilke_cases) {
        test(std::string("mixture transport: the Wilke builder rejects ") + c.name) = [&] {
            WilkeInput input = canary::air5_species_transport;
            c.change(input);
            const std::string message = error_message([&] {
                (void)canary::make_wilke_transport(input,
                                                   canary::detail::air5_species_names,
                                                   air.species,
                                                   si_scales);
            });
            expect(message == c.message) << message;
        };
    }
    "mixture transport: the Wilke builder rejects a missing species"_test = [&] {
        std::array<canary::SpeciesTransportInput<3>, 4> input {};
        std::copy_n(canary::air5_species_transport.begin(), 4, input.begin());
        const std::string message = error_message([&] {
            (void)canary::make_wilke_transport(input,
                                               canary::detail::air5_species_names,
                                               air.species,
                                               si_scales);
        });
        expect(message == "no transport data for species \"NO\"") << message;
    };
    "mixture transport: the Wilke builder rejects an unknown mixing rule"_test = [&] {
        const std::string message = error_message([&] {
            // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange): the invalid value
            const auto invalid = static_cast<canary::ConductivityMixing>(7);
            (void)canary::make_air5_wilke_transport(air.species, si_scales, invalid);
        });
        expect(message == "unknown conductivity mixing rule") << message;
    };

    // One test per error case of the Gupta-Yos builder.
    using Pairs = std::array<canary::CollisionPairInput, 15>;
    struct GuptaYosCase {
        const char* name;
        void (*change)(Pairs&);
        const char* message;
    };
    const std::array gupta_yos_cases {
      GuptaYosCase {
        .name = "an unknown species",
        .change = [](Pairs& p) { p[3].first = "Ar"; },
        .message = "collision pair 3 (Ar, N2): unknown species \"Ar\"",
      },
      GuptaYosCase {
        .name = "a species without a name",
        .change = [](Pairs& p) { p[3].second = nullptr; },
        .message = "collision pair 3 (N, ): a species has no name",
      },
      GuptaYosCase {
        .name = "a pair listed twice",
        .change =
          [](Pairs& p) {
              p[3].first = "N2";
              p[3].second = "O2";
          },
        .message = "collision pair 3 (N2, O2): the pair is listed twice, first as collision pair 1",
      },
      GuptaYosCase {
        .name = "a coefficient that is not finite",
        .change = [](Pairs& p) { p[2].pi_omega_22.C = quiet_nan; },
        .message = "collision pair 2 (O2, O2): coefficient C of pi Omega(2,2) must be finite, got "
                   "nan",
      },
      GuptaYosCase {
        .name = "a coefficient of pi Omega(1,1) that is not finite",
        .change = [](Pairs& p) { p[14].pi_omega_11.A = std::numeric_limits<double>::infinity(); },
        .message = "collision pair 14 (NO, NO): coefficient A of pi Omega(1,1) must be finite, "
                   "got inf",
      },
    };
    for (const GuptaYosCase& c : gupta_yos_cases) {
        test(std::string("mixture transport: the Gupta-Yos builder rejects ") + c.name) = [&] {
            Pairs pairs = canary::air5_collision_integrals;
            c.change(pairs);
            const std::string message = error_message([&] {
                (void)canary::make_gupta_yos_transport(pairs,
                                                       canary::detail::air5_species_names,
                                                       air.species,
                                                       canary::air5_rotational_degrees,
                                                       si_scales);
            });
            expect(message == c.message) << message;
        };
    }
    "mixture transport: the Gupta-Yos builder rejects a missing pair"_test = [&] {
        std::array<canary::CollisionPairInput, 14> pairs {};
        std::copy_n(canary::air5_collision_integrals.begin(), 14, pairs.begin());
        const std::string message = error_message([&] {
            (void)canary::make_gupta_yos_transport(pairs,
                                                   canary::detail::air5_species_names,
                                                   air.species,
                                                   canary::air5_rotational_degrees,
                                                   si_scales);
        });
        expect(message == "no collision integrals for the pair (NO, NO)") << message;
    };
    "mixture transport: the Gupta-Yos builder rejects a rotational degree not 0 or 2"_test = [&] {
        for (const int d : {-1, 1, 3, 4}) {
            const std::string message = error_message([&] {
                (void)canary::make_gupta_yos_transport(canary::air5_collision_integrals,
                                                       canary::detail::air5_species_names,
                                                       air.species,
                                                       std::array {2, 2, d, 0, 2},
                                                       si_scales);
            });
            expect(
              message ==
              std::format("species N: rotational degrees of freedom must be 0 or 2, got {}", d))
              << message;
        }
    };
    "mixture transport: the builders reject species data of another size"_test = [&] {
        auto species = air.species;
        species.n = 4;
        expect(error_message([&] {
                   (void)canary::make_air5_wilke_transport(species, si_scales);
               }) == "species data must have 5 species, got 4");
        expect(error_message([&] {
                   (void)canary::make_air5_gupta_yos_transport(species, si_scales);
               }) == "species data must have 5 species, got 4");
    };
    "mixture transport: the builders reject a species named twice"_test = [&] {
        const std::array<std::string_view, 5> names {"N2", "O2", "N", "N2", "NO"};
        expect(error_message([&] {
                   (void)canary::make_wilke_transport(canary::air5_species_transport,
                                                      names,
                                                      air.species,
                                                      si_scales);
               }) == "species \"N2\" is named twice");
        expect(error_message([&] {
                   (void)canary::make_gupta_yos_transport(canary::air5_collision_integrals,
                                                          names,
                                                          air.species,
                                                          canary::air5_rotational_degrees,
                                                          si_scales);
               }) == "species \"N2\" is named twice");
    };
    "mixture transport: the builders reject a particle mass that is not finite and positive"_test =
      [&] {
          for (const double bad : {0.0, -1.0, quiet_nan, std::numeric_limits<double>::infinity()}) {
              auto species = air.species;
              species.mass[2] = canary::ParticleMass<double, si> {bad};
              const std::string expected =
                std::format("species N: particle mass must be finite and positive, got {}", bad);
              const std::string wilke_message =
                error_message([&] { (void)canary::make_air5_wilke_transport(species, si_scales); });
              expect(wilke_message == expected) << wilke_message;
              const std::string gupta_yos_message = error_message(
                [&] { (void)canary::make_air5_gupta_yos_transport(species, si_scales); });
              expect(gupta_yos_message == expected) << gupta_yos_message;
          }
      };
    "mixture transport: the Wilke builder rejects breaks that are equal in the stored units"_test =
      [] {
          // 1000 K and 1000.00001 K are both 1 in float after division by the
          // temperature unit of 1000 K.
          const canary::UnitScales<float, scaled> scales {
            {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}};
          const auto air_float = canary::make_air5<float>(scales);
          WilkeInput input = canary::air5_species_transport;
          input[0].viscosity.segments = 1;
          input[0].viscosity.breaks[0] = 1000.0;
          input[0].viscosity.breaks[1] = 1000.00001;
          const std::string message = error_message([&] {
              (void)canary::make_wilke_transport(input,
                                                 canary::detail::air5_species_names,
                                                 air_float.species,
                                                 scales);
          });
          expect(message ==
                 "species N2: viscosity: segment boundaries 0 and 1 collapse after conversion")
            << message;
      };
    "mixture transport: the builders reject constants that overflow the parameter type"_test = [] {
        const canary::UnitScales<float, si> scales {};
        const auto air_float = canary::make_air5<float>(scales);
        // The factor of D_ij is proportional to 1 / sigma_ij^2: about 6.7e39
        // in SI for two oxygen atoms with sigma = 1e-31 m, above FLT_MAX.
        WilkeInput input = canary::air5_species_transport;
        input[3].lennard_jones.sigma = 1.0e-31;
        const std::string wilke_message = error_message([&] {
            (void)canary::make_wilke_transport(input,
                                               canary::detail::air5_species_names,
                                               air_float.species,
                                               scales);
        });
        expect(wilke_message.starts_with(
          "species O and O: binary diffusion factor must be finite and positive, got "))
          << wilke_message;
        // The constant D of a collision fit is a term of the cubic
        // coefficient d.
        Pairs pairs = canary::air5_collision_integrals;
        pairs[2].pi_omega_11.D = 1.0e39;
        const std::string gupta_yos_message = error_message([&] {
            (void)canary::make_gupta_yos_transport(pairs,
                                                   canary::detail::air5_species_names,
                                                   air_float.species,
                                                   canary::air5_rotational_degrees,
                                                   scales);
        });
        expect(gupta_yos_message == "collision pair 2 (O2, O2): cubic coefficient d of ln Delta "
                                    "from pi Omega(1,1) must be finite, got 1e+39")
          << gupta_yos_message;
    };
}
