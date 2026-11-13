// A function that requires a gas model without padding bytes rejects a model
// that has them.
#include <canary/data/air5.hpp>
#include <canary/layout.hpp>
#include <canary/thermo.hpp>

namespace
{
template <canary::GasModelType G>
    requires canary::PaddingFree<G>
int image(const G& /*gas*/)
{ return 0; }
}  // namespace

#ifdef CANARY_COMPILE_FAIL
int bad = image(canary::make_air5<float>(canary::UnitScales<float, canary::si> {}).model);
#endif
int good = image(canary::make_air5<double>(canary::UnitScales<double, canary::si> {}).model);
