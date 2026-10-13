#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <canary/status.hpp>

namespace
{

using canary::Status;

constexpr std::array all_statuses {
  Status::Valid,
  Status::CompositionViolation,
  Status::NotConverged,
  Status::Invalid,
};

}  // namespace

static_assert(canary::combine(Status::Valid, Status::Valid) == Status::Valid);
static_assert(canary::combine(Status::CompositionViolation, Status::NotConverged) ==
              Status::NotConverged);

void test_status()
{
    using namespace boost::ut;
    using canary::combine;

    "combine is commutative"_test = [] {
        for (const Status a : all_statuses) {
            for (const Status b : all_statuses) { expect(combine(a, b) == combine(b, a)); }
        }
    };

    "combine is associative"_test = [] {
        for (const Status a : all_statuses) {
            for (const Status b : all_statuses) {
                for (const Status c : all_statuses) {
                    expect(combine(combine(a, b), c) == combine(a, combine(b, c)));
                }
            }
        }
    };

    "Valid is the identity and Invalid dominates"_test = [] {
        for (const Status a : all_statuses) {
            expect(combine(Status::Valid, a) == a);
            expect(combine(Status::Invalid, a) == Status::Invalid);
        }
    };

    "combine keeps the more severe status"_test = [] {
        expect(combine(Status::CompositionViolation, Status::Valid) ==
               Status::CompositionViolation);
        expect(combine(Status::NotConverged, Status::CompositionViolation) == Status::NotConverged);
        expect(combine(Status::NotConverged, Status::Invalid) == Status::Invalid);
    };
}
