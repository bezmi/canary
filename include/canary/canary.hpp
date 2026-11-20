// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_CANARY_HPP
#define CANARY_CANARY_HPP

// Thermodynamics, transport, reaction and exchange headers, generated data
// and host builders.
// A thermally perfect device kernel can include
// <canary/models/thermally_perfect.hpp> for close and properties; add
// <canary/thermo.hpp> for the thermo_state front end. The corresponding
// model headers provide two-temperature and multi-temperature calls.
// Construct models on the host and pass them to kernels. Consumer GPU
// annotations and tested CUDA compiler options are documented in config.hpp.

#include <canary/build.hpp>                       // IWYU pragma: export
#include <canary/config.hpp>                      // IWYU pragma: export
#include <canary/data/air11.hpp>                  // IWYU pragma: export
#include <canary/data/air11_kim_jo.hpp>           // IWYU pragma: export
#include <canary/data/air11_kim_jo_exchange.hpp>  // IWYU pragma: export
#include <canary/data/air11_species.hpp>          // IWYU pragma: export
#include <canary/data/air5.hpp>                   // IWYU pragma: export
#include <canary/data/air5_electronic.hpp>        // IWYU pragma: export
#include <canary/data/air5_kim_jo.hpp>            // IWYU pragma: export
#include <canary/data/air5_park.hpp>              // IWYU pragma: export
#include <canary/data/air5_park_exchange.hpp>     // IWYU pragma: export
#include <canary/data/air5_transport.hpp>         // IWYU pragma: export
#include <canary/energy_exchange.hpp>             // IWYU pragma: export
#include <canary/layout.hpp>                      // IWYU pragma: export
#include <canary/models/calorically_perfect.hpp>  // IWYU pragma: export
#include <canary/models/multi_temperature.hpp>    // IWYU pragma: export
#include <canary/models/thermally_perfect.hpp>    // IWYU pragma: export
#include <canary/models/two_temperature.hpp>      // IWYU pragma: export
#include <canary/nasa9.hpp>                       // IWYU pragma: export
#include <canary/quantity.hpp>                    // IWYU pragma: export
#include <canary/reactions.hpp>                   // IWYU pragma: export
#include <canary/scalar.hpp>                      // IWYU pragma: export
#include <canary/species.hpp>                     // IWYU pragma: export
#include <canary/state.hpp>                       // IWYU pragma: export
#include <canary/status.hpp>                      // IWYU pragma: export
#include <canary/thermo.hpp>                      // IWYU pragma: export
#include <canary/transport/mixture.hpp>           // IWYU pragma: export
#include <canary/units.hpp>                       // IWYU pragma: export

#endif  // CANARY_CANARY_HPP
