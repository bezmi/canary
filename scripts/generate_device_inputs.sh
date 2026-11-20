#!/bin/sh
# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
# Regenerates tests/device/derived_inputs.hpp, the inputs of the device-check
# fixtures derived by forward closures (see tests/device/derive_inputs.cpp).
#
# Usage: scripts/generate_device_inputs.sh <build tree>
# The build tree must be configured with CANARY_BUILD_CUDA_CHECKS or
# CANARY_BUILD_HIP_CHECKS (the generator itself needs only the host compiler).
# The Canary commit, and whether the working tree had uncommitted changes, is
# recorded in the generated header. Review the diff before committing it.
set -eu
tree=$1
root=$(git rev-parse --show-toplevel)
cmake --build "$tree" --target canary_device_derive_inputs
commit=$(git -C "$root" rev-parse HEAD)
if git -C "$root" diff --quiet HEAD; then
  state="no uncommitted changes"
else
  state="uncommitted changes"
fi
output="$root/tests/device/derived_inputs.hpp"
"$tree/tests/device/canary_device_derive_inputs" "Canary $commit ($state)" > "$output.new"
mv "$output.new" "$output"
clang-format -i "$output"
