# Contributing to Canary

Contributions are welcome: corrections, tests, data for further species and
mechanisms, and new models. This file states the terms under which
contributions are accepted and what a change needs before it can be merged.

## Licence of contributions

Canary is licensed under the GNU Affero General Public License, version 3 or
(at your option) any later version, along with the Canary Cerberus Exception
(see [LICENSE](LICENSE) and
[LICENSE.CERBERUS-EXCEPTION](LICENSE.CERBERUS-EXCEPTION)).

By submitting a contribution to Canary you agree to the following terms.

1. **You retain the copyright** in your contribution.

2. **Project licence.** You license your contribution to every recipient of
   Canary under the GNU Affero General Public License, version 3 or (at your
   option) any later version, along with the Canary Cerberus Exception. You
   do not remove the Exception from your contribution.

3. **Licence to the maintainer.** You grant Shahzeb Imran (the "Maintainer")
   a perpetual, worldwide, non-exclusive, royalty-free, irrevocable licence
   to reproduce, modify, publicly display, publicly perform, sublicense and
   distribute your contribution and works based on it, including the right
   to license them under terms other than those of section 2. The Maintainer
   will keep making your contribution available under the terms of section 2
   in the versions of Canary that are published.

4. **Patents.** You grant the Maintainer and every recipient of Canary a
   perpetual, worldwide, non-exclusive, royalty-free, irrevocable licence
   under the patent claims that you own or control and that your
   contribution, alone or in combination with Canary, necessarily infringes,
   to make, use, sell, offer for sale, import and otherwise transfer the
   contribution and Canary.

5. **Your right to contribute.** You declare that the contribution is your
   original work, or that you have the right to submit it under these terms.
   If your employer has rights in work that you create, you declare that the
   employer has permitted the contribution or has waived those rights. You
   identify any part of the contribution that is not your own work, with its
   source and its licence.

6. **No warranty.** You provide the contribution as it is, without warranty
   of any kind.

State your agreement in the description of your first pull request with the
sentence "I agree to the terms of CONTRIBUTING.md for this and my later
contributions to Canary."

## Data and references

Every value in a data file under `database/data/` must come from a
publication, and its `reference` field must cite that publication with enough
detail to find the value (table, equation or page). Values copied from the
data files of other software are not accepted, because their licence and
their origin cannot be passed on. Tests may compare Canary with other
software.

## Before you open a pull request

1. Build and run the tests with both compilers:

   ```sh
   cmake --workflow --preset debug
   cmake --workflow --preset debug-clang
   ```

   If you have a CUDA toolkit and device, also run
   `cmake --workflow --preset cuda-runtime` for changes to per-cell code.

2. Add or extend tests for the changed behaviour. A change to an existing
   expected value or tolerance needs a stated reason in the commit message.

3. Format the changed lines with the repository's `.clang-format`:

   ```sh
   git clang-format
   ```

4. After changing a file under `database/data/`, regenerate the headers in
   `include/canary/data/` with the tools in `tools/species_database/` and
   commit both.

## Conventions

- Per-cell functions are annotated `CANARY_HOST_DEVICE`, do not allocate and
  do not throw. Host builders validate their input and throw `ParameterError`.
- Physical values use the typed quantities of `canary/quantity.hpp`.
- Prefer templates, concepts and `std::variant` to virtual functions; nothing
  in a per-cell function is dispatched at run time.
- Comments state what the code does, its units and its limits, and stand on
  their own. Commit messages state what changed and why, in plain technical
  language.
