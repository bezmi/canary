# Database validation

The JSON schemas in this directory validate the structure of Canary database
files. Passing a schema check does not establish that a file can construct a
valid gas, reaction mechanism or energy-exchange model.

Supported source importers are `tools/species_database/thermo_inp_to_json.py`
for NASA Glenn records and `tools/species_database/cantera_to_json.py` for
Cantera NASA9 YAML. GDTk Lua import is not supported. Existing GDTk-derived
Canary datasets retain their data and provenance; the JSON-to-header tools
consume those datasets directly.

## Check a file

With Python and `jsonschema` installed, run from the repository root:

```sh
python3 tools/species_database/check_species_files.py \
  database/species.schema.json database/data/air5.json
```

The command reports each file and its schema errors, and returns nonzero if
validation fails. Supply one schema and at least one file per group; separate
multiple groups with `--`. The loader also rejects duplicate JSON object keys.
A successful result establishes the structural checks described below only.

## Validation responsibilities

| Operation | Checks |
|---|---|
| JSON schema | Required properties, unknown properties, JSON value types and the array bounds and metadata constraints present in that schema. |
| C++ `parse_*` / `read_*_file` | Decoding into the database record types, supported file format and units; species parsing also checks roles and elemental metadata. Failures report `DatabaseError`. |
| Input conversion | Supported model identifiers, model-dependent required/forbidden fields, requested record counts and related array sizes. For example, `nasa9_data` checks the relationship between intervals, boundaries and blend widths; `mechanism_input` checks Arrhenius/equilibrium fields. |
| Python header generation | Shipped schemas, format/units, supported models, required fields, exact record dimensions, integer-valued fields and basic SI parameter domains before output. The generated records remain host inputs, not validated model storage. |
| Model construction | Numerical domains and stored-parameter representability for the selected model, scalar type and unit system. Failures report `ParameterError`. Consult each builder's documented checks and preconditions. |

Use the corresponding C++ input converter and model builder after reading a
file. Handle `DatabaseError` and `ParameterError` at the host-side loading
boundary. Check each converter's documented exceptions; a missing requested
species can also produce `DatabaseError` during input conversion.

Header generation does not select the caller's precision, unit scales or model
configuration. Construct the generated inputs with the corresponding C++ builder
before using a numerical kernel. C++ construction owns coupled physical rules
(such as heat-capacity positivity, oscillator ordering and cross-section interval
validity), integer storage/order bounds, species/mechanism relationships and
representability after conversion. Numerical checks occur during construction;
kernels consume the resulting strongly typed data without repeating them.

The generated schemas do not encode every supported format, unit or model name,
every exact array length, relationships between array lengths, ordered
thermodynamic intervals, or representability after unit conversion and narrowing.
For example, NASA9 coefficient rows and collision-pair arrays currently specify
upper bounds without matching lower bounds. A schema result is therefore not a
substitute for construction. Individual schemas may contain additional
constraints, including those on elemental composition; inspect the shipped
schema when integrating an editor or an external validator.

The reaction header generator requires Python's `jsonschema` package. It uses
structural validation against
`reactions.schema.json`, format/unit checks and the Arrhenius/equilibrium rules
of `mechanism_input`. The forward rate must be Arrhenius; a backward rate must
be Arrhenius or equilibrium. Arrhenius requires `A`, `n` and
`activation_temperature`; equilibrium forbids those values and
`translational_exponent`. An optional `null` value is unset, as in the C++ file
reader. Reaction files use `canary-reactions-2`; version 1 is rejected.
An Arrhenius `other_temperature` is an optional physical requirement, for example
`{"kind": "vibration", "species": "N2"}`. Equilibrium uses its separate
`temperature` field, for example `{"kind": "translation"}`, and forbids
`other_temperature`. Arrhenius forbids `temperature`. Supported kinds are
`translation`, `rotation`, `vibration`, `bound_electronic`, and `free_electron`;
only vibration and bound electronic require a species name. A model-aware
builder requires an equilibrium temperature and an Arrhenius other temperature
when its translational exponent is not one. Explicit single-temperature
construction permits omission and collapses valid requirements. Errors identify
the source and the reaction or structural field before output is written.
The header generator also rejects duplicate JSON keys. It requires a
nonempty reaction list, nonempty sides with positive integer coefficients,
finite nonnegative third-body efficiencies, positive finite Arrhenius `A`,
finite `n`, nonnegative finite activation temperatures, and translational
exponents in [0, 1]. File temperature limits must be finite and satisfy
`0 < lower < upper`. Integer storage and reaction-order bounds, representability
after unit conversion, species identity and conservation still require model
construction checks.

The energy-exchange header generator uses the schema validation helper and
requires `jsonschema`. Its model checks match
`energy_exchange_input` and `electron_exchange_input`: Park 1993 requires
`reference_temperature`, Kim 2021 requires `exponent`, and each correction
forbids the other models' optional fields. Neutral electron cross-sections
require `a`, `b`, `c` and `maximum_temperature`. Optional nulls are unset.
An absent or null correction means no correction; an explicit correction record
must name a supported physical model. The shipped schema rejects fractional
reaction indices before C++ generation; shared validation also rejects negative
indices. The optional species file supplying ionisation energies is checked for
structure, format, units and duplicate keys. Each referenced ionisation energy
must be present, finite and positive before header output.
Millikan-White coefficients and Kim exponents must be finite. Correction
cross-sections, Park reference temperatures, dissociation energies and
vibrational temperatures must be finite and positive. Neutral cross-section
coefficients must be finite and their maximum temperature positive and finite.
The constant, the exponent and the three inverse coefficients of a
`power_exponential` relaxation time must be finite, and its reference
temperature finite and positive. Mechanism-relative index
bounds, charge conservation, oscillator ordering, cross-section validity over
its interval and representability after conversion still require
model-construction checks.

Species and transport header generators also validate the shipped schemas,
file formats, SI units and selected species before output. They require exact
NASA9/CEA coefficient row widths, matching temperature-boundary and blend counts,
matching electronic degeneracy/temperature counts, and sufficient positive
array capacities. Schema checks reject fractional or out-of-int32 electronic
degeneracies before integer emission. Missing or null electronic levels still
produce zero active levels. Collision records require exactly two species names.
Species thermodynamic inputs
require positive finite molar mass and finite charge. NASA9 and CEA boundaries
must be finite, positive and strictly increasing, with finite coefficients.
NASA9 blend widths must be finite and positive; blend regions must lie within
the fit without overlap (touching regions are allowed). Electronic degeneracies
must be positive, level temperatures finite and nonnegative, and the ground
level temperature zero; no ordering of excited levels is required. Lennard-Jones
sigma and epsilon must be finite and positive. Collision-fit coefficients must
be finite. These input checks do not establish temperature-dependent physical
validity or representability after conversion, and do not replace numerical
builders.

## Tests and schema maintenance

With a configured build that enables the database and Python validation:

```sh
ctest --test-dir build/debug -R '^canary.database.schema' --output-on-failure
python3 tools/species_database/test_schemas.py .
```

`canary.database.schema` validates the shipped data. The separate
`canary.database.schema_rejection` test changes copies of those files to check
missing fields, unknown properties, wrong nested types and oversized fixed
arrays. It does not call the C++ parser or numerical builders. C++ parsing,
input conversion and model-construction checks have their own tests.

The schemas are generated from the C++ database record definitions. The
`canary.database` C++ test compares generated and shipped schema text. Update
the record definitions and generation together; do not maintain an independent
set of physical rules by hand in the generated JSON files.

The `canary.database.generated` CTest target exercises the generation-to-construction
boundary. The build runs all four header generators, then compiles and constructs
the generated species, reactions, exchange and transport inputs with float/double
and SI/scaled units. Source-valid fixtures verify that C++ construction still
rejects negative heat capacity and a mass unrepresentable in the selected storage
precision. These tests require the same Python dependencies as the converter tests.

## File replacement

The C++ file writers and Python JSON converters serialize before creating a
temporary file beside the requested destination. They check writing, flushing
and closing before renaming that file over the destination. An existing destination must be a regular file;
symlinks, including dangling symlinks, are rejected. Existing permission bits
are copied to the replacement. New files use the process's normal file-creation
permissions. Parent directories must already exist.

A failure before replacement leaves the destination untouched and attempts to
remove the temporary file. C++ errors report `DatabaseError`; Python filesystem
errors report `OSError`. Both include the destination path and failed operation.
Python JSON serialization rejects non-finite numbers with `ValueError` and
unsupported objects with `TypeError`, before touching any output file. C++ reads
report I/O failures instead of returning partial contents. Hard links to an old
destination retain its previous contents because replacement changes the
destination's directory entry.

Atomic replacement prevents readers from seeing a partially written document.
It does not provide power-loss durability: the writers do not synchronize files
or directories to persistent storage. Callers must coordinate concurrent changes
to the destination and its parent directories; inspection and replacement do
not constitute a filesystem transaction.

## Reaction construction from species definitions

Reaction construction takes one ordered set of host species definitions for
names, masses, charge numbers and elemental composition:

```cpp
const auto species = database::nasa9_species<4>(species_file, names);
const auto input = database::mechanism_input<13, 5>(reaction_file);
const auto mechanism = make_reaction_mechanism(
    input, reaction_species(species), scales);
```

`reaction_species` borrows fixed-array definitions; keep the array and its
borrowed names/composition alive through construction. For `OwnedInput`, it
retains an owning copy or accepts a moved owner, including a temporary returned
by the database converter. The numerical mechanism stores no host views.
The supplied species order must match the thermodynamic model used to evaluate
rates. Bundled `make_air5_park`, `make_air5_kim_jo` and `make_air11_kim_jo`
builders use their bundled species definitions and take `scales` first;
they no longer accept separately constructed `SpeciesData`.

Every explicitly listed reactant/product requires known elemental composition,
including species present unchanged on both sides. Each participant's charge
number must be an integer in the signed 32-bit range. `ChargeNumber` represents
that host chemistry value; it is distinct from `Charge<T,U>` in physical units.
Element counts and charge numbers are balanced with checked signed 64-bit sums.
A nonzero balance or arithmetic overflow raises `ParameterError` with reaction
context. General thermodynamic inputs may still contain fractional charges;
mixture charge neutrality is not a reaction-construction requirement.
