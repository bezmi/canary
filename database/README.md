# Database validation

The JSON schemas in this directory validate the structure of Canary database
files. Passing a schema check does not establish that a file can construct a
valid gas, reaction mechanism or energy-exchange model.

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
| Model construction | Numerical domains and stored-parameter representability for the selected model, scalar type and unit system. Failures report `ParameterError`. Consult each builder's documented checks and preconditions. |

Use the corresponding C++ input converter and model builder after reading a
file. Handle `DatabaseError` and `ParameterError` at the host-side loading
boundary. Check each converter's documented exceptions; a missing requested
species can also produce `DatabaseError` during input conversion.

The generated schemas do not encode every supported format, unit or model name,
every exact array length, relationships between array lengths, ordered
thermodynamic intervals, or representability after unit conversion and narrowing.
For example, NASA9 coefficient rows and collision-pair arrays currently specify
upper bounds without matching lower bounds. A schema result is therefore not a
substitute for construction. Individual schemas may contain additional
constraints, including those on elemental composition; inspect the shipped
schema when integrating an editor or an external validator.

The reaction JSON converter and reaction header generator require Python's
`jsonschema` package. They share structural validation against
`reactions.schema.json`, format/unit checks and the Arrhenius/equilibrium rules
of `mechanism_input`. The forward rate must be Arrhenius; a backward rate must
be Arrhenius or equilibrium. Arrhenius requires `A`, `n` and
`activation_temperature`; equilibrium forbids those values and
`translational_exponent`. An optional `null` value is unset, as in the C++ file
reader. `other_temperature` remains optional for either model. Errors identify
the source and the reaction or structural field before output is written.
The header generator also rejects duplicate JSON keys. Numerical domains and
species-dependent checks still belong to model construction.

The energy-exchange JSON converter and header generator use the same schema
validation helper and require `jsonschema`. Their shared model checks match
`energy_exchange_input` and `electron_exchange_input`: Park 1993 requires
`reference_temperature`, Kim 2021 requires `exponent`, and each correction
forbids the other models' optional fields. Neutral electron cross-sections
require `a`, `b`, `c` and `maximum_temperature`. Optional nulls are unset.
An absent or null correction means no correction; an explicit correction record
must name a supported physical model. The shipped schema rejects fractional
reaction indices before C++ generation. The optional species file supplying
ionisation energies is checked for structure, format, units and duplicate keys.
Mechanism-relative index bounds, charge conservation and numerical domains
still require model-construction checks.

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
