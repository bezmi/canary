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
