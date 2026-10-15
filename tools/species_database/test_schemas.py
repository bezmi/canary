"""Test structural rejection independently of C++ database validation.

Usage: test_schemas.py <canary source directory>
Requires jsonschema; no converters or C++ executables are called.
"""

import copy
import json
import sys
import unittest
from pathlib import Path

import jsonschema

ROOT = Path(sys.argv.pop(1))
FILES = {
    "species": ("air5", "air11"),
    "reactions": ("air5_kim_jo", "air5_park", "air11_kim_jo"),
    "energy_exchange": ("air5_park_exchange", "air11_kim_jo_exchange"),
    "collision_integrals": ("air5_collision_integrals",),
}


REQUIRED_FIELDS = {
    "species": ("format", "units", "species"),
    "reactions": ("format", "units", "temperature_limits", "reactions"),
    "energy_exchange": ("format", "units", "vibration_translation", "chemistry_vibration",
                        "electron_translation", "vibration_electron", "electron_impact_ionisation"),
    "collision_integrals": ("format", "units", "pairs"),
}


class SchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.validators = {}
        cls.documents = {}
        for kind, files in FILES.items():
            schema = json.loads((ROOT / "database" / f"{kind}.schema.json").read_text())
            jsonschema.Draft202012Validator.check_schema(schema)
            cls.validators[kind] = jsonschema.Draft202012Validator(schema)
            for name in files:
                cls.documents[name] = json.loads(
                    (ROOT / "database/data" / f"{name}.json").read_text()
                )

    def reject(self, kind, path, value, keyword):
        document = copy.deepcopy(self.documents[FILES[kind][0]])
        parent = document
        for key in path[:-1]:
            parent = parent[key]
        # Existing fields must be present so misspelled fixtures cannot silently
        # test an unrelated additional-property error.
        self.assertIn(path[-1], parent)
        parent[path[-1]] = value
        errors = list(self.validators[kind].iter_errors(document))
        self.assertTrue(
            any(error.validator == keyword and tuple(error.absolute_path) == path
                for error in errors),
            [error.message for error in errors],
        )

    def test_shipped_files(self):
        for kind, names in FILES.items():
            for name in names:
                with self.subTest(file=name):
                    self.assertEqual(list(self.validators[kind].iter_errors(self.documents[name])), [])

    def test_required_fields_and_unknown_fields(self):
        for kind in FILES:
            base = self.documents[FILES[kind][0]]
            for key in REQUIRED_FIELDS[kind]:
                with self.subTest(schema=kind, missing=key):
                    document = copy.deepcopy(base)
                    del document[key]
                    errors = list(self.validators[kind].iter_errors(document))
                    self.assertTrue(any(error.validator == "required" for error in errors))
            with self.subTest(schema=kind, unknown="unexpected"):
                document = copy.deepcopy(base)
                document["unexpected"] = 1
                errors = list(self.validators[kind].iter_errors(document))
                self.assertTrue(any(error.validator == "additionalProperties" for error in errors))

    def test_nested_field_types(self):
        cases = (
            ("species", ("species", "N2", "molar_mass"), "0.028", "type"),
            ("reactions", ("reactions", 0, "forward", "A"), "1e12", "anyOf"),
            ("energy_exchange", ("vibration_translation", 0, "relaxation_time", "a"), [], "type"),
            ("collision_integrals", ("pairs", 0, "pi_omega_11", "A"), "0", "type"),
        )
        for kind, path, value, keyword in cases:
            with self.subTest(schema=kind, path=path):
                self.reject(kind, path, value, keyword)

    def test_fixed_array_upper_bounds(self):
        cases = (
            ("species", ("species", "N2", "thermo", "coefficients"), [[0.0] * 10]),
            ("collision_integrals", ("pairs", 0, "species"), ["N2", "O2", "NO"]),
        )
        for kind, path, value in cases:
            with self.subTest(schema=kind, path=path):
                # NASA rows are nested arrays, while a collision pair is one array.
                document = copy.deepcopy(self.documents[FILES[kind][0]])
                parent = document
                for key in path[:-1]:
                    parent = parent[key]
                parent[path[-1]] = value
                errors = list(self.validators[kind].iter_errors(document))
                expected_path = path + (0,) if kind == "species" else path
                self.assertTrue(any(error.validator == "maxItems" and
                                    tuple(error.absolute_path) == expected_path for error in errors))


if __name__ == "__main__":
    unittest.main()
