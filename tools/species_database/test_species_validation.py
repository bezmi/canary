"""Check species and transport conversion dimensions before C++ emission."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import collision_json
import species_json

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools/species_database"


class SpeciesValidationTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.source = self.directory / "species.json"
        self.collisions = self.directory / "collisions.json"
        self.output = self.directory / "output.hpp"
        self.data = json.loads((ROOT / "database/data/air5.json").read_text())
        self.pairs = json.loads((ROOT / "database/data/air5_collision_integrals.json").read_text())
        self.collisions.write_text(json.dumps(self.pairs))

    def write_species(self):
        self.source.write_text(json.dumps(self.data), encoding="utf-8")

    def validate(self, kind, capacity=4):
        self.write_species()
        data = species_json.load_validated(self.source)
        record = data["species"]["N2"]
        context = f"{self.source}: species N2"
        validators = {"thermo": species_json.validate_nasa9, "electronic": species_json.validate_electronic_levels,
                      "transport": species_json.validate_transport}
        validators[kind](record, context, capacity)

    def generate(self, transport=False, options=()):
        self.write_species()
        if transport:
            command = [sys.executable, str(TOOLS / "transport_to_header.py"), str(self.source), str(self.collisions),
                       "--species", "N2", "--name", "sample", "--segments", "3"]
        else:
            command = [sys.executable, str(TOOLS / "species_to_header.py"), str(self.source), "--species", "N2",
                       "--name", "sample", "--thermo", "--segments", "4", "--electronic", "--levels", "8"]
        return subprocess.run([*command, *options, "-o", str(self.output)], capture_output=True, text=True, timeout=30)

    def test_nasa9_exact_dimensions(self):
        original = copy.deepcopy(self.data)
        for case in ("short row", "long row", "boundaries", "blends", "empty", "model", "capacity"):
            with self.subTest(case=case):
                self.data = copy.deepcopy(original)
                thermo = self.data["species"]["N2"]["thermo"]
                if case == "short row":
                    thermo["coefficients"][0].pop()
                elif case == "long row":
                    thermo["coefficients"][0].append(0.0)
                elif case == "boundaries":
                    thermo["temperature_ranges"].pop()
                elif case == "blends":
                    thermo["blend_widths"].append(0.0)
                elif case == "empty":
                    thermo.update(coefficients=[], temperature_ranges=[300.0], blend_widths=[])
                elif case == "model":
                    thermo["model"] = "unknown"
                with self.assertRaises(SystemExit) as caught:
                    self.validate("thermo", capacity=1 if case == "capacity" else 4)
                self.assertIn(str(self.source), str(caught.exception))
                self.assertIn("N2", str(caught.exception))

    def test_electronic_integer_and_dimension_checks(self):
        original = copy.deepcopy(self.data)
        for case in ("length", "capacity", 1.5, True, 2147483648):
            with self.subTest(case=case):
                self.data = copy.deepcopy(original)
                levels = self.data["species"]["N2"]["electronic_levels"]
                if case == "length":
                    levels["temperatures"].pop()
                elif case != "capacity":
                    levels["degeneracies"][0] = case
                with self.assertRaises(SystemExit) as caught:
                    self.validate("electronic", capacity=1 if case == "capacity" else 8)
                self.assertIn(str(self.source), str(caught.exception))
                self.assertIn("N2", str(caught.exception))

    def test_transport_exact_dimensions_and_required_records(self):
        original = copy.deepcopy(self.data)
        for field in ("viscosity", "thermal_conductivity"):
            for case in ("short row", "long row", "boundaries", "empty", "model", "capacity", "missing", "null"):
                with self.subTest(field=field, case=case):
                    self.data = copy.deepcopy(original)
                    transport = self.data["species"]["N2"]["transport"]
                    curve = transport[field]
                    if case == "short row":
                        curve["coefficients"][0].pop()
                    elif case == "long row":
                        curve["coefficients"][0].append(0.0)
                    elif case == "boundaries":
                        curve["temperature_ranges"].pop()
                    elif case == "empty":
                        curve.update(coefficients=[], temperature_ranges=[300.0])
                    elif case == "model":
                        curve["model"] = "unknown"
                    elif case == "missing":
                        del transport[field]
                    elif case == "null":
                        transport[field] = None
                    with self.assertRaises(SystemExit) as caught:
                        self.validate("transport", capacity=1 if case == "capacity" else 3)
                    self.assertIn(str(self.source), str(caught.exception))
                    self.assertIn("N2", str(caught.exception))
        for value in (None, {}):
            self.data = copy.deepcopy(original)
            self.data["species"]["N2"]["transport"] = value
            with self.assertRaises(SystemExit):
                self.validate("transport")

    def test_metadata_units_and_selection(self):
        original = copy.deepcopy(self.data)
        for case in ("format", "units", "role", "electron charge", "electron elements"):
            with self.subTest(case=case):
                self.data = copy.deepcopy(original)
                record = self.data["species"]["N2"]
                if case == "format":
                    self.data["format"] = "old-format"
                elif case == "units":
                    self.data["units"]["molar_mass"] = "g/mol"
                elif case == "role":
                    record["role"] = "unknown"
                else:
                    record["role"] = "electron"
                    record["charge"] = -1 if case == "electron elements" else 0
                self.write_species()
                with self.assertRaises(SystemExit) as caught:
                    species_json.load_validated(self.source)
                self.assertIn(str(self.source), str(caught.exception))
        for names in (["N2", "N2"], ["missing"]):
            with self.assertRaises(SystemExit) as caught:
                species_json.validate_selection(original["species"], names, self.source)
            self.assertIn(str(self.source), str(caught.exception))
        species_json.validate_selection(original["species"], ["O2", "N2"], self.source)

    def test_collision_pair_dimensions_units_and_duplicate_keys(self):
        original = copy.deepcopy(self.pairs)
        for names in ([], ["N2"], ["N2", "N2", "N2"]):
            with self.subTest(names=names):
                data = copy.deepcopy(original)
                data["pairs"][0]["species"] = names
                with self.assertRaises(SystemExit) as caught:
                    collision_json.validate(data, self.collisions)
                self.assertIn(str(self.collisions), str(caught.exception))
        data = copy.deepcopy(original)
        data["units"]["area"] = "cm^2"
        with self.assertRaisesRegex(SystemExit, "collision units"):
            collision_json.validate(data, self.collisions)
        text = json.dumps(original).replace('"area": "m^2"', '"area": "cm^2", "area": "m^2"')
        self.assertNotEqual(text, json.dumps(original))
        self.collisions.write_text(text)
        with self.assertRaisesRegex(SystemExit, "duplicate JSON key"):
            collision_json.load(self.collisions)

    def test_cli_rejection_preserves_output(self):
        original = copy.deepcopy(self.data)
        for case in ("NASA9", "degeneracy", "CEA", "units"):
            for existing in (False, True):
                with self.subTest(case=case, existing=existing):
                    self.data = copy.deepcopy(original)
                    record = self.data["species"]["N2"]
                    if case == "NASA9":
                        record["thermo"]["coefficients"][0].pop()
                    elif case == "degeneracy":
                        record["electronic_levels"]["degeneracies"][0] = 1.5
                    elif case == "CEA":
                        record["transport"]["viscosity"]["coefficients"][0].pop()
                    else:
                        self.data["units"]["molar_mass"] = "g/mol"
                    if existing:
                        self.output.write_bytes(b"original\n")
                    elif self.output.exists():
                        self.output.unlink()
                    result = self.generate(transport=case == "CEA")
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(str(self.source), result.stderr)
                    self.assertNotIn("Traceback", result.stderr)
                    if existing:
                        self.assertEqual(self.output.read_bytes(), b"original\n")
                    else:
                        self.assertFalse(self.output.exists())

    def test_positive_capacities_required(self):
        for transport, option in ((False, "--segments"), (False, "--levels"), (True, "--segments")):
            for value in ("0", "-1"):
                with self.subTest(transport=transport, option=option, value=value):
                    result = self.generate(transport, [option, value])
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(f"{option} must be positive", result.stderr)
                    self.assertFalse(self.output.exists())

    def test_absent_null_and_empty_electronic_levels_are_equivalent(self):
        record = self.data["species"]["N2"]
        del record["electronic_levels"]
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        original = self.output.read_bytes()
        for value in (None, {"degeneracies": [], "temperatures": []}):
            record["electronic_levels"] = value
            result = self.generate()
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(self.output.read_bytes(), original)

    def test_valid_inputs_preserve_data(self):
        self.write_species()
        self.assertEqual(species_json.load_validated(self.source), self.data)
        for kind in ("thermo", "electronic", "transport"):
            self.validate(kind, capacity=8)
        self.assertEqual(collision_json.load(self.collisions), self.pairs)


if __name__ == "__main__":
    unittest.main()
