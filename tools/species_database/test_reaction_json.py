"""Check reaction model validation before JSON and header output."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import reaction_json

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools" / "species_database"


class ReactionJsonTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.source = self.directory / "reactions.json"
        self.output = self.directory / "reactions.hpp"
        self.data = {
            "format": "canary-reactions-1",
            "units": {"quantity": "mol", "length": "m", "time": "s", "temperature": "K"},
            "temperature_limits": {"lower": 300.0, "upper": 50000.0},
            "reactions": [{"equation": "N2 <=> N + N", "reactants": {"N2": 1.0}, "products": {"N": 2.0},
                           "forward": {"model": "arrhenius", "A": 1.25e9, "n": 0.0,
                                       "activation_temperature": 0.0},
                           "backward": {"model": "equilibrium"}}],
        }

    def assert_rejected(self, data, reason):
        with self.assertRaises(SystemExit) as caught:
            reaction_json.validate(data, self.source)
        message = str(caught.exception)
        self.assertIn(str(self.source), message)
        self.assertIn("reaction 0 (N2 <=> N + N)", message)
        self.assertIn(reason, message)

    def generate(self):
        return subprocess.run([sys.executable, str(TOOLS / "reactions_to_header.py"), str(self.source),
                               "--name", "sample", "--alias", "Sample",
                               "--species-header", "canary/data/air5.hpp", "--species-array", "air5_species",
                               "-o", str(self.output)], capture_output=True, text=True, timeout=30)

    def test_supported_rates_and_zero_parameters(self):
        for backward in ({"model": "equilibrium"}, {"model": "equilibrium", "other_temperature": "electron"},
                         copy.deepcopy(self.data["reactions"][0]["forward"])):
            with self.subTest(backward=backward):
                self.data["reactions"][0]["backward"] = backward
                original = copy.deepcopy(self.data)
                reaction_json.validate(self.data, self.source)
                self.assertEqual(self.data, original)

    def test_missing_and_null_arrhenius_fields(self):
        for role in ("forward", "backward"):
            for field in ("A", "n", "activation_temperature"):
                for null in (False, True):
                    with self.subTest(role=role, field=field, null=null):
                        data = copy.deepcopy(self.data)
                        rate = copy.deepcopy(data["reactions"][0]["forward"])
                        if null:
                            rate[field] = None
                        else:
                            del rate[field]
                        data["reactions"][0][role] = rate
                        self.assert_rejected(data, f"{role} Arrhenius rate requires {field}")

    def test_equilibrium_forbidden_fields_include_zero(self):
        for field in ("A", "n", "activation_temperature", "translational_exponent"):
            for value in (0.0, 1.0):
                with self.subTest(field=field, value=value):
                    data = copy.deepcopy(self.data)
                    data["reactions"][0]["backward"][field] = value
                    self.assert_rejected(data, f"equilibrium backward rate forbids {field}")

    def test_unknown_models_are_not_reinterpreted(self):
        for role in ("forward", "backward"):
            for model in ("power", "", "Arrhenius"):
                with self.subTest(role=role, model=model):
                    data = copy.deepcopy(self.data)
                    data["reactions"][0][role] = dict(data["reactions"][0]["forward"], model=model)
                    self.assert_rejected(data, f"{role} rate model must be")
        self.data["reactions"][0]["forward"] = {"model": "equilibrium"}
        self.assert_rejected(self.data, "forward rate model must be arrhenius")

    def test_structural_format_and_unit_errors(self):
        cases = []
        data = copy.deepcopy(self.data)
        data["reactions"][0]["forward"]["A"] = True
        cases.append((data, "reactions/0/forward/A"))
        data = copy.deepcopy(self.data)
        data["reactions"][0]["forward"]["unexpected"] = 1
        cases.append((data, "unexpected"))
        data = copy.deepcopy(self.data)
        del data["reactions"][0]["backward"]["model"]
        cases.append((data, "model"))
        data = copy.deepcopy(self.data)
        data["format"] = "other"
        cases.append((data, "format must be"))
        data = copy.deepcopy(self.data)
        data["units"]["quantity"] = "kmol"
        cases.append((data, "units must be"))
        for data, reason in cases:
            with self.subTest(reason=reason):
                with self.assertRaises(SystemExit) as caught:
                    reaction_json.validate(data, self.source)
                self.assertIn(str(self.source), str(caught.exception))
                self.assertIn(reason, str(caught.exception))

    def test_optional_nulls_match_omitted_fields_in_header(self):
        self.source.write_text(json.dumps(self.data), encoding="utf-8")
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        original = self.output.read_bytes()
        record = self.data["reactions"][0]
        record["forward"].update(translational_exponent=None, other_temperature=None)
        record["backward"].update(A=None, n=None, activation_temperature=None,
                                  translational_exponent=None, other_temperature=None)
        record["third_body"] = None
        self.source.write_text(json.dumps(self.data), encoding="utf-8")
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_bytes(), original)

    def test_header_rejection_preserves_output(self):
        for backward in ({"model": "power", "A": 1.0, "n": 0.0, "activation_temperature": 0.0},
                         {"model": "equilibrium", "n": 0.0},
                         {"model": "arrhenius", "A": 1.0, "n": 0.0}):
            for existing in (False, True):
                with self.subTest(backward=backward, existing=existing):
                    if existing:
                        self.output.write_bytes(b"original output\n")
                    elif self.output.exists():
                        self.output.unlink()
                    self.data["reactions"][0]["backward"] = backward
                    self.source.write_text(json.dumps(self.data), encoding="utf-8")
                    result = self.generate()
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(str(self.source), result.stderr)
                    self.assertIn("reaction 0 (N2 <=> N + N)", result.stderr)
                    self.assertNotIn("Traceback", result.stderr)
                    if existing:
                        self.assertEqual(self.output.read_bytes(), b"original output\n")
                    else:
                        self.assertFalse(self.output.exists())

    def test_duplicate_json_keys_are_rejected(self):
        text = json.dumps(self.data).replace('"model": "equilibrium"',
                                            '"model": "power", "model": "equilibrium"')
        self.assertNotEqual(text, json.dumps(self.data))
        self.source.write_text(text, encoding="utf-8")
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("duplicate JSON key", result.stderr)
        self.assertIn(str(self.source), result.stderr)
        self.assertFalse(self.output.exists())

    def test_all_shipped_mechanisms_pass(self):
        for name in ("air5_kim_jo", "air5_park", "air11_kim_jo"):
            with self.subTest(name=name):
                reaction_json.load(ROOT / "database" / "data" / f"{name}.json")


if __name__ == "__main__":
    unittest.main()
