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
            "format": "canary-reactions-2",
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
        for backward in ({"model": "equilibrium"}, {"model": "equilibrium", "temperature": {"kind": "free_electron"}},
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
        record["forward"].update(translational_exponent=None, other_temperature=None, temperature=None)
        record["backward"].update(A=None, n=None, activation_temperature=None,
                                  translational_exponent=None, other_temperature=None, temperature=None)
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

    def test_physical_temperature_requirements(self):
        valid = [{"kind": k} for k in ("translation", "rotation", "free_electron")]
        valid += [{"kind": k, "species": "N2"} for k in ("vibration", "bound_electronic")]
        for role, field in (("forward", "other_temperature"), ("backward", "temperature")):
            for requirement in valid:
                data = copy.deepcopy(self.data)
                data["reactions"][0][role][field] = requirement
                reaction_json.validate(data, self.source)
            for requirement, reason in (
                ({"kind": "bad"}, "unknown temperature kind"),
                ({"kind": "vibration"}, "needs a species"),
                ({"kind": "bound_electronic", "species": ""}, "needs a species"),
                ({"kind": "translation", "species": "N2"}, "does not take a species"),
                ({"kind": "rotation", "species": "N2"}, "does not take a species"),
                ({"kind": "free_electron", "species": "e-"}, "does not take a species")):
                data = copy.deepcopy(self.data)
                data["reactions"][0][role][field] = requirement
                self.assert_rejected(data, reason)
        for role, field in (("forward", "temperature"), ("backward", "other_temperature")):
            data = copy.deepcopy(self.data)
            data["reactions"][0][role][field] = {"kind": "translation"}
            self.assert_rejected(data, f"forbids {field}")

    def test_old_format_and_string_requirements_are_rejected(self):
        self.data["format"] = "canary-reactions-1"
        with self.assertRaises(SystemExit):
            reaction_json.validate(self.data, self.source)
        self.data["format"] = "canary-reactions-2"
        self.data["reactions"][0]["forward"]["other_temperature"] = "electron"
        self.source.write_text(json.dumps(self.data), encoding="utf-8")
        self.output.write_bytes(b"original output")
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.output.read_bytes(), b"original output")

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

    def test_arrhenius_input_domains(self):
        invalid = {"A": (0.0, -1.0, float("nan"), float("inf")),
                   "n": (float("nan"), float("inf")),
                   "activation_temperature": (-1.0, float("nan"), float("inf")),
                   "translational_exponent": (-0.01, 1.01, float("nan"), float("inf"))}
        for role in ("forward", "backward"):
            for field, values in invalid.items():
                for value in values:
                    with self.subTest(role=role, field=field, value=value):
                        data = copy.deepcopy(self.data)
                        rate = copy.deepcopy(data["reactions"][0]["forward"])
                        rate[field] = value
                        data["reactions"][0][role] = rate
                        with self.assertRaises(SystemExit) as caught:
                            reaction_json.validate(data, self.source)
                        message = str(caught.exception)
                        self.assertIn(str(self.source), message)
                        self.assertIn(role, message)
                        self.assertIn(field, message)
        for role in ("forward", "backward"):
            for exponent in (None, 0.0, 1.0):
                data = copy.deepcopy(self.data)
                data["reactions"][0][role] = {"model": "arrhenius", "A": 1e-300, "n": -2.0,
                                              "activation_temperature": -0.0,
                                              "translational_exponent": exponent}
                reaction_json.validate(data, self.source)

    def test_stoichiometry_and_efficiency_domains(self):
        for role in ("reactants", "products"):
            for value in (0.0, -1.0, 0.5, float("nan"), float("inf")):
                with self.subTest(role=role, value=value):
                    data = copy.deepcopy(self.data)
                    data["reactions"][0][role] = {"N2": value}
                    with self.assertRaises(SystemExit) as caught:
                        reaction_json.validate(data, self.source)
                    self.assertIn(str(self.source), str(caught.exception))
                    self.assertIn(role, str(caught.exception))
            data = copy.deepcopy(self.data)
            data["reactions"][0][role] = {}
            self.assert_rejected(data, f"{role} must contain at least one species")
        for value in (-1.0, float("nan"), float("inf")):
            with self.subTest(efficiency=value):
                data = copy.deepcopy(self.data)
                data["reactions"][0]["third_body"] = {"efficiencies": {"N2": value}}
                with self.assertRaises(SystemExit) as caught:
                    reaction_json.validate(data, self.source)
                self.assertIn(str(self.source), str(caught.exception))
                self.assertIn("N2", str(caught.exception))
        for coefficients in ({"N2": 1}, {"N2": 1.0, "N": 2.0}):
            for efficiencies in ({}, {"N2": -0.0}, {"N2": 1.0}):
                data = copy.deepcopy(self.data)
                data["reactions"][0]["reactants"] = coefficients
                data["reactions"][0]["third_body"] = {"efficiencies": efficiencies}
                reaction_json.validate(data, self.source)

    def test_temperature_limits_distinct_as_double_inputs(self):
        self.data["temperature_limits"] = {"lower": 2**53, "upper": 2**53 + 1}
        self.source.write_text(json.dumps(self.data))
        self.output.write_bytes(b"original\n")
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(str(self.source), result.stderr)
        self.assertIn("temperature limits", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.assertEqual(self.output.read_bytes(), b"original\n")

    def test_temperature_limits_and_empty_mechanisms(self):
        for lower, upper in ((0.0, 300.0), (-1.0, 300.0), (300.0, 300.0), (301.0, 300.0),
                             (float("nan"), 300.0), (300.0, float("nan")), (300.0, float("inf"))):
            with self.subTest(lower=lower, upper=upper):
                data = copy.deepcopy(self.data)
                data["temperature_limits"] = {"lower": lower, "upper": upper}
                with self.assertRaises(SystemExit) as caught:
                    reaction_json.validate(data, self.source)
                self.assertIn(str(self.source), str(caught.exception))
                self.assertIn("temperature", str(caught.exception))
        self.data["reactions"] = []
        with self.assertRaisesRegex(SystemExit, "at least one reaction"):
            reaction_json.validate(self.data, self.source)

    def test_invalid_domains_preserve_header_output(self):
        original = copy.deepcopy(self.data)
        for case in ("coefficient", "efficiency", "rate", "limits", "empty"):
            for existing in (False, True):
                with self.subTest(case=case, existing=existing):
                    self.data = copy.deepcopy(original)
                    if case == "coefficient":
                        self.data["reactions"][0]["reactants"]["N2"] = 0.5
                    elif case == "efficiency":
                        self.data["reactions"][0]["third_body"] = {"efficiencies": {"N2": -1.0}}
                    elif case == "rate":
                        self.data["reactions"][0]["forward"]["A"] = -1.0
                    elif case == "limits":
                        self.data["temperature_limits"] = {"lower": 50000.0, "upper": 300.0}
                    else:
                        self.data["reactions"] = []
                    self.source.write_text(json.dumps(self.data))
                    if existing:
                        self.output.write_bytes(b"original\n")
                    elif self.output.exists():
                        self.output.unlink()
                    result = self.generate()
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(str(self.source), result.stderr)
                    self.assertNotIn("Traceback", result.stderr)
                    if existing:
                        self.assertEqual(self.output.read_bytes(), b"original\n")
                    else:
                        self.assertFalse(self.output.exists())

    def test_all_shipped_mechanisms_pass(self):
        for name in ("air5_kim_jo", "air5_park", "air11_kim_jo"):
            with self.subTest(name=name):
                reaction_json.load(ROOT / "database" / "data" / f"{name}.json")


if __name__ == "__main__":
    unittest.main()
