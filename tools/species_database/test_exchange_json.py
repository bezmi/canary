"""Check exchange input validation before JSON and C++ header output."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import exchange_json

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools" / "species_database"
LISTS = ("vibration_translation", "chemistry_vibration", "electron_translation",
         "vibration_electron", "electron_impact_ionisation")


class ExchangeJsonTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.source = self.directory / "exchange.json"
        self.output = self.directory / "exchange.hpp"
        self.species_source = self.directory / "species.json"
        self.species_source.write_bytes((ROOT / "database/data/air11.json").read_bytes())
        self.data = json.loads((ROOT / "database/data/air11_kim_jo_exchange.json").read_text())
        for name in LISTS:
            self.data[name] = self.data[name][:1]
            self.assertTrue(self.data[name], name)
        self.assertEqual(self.correction(self.data)["model"], "kim_2021")
        self.assertEqual(self.data["electron_translation"][0]["cross_section"]["model"], "gnoffo_neutral")

    @staticmethod
    def correction(data):
        return data["vibration_translation"][0]["relaxation_time"]["high_temperature_correction"]

    def assert_rejected(self, data, context, reason):
        with self.assertRaises(SystemExit) as caught:
            exchange_json.validate(data, self.source)
        message = str(caught.exception)
        self.assertIn(str(self.source), message)
        self.assertIn(context, message)
        self.assertIn(reason, message)

    def generate(self):
        self.source.write_text(json.dumps(self.data), encoding="utf-8")
        return subprocess.run([sys.executable, str(TOOLS / "exchange_to_header.py"), str(self.source),
                               "--name", "sample", "--alias", "Sample",
                               "--species-header", "canary/data/air11.hpp", "--species-array", "air11_species",
                               "--mechanism-header", "canary/data/air11_kim_jo.hpp", "--mechanism-alias", "Air11KimJo",
                               "--model-header", "canary/data/air11.hpp", "--model-alias", "Air11ThreeTemperature",
                               "--species-file", str(self.species_source), "-o", str(self.output)],
                              capture_output=True, text=True, timeout=30)

    def test_correction_required_and_forbidden_fields(self):
        for model in ("park_1985", "park_1993", "kim_2021"):
            base = {"model": model, "cross_section": 1e-20}
            if model == "park_1993":
                base["reference_temperature"] = 50000.0
            elif model == "kim_2021":
                base["exponent"] = 0.0
            data = copy.deepcopy(self.data)
            data["vibration_translation"][0]["relaxation_time"]["high_temperature_correction"] = base
            exchange_json.validate(data, self.source)
            for field, required_model in (("reference_temperature", "park_1993"), ("exponent", "kim_2021")):
                for value in (None, 0.0):
                    with self.subTest(model=model, field=field, value=value):
                        changed = copy.deepcopy(data)
                        self.correction(changed)[field] = value
                        required = model == required_model
                        if (value is not None) == required:
                            exchange_json.validate(changed, self.source)
                        else:
                            self.assert_rejected(changed, "vibration_translation 0", f"{'requires' if required else 'forbids'} {field}")
                if model == required_model:
                    changed = copy.deepcopy(data)
                    del self.correction(changed)[field]
                    self.assert_rejected(changed, "vibration_translation 0", f"requires {field}")

    def test_neutral_cross_section_fields_are_required_even_when_zero(self):
        for field in ("a", "b", "c", "maximum_temperature"):
            for state in ("missing", "null", "zero"):
                with self.subTest(field=field, state=state):
                    data = copy.deepcopy(self.data)
                    section = data["electron_translation"][0]["cross_section"]
                    if state == "missing":
                        del section[field]
                    else:
                        section[field] = None if state == "null" else 0.0
                    if state == "zero":
                        exchange_json.validate(data, self.source)
                    else:
                        self.assert_rejected(data, "electron_translation 0", f"requires {field}")

    def test_unsupported_models(self):
        paths = [("vibration_translation",), ("vibration_translation", "relaxation_time"),
                 ("vibration_translation", "relaxation_time", "high_temperature_correction"),
                 ("chemistry_vibration",), ("electron_translation",),
                 ("electron_translation", "cross_section"), ("vibration_electron",),
                 ("vibration_electron", "relaxation_time")]
        for path in paths:
            with self.subTest(path=path):
                data = copy.deepcopy(self.data)
                record = data[path[0]][0]
                for key in path[1:]:
                    record = record[key]
                record["model"] = "unsupported"
                self.assert_rejected(data, f"{path[0]} 0", "model must be")
        self.correction(self.data)["model"] = "none"
        self.assert_rejected(self.data, "vibration_translation 0", "model must be")

    def test_reaction_indices_are_not_truncated(self):
        for name in ("chemistry_vibration", "electron_impact_ionisation"):
            for value in (0.5, True, "0", 2147483648):
                with self.subTest(name=name, value=value):
                    data = copy.deepcopy(self.data)
                    data[name][0]["reaction"] = value
                    self.assert_rejected(data, f"{name}/0/reaction", "")

    def test_format_units_and_missing_lists(self):
        data = copy.deepcopy(self.data)
        data["format"] = "other"
        self.assert_rejected(data, "format", "must be")
        data = copy.deepcopy(self.data)
        data["units"]["energy"] = "J/kmol"
        self.assert_rejected(data, "units", "must be")
        for name in LISTS:
            with self.subTest(name=name):
                data = copy.deepcopy(self.data)
                del data[name]
                self.assert_rejected(data, name, "required property")
        data = copy.deepcopy(self.data)
        for name in LISTS:
            data[name] = []
        exchange_json.validate(data, self.source)

    def test_optional_nulls_match_omitted_values_in_header(self):
        self.data["electron_translation"][0]["cross_section"] = {"model": "coulomb"}
        self.data["vibration_translation"][0]["relaxation_time"]["high_temperature_correction"] = {
            "model": "park_1985", "cross_section": 1e-20}
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = self.output.read_bytes()
        self.correction(self.data).update(reference_temperature=None, exponent=None)
        self.data["electron_translation"][0]["cross_section"].update(a=None, b=None, c=None, maximum_temperature=None)
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_bytes(), expected)
        time = self.data["vibration_translation"][0]["relaxation_time"]
        del time["high_temperature_correction"]
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = self.output.read_bytes()
        time["high_temperature_correction"] = None
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_bytes(), expected)

    def test_cli_rejection_preserves_existing_output(self):
        cases = []
        data = copy.deepcopy(self.data)
        del self.correction(data)["exponent"]
        cases.append((data, "vibration_translation 0"))
        data = copy.deepcopy(self.data)
        del data["electron_translation"][0]["cross_section"]["b"]
        cases.append((data, "electron_translation 0"))
        data = copy.deepcopy(self.data)
        data["chemistry_vibration"][0]["reaction"] = 0.5
        cases.append((data, "chemistry_vibration/0/reaction"))
        for data, context in cases:
            for existing in (False, True):
                with self.subTest(context=context, existing=existing):
                    self.data = data
                    if existing:
                        self.output.write_bytes(b"original\n")
                    elif self.output.exists():
                        self.output.unlink()
                    result = self.generate()
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(str(self.source), result.stderr)
                    self.assertIn(context, result.stderr)
                    self.assertNotIn("Traceback", result.stderr)
                    if existing:
                        self.assertEqual(self.output.read_bytes(), b"original\n")
                    else:
                        self.assertFalse(self.output.exists())

    def test_duplicate_exchange_keys_and_species_file_units(self):
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        previous = self.output.read_bytes()
        text = self.source.read_text().replace('"format":', '"format": "duplicate", "format":', 1)
        self.assertNotEqual(text, self.source.read_text())
        self.source.write_text(text)
        with self.assertRaises(SystemExit) as caught:
            exchange_json.load(self.source)
        self.assertIn("duplicate JSON key", str(caught.exception))
        data = json.loads(self.species_source.read_text())
        data["units"]["energy"] = "J/kmol"
        self.species_source.write_text(json.dumps(data))
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(str(self.species_source), result.stderr)
        self.assertIn("SI units", result.stderr)
        self.assertEqual(self.output.read_bytes(), previous)

    def test_shipped_files_pass_without_mutation(self):
        for name in ("air5_park_exchange", "air11_kim_jo_exchange"):
            with self.subTest(name=name):
                path = ROOT / "database/data" / f"{name}.json"
                original = json.loads(path.read_text())
                self.assertEqual(exchange_json.load(path), original)


if __name__ == "__main__":
    unittest.main()
