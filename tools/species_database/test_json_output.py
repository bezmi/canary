# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
"""Check finite JSON serialization and atomic file replacement failures."""

import json
import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest import mock

import species_json


class JsonOutputTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.path = self.directory / "output.json"
        self.original = b'{"original": true}\n'
        self.path.write_bytes(self.original)
        self.data = {"nested": [1.25, {"label": "temperature \u0398"}], "zero": -0.0}

    def assert_preserved(self):
        self.assertEqual(self.path.read_bytes(), self.original)
        self.assertEqual(list(self.directory.glob(".canary-*.tmp")), [])

    def test_complete_serialization_precedes_filesystem_access(self):
        cycle = []
        cycle.append(cycle)
        for value, error in [(float("nan"), ValueError), (float("inf"), ValueError),
                             (-float("inf"), ValueError), (object(), TypeError),
                             (cycle, ValueError)]:
            with self.subTest(value=type(value).__name__, error=error):
                with mock.patch.object(Path, "lstat", side_effect=AssertionError("filesystem inspected")):
                    with self.assertRaises(error):
                        species_json.dump({"nested": [value]}, self.path)
                self.assert_preserved()
        self.path.unlink()
        with self.assertRaises(ValueError):
            species_json.dump({"nested": [float("nan")]}, self.path)
        self.assertEqual(list(self.directory.iterdir()), [])

    def test_replacement_preserves_format_and_old_hard_link(self):
        link = self.directory / "old.json"
        os.link(self.path, link)
        species_json.dump(self.data, str(self.path))
        self.assertEqual(self.path.read_bytes(), (json.dumps(self.data, indent=3) + "\n").encode())
        self.assertEqual(link.read_bytes(), self.original)
        self.assertEqual(list(self.directory.glob(".canary-*.tmp")), [])

    @unittest.skipUnless(os.name == "posix", "requires POSIX permission bits")
    def test_existing_permissions_apply_before_contents_and_after_close(self):
        self.path.chmod(0o640)
        original_open = species_json._open_temporary
        observed = []

        def open_temporary(path):
            temporary, stream = original_open(path)
            wrapper = mock.Mock(wraps=stream)
            wrapper.closed = False

            def write(text):
                observed.append(stat.S_IMODE(temporary.stat().st_mode))
                return stream.write(text)

            wrapper.write.side_effect = write
            return temporary, wrapper

        with mock.patch.object(species_json, "_open_temporary", side_effect=open_temporary):
            species_json.dump(self.data, self.path)
        self.assertEqual(observed, [0o640])
        self.assertEqual(stat.S_IMODE(self.path.stat().st_mode), 0o640)

    @unittest.skipUnless(os.name == "posix", "requires POSIX creation permissions")
    def test_new_file_uses_normal_creation_permissions(self):
        self.path.unlink()
        old_mask = os.umask(0o027)
        try:
            species_json.dump(self.data, self.path)
        finally:
            os.umask(old_mask)
        self.assertEqual(stat.S_IMODE(self.path.stat().st_mode), 0o640)
        self.assertEqual(json.loads(self.path.read_text()), self.data)

    def test_symlinks_are_rejected_including_dangling_links(self):
        link = self.directory / "link.json"
        for target in (self.path, self.directory / "missing.json"):
            with self.subTest(target=target.name):
                link.symlink_to(target)
                with self.assertRaisesRegex(OSError, "symlinks are not supported"):
                    species_json.dump(self.data, link)
                self.assertTrue(link.is_symlink())
                self.assert_preserved()
                link.unlink()
        self.assertFalse((self.directory / "missing.json").exists())

    def test_invalid_destinations(self):
        for path in (self.directory, self.directory / "missing" / "out.json", self.path / "out.json"):
            with self.subTest(path=path):
                with self.assertRaises(OSError) as caught:
                    species_json.dump(self.data, path)
                self.assertIn(str(path), str(caught.exception))
                self.assert_preserved()

    def test_write_flush_and_close_failures(self):
        original_open = species_json._open_temporary
        for operation in ("write", "short write", "flush", "close"):
            with self.subTest(operation=operation):
                streams = []

                def open_temporary(path):
                    temporary, stream = original_open(path)
                    streams.append(stream)
                    wrapper = mock.Mock(wraps=stream)
                    wrapper.closed = False
                    if operation == "short write":
                        wrapper.write.side_effect = lambda text: stream.write(text[:3])
                    else:
                        getattr(wrapper, operation).side_effect = OSError("injected failure")
                    # Permit cleanup after the checked close reports its failure.
                    if operation == "close":
                        def close():
                            stream.close()
                            raise OSError("injected failure")
                        wrapper.close.side_effect = close
                    return temporary, wrapper

                with mock.patch.object(species_json, "_open_temporary", side_effect=open_temporary):
                    with self.assertRaises(OSError) as caught:
                        species_json.dump(self.data, self.path)
                self.assertIn(f"cannot {operation.replace('short ', '')} output", str(caught.exception))
                self.assertIn(str(self.path), str(caught.exception))
                self.assertTrue(all(stream.closed for stream in streams))
                self.assert_preserved()

    def test_filesystem_failures(self):
        for name, operation in (("_open_temporary", "create temporary file"),
                                ("os.chmod", "set temporary file permissions"),
                                ("os.replace", "replace destination")):
            with self.subTest(operation=operation):
                with mock.patch(f"species_json.{name}", side_effect=OSError("injected failure")):
                    with self.assertRaisesRegex(OSError, f"cannot {operation}"):
                        species_json.dump(self.data, self.path)
                self.assert_preserved()

    def test_late_destination_inspection_and_permission_failures(self):
        for name, operation in (("_destination_permissions", "inspect destination"),
                                ("os.chmod", "set replacement permissions")):
            with self.subTest(operation=operation):
                with mock.patch(f"species_json.{name}", side_effect=[0o644, OSError("injected failure")]):
                    with self.assertRaisesRegex(OSError, f"cannot {operation}"):
                        species_json.dump(self.data, self.path)
                self.assert_preserved()

    def test_temporary_collisions_do_not_overwrite_other_files(self):
        occupied = self.directory / ".canary-occupied.tmp"
        occupied.write_text("unrelated")
        with mock.patch.object(species_json.secrets, "token_hex", side_effect=["occupied", "available"]):
            species_json.dump(self.data, self.path)
        self.assertEqual(occupied.read_text(), "unrelated")
        self.assertFalse((self.directory / ".canary-available.tmp").exists())
        self.path.write_bytes(self.original)
        with mock.patch.object(species_json.secrets, "token_hex", return_value="occupied"):
            with self.assertRaisesRegex(OSError, "cannot create a unique temporary file"):
                species_json.dump(self.data, self.path)
        self.assertEqual(occupied.read_text(), "unrelated")
        occupied.unlink()
        self.assert_preserved()


if __name__ == "__main__":
    unittest.main()
