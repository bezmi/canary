# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
"""Shared output structure of canary species files (format canary-species-2).

Numbers are written by Python's json module in the shortest form that reads
back to the same double. Decimal strings from source files are converted with
float(), and molar masses given in g/mol are divided by 1000 in decimal before
rounding, so that "28.0134" g/mol and "28.0134e-3" kg/mol give the same double.
"""

import json
import math
import os
import secrets
import stat
from decimal import Decimal
from pathlib import Path

FORMAT = "canary-species-2"
UNITS = {"molar_mass": "kg/mol", "energy": "J/mol", "temperature": "K", "length": "m"}


def kg_per_mol(grams_per_mol: str) -> float:
    return float(Decimal(grams_per_mol.replace("D", "E").replace("d", "e")) / Decimal(1000))


def number(text: str) -> float:
    return float(text.strip().replace("D", "E").replace("d", "e"))


def nasa9_thermo(ranges, coefficients, blend_widths, reference=None, notes=None):
    if len(ranges) < 2:
        raise ValueError("NASA9 requires at least two temperature boundaries")
    if len(coefficients) != len(ranges) - 1:
        raise ValueError("NASA9 coefficient rows must match the temperature intervals")
    if any(len(row) != 9 for row in coefficients):
        raise ValueError("NASA9 requires nine coefficients per temperature interval")
    if len(blend_widths) != len(ranges) - 2:
        raise ValueError("NASA9 blend widths must match the interior temperature boundaries")
    thermo = {
        "model": "nasa9",
        "temperature_ranges": ranges,
        "blend_widths": blend_widths,
        "coefficients": coefficients,
    }
    if reference:
        thermo["reference"] = reference
    if notes:
        thermo["notes"] = notes
    return thermo


def blend_widths_for(name, ranges, widths_by_boundary):
    """Blend widths for the interior boundaries of ranges, looked up by
    boundary temperature in widths_by_boundary ({temperature: width})."""
    widths = []
    for boundary in ranges[1:-1]:
        if boundary not in widths_by_boundary:
            raise SystemExit(f"{name}: no blend width given for the boundary at {boundary} K")
        widths.append(widths_by_boundary[boundary])
    return widths


def parse_blend_widths(items):
    """'1000=400' style arguments to {1000.0: 400.0}."""
    widths = {}
    for item in items:
        boundary, width = item.split("=")
        widths[float(boundary)] = float(width)
    return widths


def parse_roles(items, names):
    """Require exactly one explicit heavy/electron role per selected species."""
    selected = set(names)
    if len(selected) != len(names):
        raise SystemExit("selected species must be distinct")
    roles = {}
    for item in items:
        name, separator, role = item.partition("=")
        if not separator or not name or role not in ("heavy", "electron"):
            raise SystemExit(f"invalid role declaration {item!r}: expected NAME=heavy or NAME=electron")
        if name not in selected:
            raise SystemExit(f"role for {name}: species is not selected")
        if name in roles:
            raise SystemExit(f"role for {name}: given twice")
        roles[name] = role
    for name in names:
        if name not in roles:
            raise SystemExit(f"role for {name}: missing declaration")
    return roles


def element_counts(entries, name):
    """Validate nuclear counts and retain decoded element identifiers."""
    result = {}
    for element, count in entries:
        if not isinstance(element, str) or not element:
            raise SystemExit(f"{name}: element identifier must not be empty")
        if element in result:
            raise SystemExit(f"{name}: element {element}: given twice")
        if isinstance(count, bool) or not isinstance(count, (int, float, Decimal)):
            raise SystemExit(f"{name}: element {element}: count must be a positive unsigned integer")
        if ((isinstance(count, Decimal) and not count.is_finite())
                or not 0 < count <= 4294967295 or count != int(count)):
            raise SystemExit(f"{name}: element {element}: count must be a positive unsigned integer")
        result[element] = int(count)
    return result


def validate_elements(record, name):
    elements = record.get("elements")
    if elements is None:
        return
    if not isinstance(elements, dict):
        raise SystemExit(f"{name}: elements must be an object")
    validated = element_counts(elements.items(), name)
    if record.get("role") == "electron" and validated:
        raise SystemExit(f"{name}: free-electron elemental composition must be empty")


def load(path):
    """Read JSON without silently replacing duplicate object keys."""
    def unique_object(entries):
        out = {}
        for key, value in entries:
            if key in out:
                raise SystemExit(f"{path}: duplicate JSON key {key!r}")
            out[key] = value
        return out
    with open(path, encoding="utf-8") as stream:
        return json.load(stream, object_pairs_hook=unique_object)


def load_validated(path):
    """Read a species file and check its schema, format, units and metadata.

    Conversion-specific dimensions and numerical domains require subsequent
    validation. Import schema support here so numeric/output helpers retain
    their standard-library-only dependency.
    """
    from _schema_validation import validate_structure

    data = load(path)
    if isinstance(data, dict) and "format" in data and data["format"] != FORMAT:
        raise SystemExit(f"{path}: format must be {FORMAT}; regenerate with explicit species roles")
    validate_structure(data, "species.schema.json", path)
    if data["units"] != UNITS:
        raise SystemExit(f"{path}: not a {FORMAT} file in SI units")
    for name, record in data["species"].items():
        context = f"{path}: species {name}"
        if record["role"] == "electron" and record["charge"] != -1:
            raise SystemExit(f"{context}: electron role requires charge -1")
        validate_elements(record, context)
    return data


def validate_selection(species, names, source):
    """Require distinct selected species that are present in the source file."""
    if len(set(names)) != len(names):
        raise SystemExit(f"{source}: selected species must be distinct")
    for name in names:
        if name not in species:
            raise SystemExit(f"{source}: species {name!r} is not in the file")


def _validate_boundaries(values, context):
    # The C++ host input stores boundaries as double, including integer literals.
    boundaries = [float(value) for value in values]
    for index, value in enumerate(boundaries):
        if not math.isfinite(value) or value <= 0 or (index and value <= boundaries[index - 1]):
            raise SystemExit(f"{context}: temperature boundary {index} must be finite, positive and increasing")
    return boundaries


def _validate_coefficients(rows, context):
    for row_index, row in enumerate(rows):
        for column, value in enumerate(row):
            if not math.isfinite(value):
                raise SystemExit(f"{context}: coefficient {column} of row {row_index} must be finite")


def validate_nasa9(record, context, segments):
    """Validate species SI inputs and NASA9 dimensions, boundaries and blends."""
    if not math.isfinite(record["molar_mass"]) or record["molar_mass"] <= 0:
        raise SystemExit(f"{context}: molar_mass must be finite and positive")
    if not math.isfinite(record["charge"]):
        raise SystemExit(f"{context}: charge must be finite")
    thermo = record["thermo"]
    if thermo["model"] != "nasa9":
        raise SystemExit(f"{context}: thermo model must be nasa9")
    try:
        nasa9_thermo(thermo["temperature_ranges"], thermo["coefficients"], thermo["blend_widths"])
    except ValueError as error:
        raise SystemExit(f"{context}: {error}") from None
    if not 1 <= len(thermo["coefficients"]) <= segments:
        raise SystemExit(f"{context}: NASA9 segments must fit capacity {segments}")
    boundaries = _validate_boundaries(thermo["temperature_ranges"], f"{context}: NASA9")
    previous_end = boundaries[0]
    for index, source_width in enumerate(thermo["blend_widths"]):
        width = float(source_width)
        if not math.isfinite(width) or width <= 0:
            raise SystemExit(f"{context}: blend width {index} must be finite and positive")
        boundary = boundaries[index + 1]
        if boundary - 0.5 * width < previous_end:
            raise SystemExit(f"{context}: blend region {index} starts before the preceding region or fit boundary")
        previous_end = boundary + 0.5 * width
    if previous_end > boundaries[-1]:
        raise SystemExit(f"{context}: final blend region extends beyond the fit boundary")
    _validate_coefficients(thermo["coefficients"], f"{context}: NASA9")


def validate_electronic_levels(record, context, levels):
    """Validate electronic level dimensions and SI domains after schema checks."""
    data = record.get("electronic_levels")
    if data is None:
        return
    if len(data["degeneracies"]) != len(data["temperatures"]):
        raise SystemExit(f"{context}: electronic degeneracies and temperatures must have equal lengths")
    if len(data["degeneracies"]) > levels:
        raise SystemExit(f"{context}: electronic levels exceed capacity {levels}")
    for index, (degeneracy, theta) in enumerate(zip(data["degeneracies"], data["temperatures"])):
        if degeneracy <= 0:
            raise SystemExit(f"{context}: electronic level {index}: degeneracy must be positive")
        if not math.isfinite(theta) or theta < 0:
            raise SystemExit(f"{context}: electronic level {index}: temperature must be finite and non-negative")
        if index == 0 and theta != 0:
            raise SystemExit(f"{context}: electronic level 0: ground temperature must be zero")


def validate_transport(record, context, segments):
    """Validate required transport records, CEA dimensions and SI input domains."""
    transport = record.get("transport")
    fields = ("viscosity", "thermal_conductivity", "lennard_jones")
    if transport is None or any(transport.get(field) is None for field in fields):
        raise SystemExit(f"{context}: needs viscosity, thermal_conductivity and lennard_jones transport data")
    for field in ("sigma", "epsilon"):
        value = transport["lennard_jones"][field]
        if not math.isfinite(value) or value <= 0:
            raise SystemExit(f"{context}: Lennard-Jones {field} must be finite and positive")
    for field in ("viscosity", "thermal_conductivity"):
        curve = transport[field]
        location = f"{context}: {field}"
        if curve["model"] != "cea":
            raise SystemExit(f"{location}: transport curve model must be cea")
        rows = curve["coefficients"]
        if not 1 <= len(rows) <= segments:
            raise SystemExit(f"{location}: CEA segments must fit capacity {segments}")
        if len(curve["temperature_ranges"]) != len(rows) + 1:
            raise SystemExit(f"{location}: temperature boundaries must number one more than coefficient rows")
        for index, row in enumerate(rows):
            if len(row) != 4:
                raise SystemExit(f"{location}: coefficient row {index} requires exactly four values")
        _validate_boundaries(curve["temperature_ranges"], location)
        _validate_coefficients(rows, location)


def species_file(species, reference=None, notes=None):
    for name, record in species.items():
        if record.get("role") not in ("heavy", "electron"):
            raise SystemExit(f"{name}: role must be heavy or electron")
        if record["role"] == "electron" and record["charge"] != -1:
            raise SystemExit(f"{name}: electron role requires charge -1")
        validate_elements(record, name)
    out = {"format": FORMAT, "units": UNITS, "species": species}
    if reference:
        out["reference"] = reference
    if notes:
        out["notes"] = notes
    return out


def _destination_permissions(path):
    try:
        mode = path.lstat().st_mode
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(mode):
        raise OSError("destination must be a regular file; symlinks are not supported")
    return stat.S_IMODE(mode)


def _open_temporary(path):
    for _ in range(128):
        temporary = path.parent / f".canary-{secrets.token_hex(16)}.tmp"
        try:
            return temporary, open(temporary, "x", encoding="utf-8")
        except FileExistsError:
            continue
    raise OSError("cannot create a unique temporary file")


def dump(data, path):
    """Serialize finite JSON, then atomically replace a regular-file destination.

    Preserve existing permission bits and reject symlinks. Failures before
    replacement leave the destination unchanged and attempt temporary cleanup.
    Serialization raises TypeError or ValueError; filesystem failures raise
    OSError with the destination and operation. Parent directories must exist.
    Callers coordinate concurrent path changes; no power-loss durability is
    provided.
    """
    text = json.dumps(data, indent=3, allow_nan=False) + "\n"
    path = Path(path)
    temporary = None
    stream = None
    operation = "inspect destination"
    try:
        permissions = _destination_permissions(path)
        operation = "create temporary file"
        temporary, stream = _open_temporary(path)
        # Restrict access before writing contents; writes can clear special bits.
        operation = "set temporary file permissions"
        if permissions is not None:
            os.chmod(temporary, permissions)
        operation = "write output"
        if stream.write(text) != len(text):
            raise OSError("incomplete write")
        operation = "flush output"
        stream.flush()
        operation = "close output"
        stream.close()
        operation = "inspect destination"
        permissions = _destination_permissions(path)
        operation = "set replacement permissions"
        if permissions is not None:
            os.chmod(temporary, permissions)
        operation = "replace destination"
        os.replace(temporary, path)
        temporary = None
    except OSError as error:
        raise OSError(f"{path}: cannot {operation}: {error}") from error
    finally:
        # Cleanup must not obscure the failed operation.
        if stream is not None and not stream.closed:
            try:
                stream.close()
            except OSError:
                pass
        if temporary is not None:
            try:
                temporary.unlink()
            except OSError:
                pass
