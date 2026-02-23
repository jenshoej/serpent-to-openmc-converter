
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple
import re
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from src.ir_openmc import materials_from_ir
from src.serpent_mat_ir import parse_materials_to_ir, parse_therm_cards
from src.serpent_preprocess import load_serpent_lines

ComponentKey = Tuple[str, str]
ComponentValue = Tuple[float, str | None]
ComponentEntry = Tuple[ComponentKey, ComponentValue]


def get_entry_value(entry: Any, index: int, attr: str) -> Any:
    if hasattr(entry, attr):
        return getattr(entry, attr)
    if isinstance(entry, dict):
        return entry.get(attr)
    try:
        return entry[index]
    except (TypeError, IndexError):
        return None


def collect_openmc_components(material: Any) -> Dict[ComponentKey, ComponentValue]:
    components: Dict[ComponentKey, ComponentValue] = {}

    for entry in getattr(material, "nuclides", []) or getattr(material, "_nuclides", []):
        name = get_entry_value(entry, 0, "name")
        percent = get_entry_value(entry, 1, "percent")
        percent_type = get_entry_value(entry, 2, "percent_type")
        if name is None:
            continue
        components[("nuclide", str(name))] = (float(percent), percent_type)

    for entry in getattr(material, "elements", []) or getattr(material, "_elements", []):
        name = get_entry_value(entry, 0, "name")
        percent = get_entry_value(entry, 1, "percent")
        percent_type = get_entry_value(entry, 2, "percent_type")
        if name is None:
            continue
        components[("element", str(name))] = (float(percent), percent_type)

    return components


def element_from_nuclide(nuclide: str) -> str | None:
    match = re.match(r"^([A-Z][a-z]?)(?:-?\d+).*", nuclide)
    if not match:
        return None
    return match.group(1)


def collapse_element_components(
    ir_elements: Iterable[ComponentEntry],
    omc_components: Dict[ComponentKey, ComponentValue],
    rtol: float,
    atol: float,
) -> Tuple[List[ComponentKey], List[ComponentKey]]:
    matched_missing: List[ComponentKey] = []
    matched_extra: List[ComponentKey] = []

    omc_nuclides_by_element: Dict[str, List[ComponentEntry]] = {}
    for key, value in omc_components.items():
        if key[0] != "nuclide":
            continue
        element = element_from_nuclide(key[1])
        if element is None:
            continue
        omc_nuclides_by_element.setdefault(element, []).append((key, value))

    for (kind, element), (ir_percent, ir_mode) in ir_elements:
        if kind != "element":
            continue
        entries = omc_nuclides_by_element.get(element, [])
        if not entries:
            continue
        modes = {mode for _, (_, mode) in entries}
        if len(modes) != 1 or ir_mode not in modes:
            continue
        total = sum(percent for _, (percent, _) in entries)
        if not close(float(total), float(ir_percent), rtol, atol):
            continue
        matched_missing.append((kind, element))
        matched_extra.extend(key for key, _ in entries)

    return matched_missing, matched_extra


def collect_ir_components(record: Dict[str, Any]) -> Dict[ComponentKey, ComponentValue]:
    components: Dict[ComponentKey, ComponentValue] = {}
    for entry in record["nuclides"]:
        if entry["A"] is None:
            key = ("element", entry["element"])
        else:
            key = ("nuclide", entry["name"])
        components[key] = (abs(entry["percent"]), entry["percent_mode"])
    return components


def sab_names(material: Any) -> List[str]:
    sab_entries = []
    if hasattr(material, "s_alpha_beta"):
        sab_entries = material.s_alpha_beta
    elif hasattr(material, "sab"):
        sab_entries = material.sab

    names = []
    for entry in sab_entries:
        if isinstance(entry, str):
            names.append(entry)
        else:
            name = get_entry_value(entry, 0, "name")
            if name is not None:
                names.append(str(name))
    return names


def close(a: float, b: float, rtol: float, atol: float) -> bool:
    return abs(a - b) <= max(atol, rtol * abs(b))


def compare_material(
    name: str,
    record: Dict[str, Any],
    material: Any,
    rtol: float,
    atol: float,
) -> List[str]:
    issues: List[str] = []

    if record["density_mode"] != "sum":
        density = getattr(material, "density", None)
        if density is None:
            issues.append("missing density")
        else:
            if not close(float(density), float(record["density"]), rtol, atol):
                issues.append(
                    f"density mismatch (ir={record['density']}, openmc={density})"
                )
        units = getattr(material, "density_units", None)
        if units is not None and units != record["density_units"]:
            issues.append(
                f"density units mismatch (ir={record['density_units']}, openmc={units})"
            )

    if record["temperature"] is not None:
        temperature = getattr(material, "temperature", None)
        if temperature is None:
            issues.append("missing temperature")
        else:
            if not close(float(temperature), float(record["temperature"]), rtol, atol):
                issues.append(
                    f"temperature mismatch (ir={record['temperature']}, openmc={temperature})"
                )

    ir_sab = set(record["sab"])
    if ir_sab:
        omc_sab = set(sab_names(material))
        if ir_sab != omc_sab:
            issues.append(f"sab mismatch (ir={sorted(ir_sab)}, openmc={sorted(omc_sab)})")

    ir_components = collect_ir_components(record)
    omc_components = collect_openmc_components(material)

    missing = sorted(set(ir_components) - set(omc_components))
    extra = sorted(set(omc_components) - set(ir_components))
    if missing and extra:
        matched_missing, matched_extra = collapse_element_components(
            [(key, ir_components[key]) for key in missing],
            omc_components,
            rtol,
            atol,
        )
        missing = [key for key in missing if key not in matched_missing]
        extra = [key for key in extra if key not in matched_extra]
    if missing:
        issues.append(f"missing components: {missing}")
    if extra:
        issues.append(f"extra components: {extra}")

    for key, (percent, percent_mode) in ir_components.items():
        if key not in omc_components:
            continue
        omc_percent, omc_mode = omc_components[key]
        if omc_mode is not None and omc_mode != percent_mode:
            issues.append(
                f"{key} percent mode mismatch (ir={percent_mode}, openmc={omc_mode})"
            )
        if not close(float(omc_percent), float(percent), rtol, atol):
            issues.append(
                f"{key} percent mismatch (ir={percent}, openmc={omc_percent})"
            )

    return issues


def run_material_checks(
    input_file: Path | None = None,
    rtol: float = 1e-8,
    atol: float = 1e-12,
    show_all: bool = False,
) -> int:
    if input_file is None:
        input_file = SCRIPT_DIR / "SPX_case" / "SPX"

    lines = load_serpent_lines(input_file)
    therm = parse_therm_cards(lines)
    materials_ir = parse_materials_to_ir(lines, therm)
    openmc_materials = materials_from_ir(materials_ir)

    issues_found = 0
    skipped = 0
    for name, record in materials_ir.items():
        if record["type"] != "material":
            skipped += 1
            continue
        material = openmc_materials.get(name)
        if material is None:
            print(f"[missing] {name}")
            issues_found += 1
            continue
        issues = compare_material(name, record, material, rtol, atol)
        if issues:
            issues_found += 1
            print(f"[mismatch] {name}")
            for issue in issues:
                print(f"  - {issue}")
        elif show_all:
            print(f"[ok] {name}")

    if skipped:
        print(f"Skipped {skipped} mix materials (no direct IR composition to compare).")

    if issues_found:
        print(f"Found {issues_found} material mismatches.")
        return 1

    print("All materials matched.")
    return 0
