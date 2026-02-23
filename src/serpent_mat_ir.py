
from __future__ import annotations

from typing import Any, Dict, List, Tuple

from openmc.data import get_thermal_name
from openmc.data.ace import get_metadata

from .serpent_preprocess import (
    INPUT_KEYWORDS,
    first_word,
    join_lines,
    remove_comments,
)

MATERIAL_KEYWORD_PARAMS: Dict[str, int] = {
    "tmp": 2,
    "tms": 2,
    "tft": 3,
    "rgb": 4,
    "vol": 2,
    "mass": 2,
    "burn": 2,
    "fix": 3,
    "moder": 3,
}


def parse_therm_cards(lines: List[str]) -> Dict[str, str]:
    """Parse thermal scattering cards into a name -> SAB mapping."""
    therm_materials: Dict[str, str] = {}
    for line in lines:
        words = line.split()
        if not words or first_word(words) != "therm":
            continue

        if len(words) > 3:
            therm_data = {"name": words[1], "temp": words[2], "lib": words[3:]}
        else:
            therm_data = {"name": words[1], "lib": words[2:]}

        name = therm_data["lib"][0]
        if "." in name:
            name, _ = name.split(".")
        therm_materials[therm_data["name"]] = get_thermal_name(name)

    return therm_materials


def parse_zaid(nuclide: str) -> Tuple[str, str | None, str | None, int | None, bool]:
    """Parse a Serpent nuclide token into metadata fields.

    Returns (name, element, zaid, A, metastable_ignored).
    """
    metastable_ignored = False
    if "." in nuclide:
        zaid, _ = nuclide.split(".")
    else:
        zaid = nuclide

    if "-" in zaid:
        element, A = zaid.split("-")
        if A == "nat":
            name = zaid.translate({ord("-"): None})
            return name, element, zaid, None, metastable_ignored
        if "m" in A[-1]:
            A = A.translate({ord("m"): None})
            zaid = zaid[:-1]
            metastable_ignored = True
        A = int(A)
        name = zaid.translate({ord("-"): None})
        return name, element, zaid, A, metastable_ignored

    name, element, _, A, _ = get_metadata(int(zaid))
    return name, element, zaid, A, metastable_ignored


def parse_materials_to_ir(
    lines: List[str],
    therm_materials: Dict[str, str] | None = None,
    preprocess: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """Parse Serpent material and mix cards into a dict-based IR.

    Returns
    -------
    dict
        Mapping from Serpent material/mix names to structured IR records.
    """
    if therm_materials is None:
        therm_materials = {}

    if preprocess:
        lines = join_lines(remove_comments(lines))

    materials: Dict[str, Dict[str, Any]] = {}

    index = 0
    while index < len(lines):
        words = lines[index].split()
        if not words:
            index += 1
            continue

        keyword = first_word(words)

        if keyword == "mat":
            tokens = list(words)
            scan_index = index + 1
            while scan_index < len(lines):
                next_words = lines[scan_index].split()
                if not next_words:
                    scan_index += 1
                    continue
                if first_word(next_words) in INPUT_KEYWORDS:
                    break
                tokens.extend(next_words)
                scan_index += 1

            if len(tokens) < 3:
                index = scan_index
                continue

            name = tokens[1]
            material_id = int(name) if name.isnumeric() else None
            density_token = tokens[2]
            if density_token == "sum":
                density = None
                density_units = None
                density_mode = "sum"
            else:
                density = float(density_token)
                if density < 0:
                    density = abs(density)
                    density_units = "g/cm3"
                    density_mode = "mass"
                else:
                    density_units = "atom/b-cm"
                    density_mode = "atom"

            record: Dict[str, Any] = {
                "type": "material",
                "name": name,
                "material_id": material_id,
                "density": density,
                "density_units": density_units,
                "density_mode": density_mode,
                "temperature": None,
                "sab": [],
                "keywords": {},
                "nuclides": [],
            }

            token_index = 3
            while token_index < len(tokens):
                mat_keyword = tokens[token_index].lower()
                if mat_keyword in MATERIAL_KEYWORD_PARAMS:
                    end_index = token_index + MATERIAL_KEYWORD_PARAMS[mat_keyword]
                    params = tokens[token_index + 1 : end_index]
                    record["keywords"][mat_keyword] = params
                    if mat_keyword == "tmp" and params:
                        record["temperature"] = float(params[0])
                    elif mat_keyword == "moder" and params:
                        sab = therm_materials.get(params[0], params[0])
                        record["sab"].append(sab)
                    token_index = end_index
                else:
                    break

            for nuclide, percent in zip(
                tokens[token_index::2], tokens[token_index + 1 :: 2]
            ):
                percent_value = float(percent)
                nuclide_name, element, zaid, A, metastable_ignored = parse_zaid(
                    nuclide
                )
                record["nuclides"].append(
                    {
                        "token": nuclide,
                        "zaid": zaid,
                        "name": nuclide_name,
                        "element": element,
                        "A": A,
                        "percent": percent_value,
                        "percent_mode": "wo" if percent_value < 0 else "ao",
                        "metastable_ignored": metastable_ignored,
                    }
                )

            materials[name] = record
            index = scan_index
            continue

        if keyword == "mix":
            name = words[1]
            components = []
            for mix_id, percent in zip(words[2::2], words[3::2]):
                components.append({"material": mix_id, "percent": float(percent) / 100.0})
            mix_type = "vo" if components and components[0]["percent"] > 0 else "wo"
            materials[name] = {
                "type": "mix",
                "name": name,
                "material_id": int(name) if name.isnumeric() else None,
                "components": components,
                "mix_type": mix_type,
                "keywords": {},
            }

        index += 1

    return materials
