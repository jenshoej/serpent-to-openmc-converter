# SPDX-FileCopyrightText: 2024 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

from .serpent_preprocess import (
    INPUT_KEYWORDS,
    first_word,
    load_serpent_lines,
)

CardTokens = Tuple[List[str], int]


def gather_card_tokens(lines: List[str], index: int) -> CardTokens:
    tokens = lines[index].split()
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
    return tokens, scan_index


def iter_card_tokens(lines: List[str], keyword: str) -> List[List[str]]:
    tokens_list: List[List[str]] = []
    index = 0
    while index < len(lines):
        words = lines[index].split()
        if not words:
            index += 1
            continue
        if first_word(words) != keyword:
            index += 1
            continue
        tokens, next_index = gather_card_tokens(lines, index)
        tokens_list.append(tokens)
        index = next_index
    return tokens_list


def parse_surface_cards(lines: List[str]) -> Dict[str, Dict[str, Any]]:
    """Parse Serpent surface cards into a dict-based IR."""

    surfaces: Dict[str, Dict[str, Any]] = {}
    for tokens in iter_card_tokens(lines, "surf"):
        if len(tokens) < 3:
            continue
        name = tokens[1]
        coefficients: List[float] = []
        for value in tokens[3:]:
            try:
                coefficients.append(float(value))
            except ValueError as exc:
                raise ValueError(
                    f"Surface '{name}' has non-numeric parameter '{value}'."
                ) from exc
        surfaces[name] = {
            "type": "surface",
            "name": name,
            "surface_type": tokens[2],
            "coefficients": coefficients,
        }

    return surfaces


def parse_cell_cards(lines: List[str]) -> Dict[str, Dict[str, Any]]:
    """Parse Serpent cell cards into a dict-based IR."""

    cells: Dict[str, Dict[str, Any]] = {}
    for tokens in iter_card_tokens(lines, "cell"):
        if len(tokens) < 3:
            continue
        name = tokens[1]
        record: Dict[str, Any] = {
            "type": "cell",
            "name": name,
            "universe": tokens[2],
        }
        fill_mode = None
        region_tokens: List[str] = []
        if len(tokens) > 3:
            keyword = tokens[3]
            if keyword == "fill":
                fill_mode = "fill"
                if len(tokens) > 4:
                    record["fill"] = tokens[4]
                region_tokens = tokens[5:]
            elif keyword == "void":
                fill_mode = "void"
                region_tokens = tokens[4:]
            elif keyword == "outside":
                fill_mode = "outside"
                region_tokens = tokens[4:]
            else:
                fill_mode = "material"
                record["material"] = keyword
                region_tokens = tokens[4:]

        record["fill_mode"] = fill_mode
        record["region"] = " ".join(region_tokens) if region_tokens else ""
        cells[name] = record

    return cells


def parse_pin_cards(lines: List[str]) -> Dict[str, Dict[str, Any]]:
    """Parse Serpent pin cards into a dict-based IR."""

    pins: Dict[str, Dict[str, Any]] = {}
    for tokens in iter_card_tokens(lines, "pin"):
        if len(tokens) < 2:
            continue
        name = tokens[1]
        layers: List[Dict[str, Any]] = []
        index = 2
        while index < len(tokens):
            token = tokens[index]
            fill_mode = "material"
            if token == "void":
                fill_mode = "void"
                fill_name = None
                index += 1
            elif token == "fill":
                if index + 1 >= len(tokens):
                    raise ValueError(f"Pin '{name}' missing fill target.")
                fill_mode = "fill"
                fill_name = tokens[index + 1]
                index += 2
            else:
                fill_name = token
                index += 1

            radius = None
            if index < len(tokens):
                try:
                    radius = float(tokens[index])
                except ValueError:
                    radius = None
                else:
                    index += 1

            layers.append(
                {"fill": fill_name, "fill_mode": fill_mode, "radius": radius}
            )
        pins[name] = {
            "type": "pin",
            "name": name,
            "layers": layers,
        }

    return pins


def parse_lattice_cards(lines: List[str]) -> Dict[str, Dict[str, Any]]:
    """Parse Serpent lattice cards into a dict-based IR."""

    lattices: Dict[str, Dict[str, Any]] = {}
    for tokens in iter_card_tokens(lines, "lat"):
        if len(tokens) < 2:
            continue
        name = tokens[1]
        lat_type_raw = tokens[2] if len(tokens) > 2 else None
        try:
            lat_type = int(lat_type_raw) if lat_type_raw is not None else None
        except ValueError:
            lat_type = None
        record: Dict[str, Any] = {
            "type": "lattice",
            "name": name,
            "lat_type": lat_type_raw,
            "lat_type_int": lat_type,
        }

        if lat_type in (1, 2, 3, 14) and len(tokens) >= 8:
            record.update(
                {
                    "x0": float(tokens[3]),
                    "y0": float(tokens[4]),
                    "nx": int(tokens[5]),
                    "ny": int(tokens[6]),
                    "pitch": float(tokens[7]),
                    "entries": tokens[8:],
                }
            )
        elif lat_type in (6, 7, 8) and len(tokens) >= 7:
            record.update(
                {
                    "x0": float(tokens[3]),
                    "y0": float(tokens[4]),
                    "pitch": float(tokens[5]),
                    "entries": [tokens[6]],
                }
            )
        elif lat_type == 9 and len(tokens) >= 6:
            record.update(
                {
                    "x0": float(tokens[3]),
                    "y0": float(tokens[4]),
                    "n": int(tokens[5]),
                    "z_values": [float(value) for value in tokens[6::2]],
                    "stack_universes": tokens[7::2],
                }
            )
        lattices[name] = record

    return lattices


def parse_transformation_cards(lines: List[str]) -> List[Dict[str, Any]]:
    """Parse Serpent transformation cards into a list-based IR."""

    transformations: List[Dict[str, Any]] = []
    for keyword in ("trans", "transa", "transb", "transv"):
        for tokens in iter_card_tokens(lines, keyword):
            if len(tokens) < 2:
                continue
            record: Dict[str, Any] = {"type": "transformation", "card": keyword}
            if keyword == "trans":
                record["trans_type"] = tokens[1] if len(tokens) > 1 else None
                record["target"] = tokens[2] if len(tokens) > 2 else None
                record["values"] = tokens[3:] if len(tokens) > 3 else []
            else:
                record["trans_type"] = keyword
                record["target"] = tokens[1] if len(tokens) > 1 else None
                record["values"] = tokens[2:] if len(tokens) > 2 else []
            transformations.append(record)

    return transformations


def parse_set_cards(lines: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """Parse Serpent set cards needed for geometry expansion."""

    settings: Dict[str, List[Dict[str, Any]]] = {"usym": []}
    for line in lines:
        words = line.split()
        if not words or first_word(words) != "set":
            continue
        if len(words) < 2:
            continue
        if words[1].lower() != "usym":
            continue

        if len(words) < 3:
            continue

        record: Dict[str, Any] = {
            "type": "usym",
            "universe": words[2],
        }
        if len(words) > 3:
            record["symmetry"] = int(words[3])
        if len(words) > 4:
            record["symmetry_mode"] = int(words[4])
        if len(words) > 6:
            record["x0"] = float(words[5])
            record["y0"] = float(words[6])
        if len(words) > 8:
            record["angle0"] = float(words[7])
            record["angle"] = float(words[8])

        settings["usym"].append(record)

    return settings


def parse_geometry_ir(lines: List[str]) -> Dict[str, Any]:
    """Parse geometry-related cards into a lossless IR bundle.

    Notes
    -----
    ``lines`` should already be preprocessed (includes expanded, comments removed,
    continuation lines joined). Use ``geometry_ir_from_file`` for raw input files.
    """
    return {
        "surfaces": parse_surface_cards(lines),
        "cells": parse_cell_cards(lines),
        "pins": parse_pin_cards(lines),
        "lattices": parse_lattice_cards(lines),
        "transformations": parse_transformation_cards(lines),
        "settings": parse_set_cards(lines),
    }


def geometry_ir_from_file(path: str) -> Dict[str, Any]:
    """Load a Serpent input file and parse geometry cards into IR."""
    lines = load_serpent_lines(Path(path))
    return parse_geometry_ir(lines)
