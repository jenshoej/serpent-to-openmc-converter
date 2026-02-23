# SPDX-FileCopyrightText: 2024 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

from __future__ import annotations

from pathlib import Path
from typing import List
import re
import shlex

INPUT_KEYWORDS = [
    "branch",
    "casematrix",
    "cell",
    "coef",
    "datamesh",
    "dep",
    "det",
    "div",
    "dtrans",
    "ene",
    "ftrans",
    "fun",
    "hisv",
    "ifc",
    "include",
    "lat",
    "mat",
    "mesh",
    "mflow",
    "mix",
    "nest",
    "particle",
    "pbed",
    "phb",
    "pin",
    "plot",
    "rep",
    "sample",
    "sens",
    "set",
    "solid",
    "src",
    "strans",
    "surf",
    "therm",
    "thermstoch",
    "tme",
    "trans",
    "transa",
    "transb",
    "transv",
    "umsh",
    "utrans",
    "wwgen",
    "wwin",
]


def first_word(input: str | List[str]) -> str:
    """Return the first token (lowercased) from a line or token list."""
    words = input.split() if isinstance(input, str) else input
    return words[0].lower()


def expand_include_cards(lines: List[str], basedir: Path) -> List[str]:
    """Inline Serpent 'include' cards by splicing in referenced file contents."""
    index = 0
    while True:
        if index >= len(lines):
            return lines

        words = lines[index].split()
        if words and first_word(words) == "include":
            include_path = Path(shlex.split(lines[index])[1])
            if not include_path.is_absolute():
                include_path = basedir / include_path
            with include_path.open("r") as fh:
                insert_lines = fh.readlines()
            lines[index : index + 1] = insert_lines
        else:
            index += 1


def remove_comments(lines: List[str]) -> List[str]:
    """Strip Serpent comments and blank lines from input."""
    text = "\n".join(line.rstrip("\n") for line in lines)
    text = re.sub(r"%.*$", "", text, flags=re.MULTILINE)
    return [line for line in text.splitlines(keepends=True) if line.strip()]


def join_lines(lines: List[str]) -> List[str]:
    """Join multi-line Serpent cards into a single logical line."""
    index = 0
    while True:
        if index >= len(lines):
            return lines

        words = lines[index].split()
        if not words:
            index += 1
            continue
        if first_word(words) in INPUT_KEYWORDS:
            while index + 1 < len(lines):
                if first_word(lines[index + 1]) in INPUT_KEYWORDS:
                    break
                lines[index] += lines.pop(index + 1)
        index += 1


def load_serpent_lines(path: Path) -> List[str]:
    """Load and preprocess Serpent input lines."""
    with path.open("r") as fh:
        all_lines = fh.readlines()
    all_lines = expand_include_cards(all_lines, path.parent)
    all_lines = remove_comments(all_lines)
    return join_lines(all_lines)
