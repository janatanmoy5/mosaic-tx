"""Strict CSV input and output helpers."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Iterable


class DataContractError(ValueError):
    """Raised when an input violates a documented table contract."""


def read_rows(path: str | Path, required: Iterable[str]) -> list[dict[str, str]]:
    source = Path(path)
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = set(required) - fields
        if missing:
            raise DataContractError(f"{source}: missing columns {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise DataContractError(f"{source}: table is empty")
    return rows


def finite_number(row: dict[str, str], field: str, context: str) -> float:
    import math

    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise DataContractError(f"{context}: {field} must be numeric") from exc
    if not math.isfinite(value):
        raise DataContractError(f"{context}: {field} must be finite")
    return value


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: str | Path, rows: list[dict[str, object]]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty result")
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

