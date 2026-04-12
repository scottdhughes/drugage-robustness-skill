"""Utilities."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


def read_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def now_timestamp() -> str:
    return datetime.now(UTC).isoformat()


def collapse_whitespace(value: str) -> str:
    return " ".join(value.strip().split())


def normalize_key(value: str) -> str:
    return collapse_whitespace(value).casefold()


def trimmed_mean(values: list[float] | np.ndarray, proportion: float) -> float:
    array = np.sort(np.asarray(values, dtype=float))
    if array.size == 0:
        return math.nan
    trim = int(math.floor(array.size * proportion))
    if trim * 2 >= array.size:
        return float(array.mean())
    trimmed = array[trim : array.size - trim]
    return float(trimmed.mean())


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def package_versions(names: list[str]) -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in names:
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def runtime_environment(package_names: list[str], pythonhashseed: int) -> dict[str, Any]:
    return {
        "python_version": platform.python_version(),
        "package_versions": package_versions(package_names),
        "pythonhashseed": pythonhashseed,
    }


def set_runtime_environment(pythonhashseed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(pythonhashseed)


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
