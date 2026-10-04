"""Reading and writing the small subset of TOML a setup file uses.

Reading uses :mod:`tomllib` (Python 3.11+) or the ``tomli`` package it came from. Writing is done here so no
extra package is needed: tables, arrays of tables, strings, numbers, booleans and flat arrays.
"""

from __future__ import annotations

import json
from typing import Any

from .quantity import format_number

try:
    import tomllib as _toml_reader
except ImportError:  # Python < 3.11
    try:
        import tomli as _toml_reader  # type: ignore[no-redef]
    except ImportError:
        _toml_reader = None  # type: ignore[assignment]


class TOMLError(ValueError):
    """The text is not valid TOML."""


def loads(text: str) -> dict[str, Any]:
    if _toml_reader is None:
        raise ImportError("reading setup files needs Python 3.11+ or the 'tomli' package (pip install tomli)")
    try:
        return _toml_reader.loads(text)
    except _toml_reader.TOMLDecodeError as e:
        raise TOMLError(str(e)) from None


def _value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return format_number(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_value(x) for x in v) + "]"
    raise TypeError(f"cannot write {type(v).__name__} to TOML")


def _is_table_array(v: Any) -> bool:
    return isinstance(v, list) and bool(v) and all(isinstance(x, dict) for x in v)


def _table(name: str, d: dict[str, Any], out: list[str]) -> None:
    scalars = [(k, v) for k, v in d.items() if not isinstance(v, dict) and not _is_table_array(v)]
    if name:
        out.append("")
        out.append(f"[{name}]")
    out.extend(f"{k} = {_value(v)}" for k, v in scalars)
    for k, v in d.items():
        full = f"{name}.{k}" if name else k
        if isinstance(v, dict):
            _table(full, v, out)
        elif _is_table_array(v):
            for item in v:
                out.append("")
                out.append(f"[[{full}]]")
                out.extend(f"{ik} = {_value(iv)}" for ik, iv in item.items())


def dumps(d: dict[str, Any]) -> str:
    out: list[str] = []
    _table("", d, out)
    return "\n".join(out).lstrip("\n") + "\n"
