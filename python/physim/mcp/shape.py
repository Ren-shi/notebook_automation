"""What goes back to the model: the planner's results made small, plain and JSON-safe.

The :class:`~physim.nuclear.planner.Planner` returns NumPy arrays, dataclasses, paths and quantities. A model reads
JSON, and a spectrum of hundreds of bins fills its context while saying little, so :func:`plain` applies the same
rules to every result:

- arrays and number lists longer than :data:`LIMIT` become a summary (length, range, mean);
- floats keep 6 significant figures; NaN and infinity become None;
- dataclasses become dicts, paths and quantities text.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path
from typing import Any

import numpy as np

#: The longest array or list of numbers sent as it is.
LIMIT = 50


def _float(x: float):
    x = float(x)
    if not math.isfinite(x):
        return None
    return float(f"{x:.6g}")


def _summary(values: np.ndarray) -> dict:
    finite = values[np.isfinite(values)] if values.dtype.kind in "fc" else values
    out: dict[str, Any] = {"summary": "array", "length": int(values.size)}
    if values.ndim > 1:
        out["shape"] = list(values.shape)
    if finite.size:
        out.update(min=_float(finite.min()), max=_float(finite.max()), mean=_float(finite.mean()))
    return out


def plain(x: Any, limit: int = LIMIT) -> Any:
    """``x`` as JSON-safe data, with long arrays summarised (see the module's rules)."""
    if x is None or isinstance(x, (bool, str)):
        return x
    if isinstance(x, (int, np.integer)) and not isinstance(x, np.bool_):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return _float(x)
    if isinstance(x, np.ndarray):
        if x.dtype.kind in "biuf":
            return _summary(x) if x.size > limit else [plain(v, limit) for v in x.tolist()]
        return [plain(v, limit) for v in x.tolist()[:limit]]
    if isinstance(x, dict):
        return {str(k): plain(v, limit) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        items = list(x)
        if len(items) > limit and all(isinstance(v, (int, float, np.number)) and not isinstance(v, bool)
                                      for v in items):
            return _summary(np.asarray(items, dtype=float))
        return [plain(v, limit) for v in items]
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return plain(dataclasses.asdict(x), limit)
    if isinstance(x, Path):
        return str(x)
    return str(x)
