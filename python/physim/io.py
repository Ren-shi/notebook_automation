"""Saving and loading checkpoints and trajectories.

Checkpoints (``.npz``) hold everything needed to continue a run bit for bit:
``save_checkpoint(world, path)`` now, ``load_checkpoint(path)`` later.

Trajectories go to ``.npz`` (always available) or ``.h5``/``.hdf5`` (needs ``h5py``).
:func:`save_trajectory` writes one you already have; :class:`TrajectoryWriter` is a
``sink`` for :meth:`World.run` that writes frames as they are produced, so memory stays
flat however long the run::

    with ps.TrajectoryWriter("orbit.h5") as out:
        world.run(dt=1e-3, steps=10_000_000, sink=out)
    traj = ps.load_trajectory("orbit.h5")

Both file kinds store the run's metadata (engine version, integrator, dt, forces...) as
JSON, readable without physim.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Mapping, Optional, Union

import numpy as np

from ._core import Trajectory, World

PathLike = Union[str, "os.PathLike[str]"]

_CHECKPOINT_ARRAYS = ("positions", "velocities", "masses", "pinned")
# Arrays added in later versions; files without them still load.
_OPTIONAL_CHECKPOINT_ARRAYS = ("charges",)

# name: (dtype, per-row shape given n_particles and n_constraints)
_FIELDS = {
    "t": (np.float64, lambda n, c: ()),
    "pos": (np.float64, lambda n, c: (n, 3)),
    "vel": (np.float64, lambda n, c: (n, 3)),
    "kinetic": (np.float64, lambda n, c: ()),
    "potential": (np.float64, lambda n, c: ()),
    "tension": (np.float64, lambda n, c: (c,)),
    "event_t": (np.float64, lambda n, c: ()),
    "event_index": (np.int64, lambda n, c: ()),
    "event_pos": (np.float64, lambda n, c: (n, 3)),
    "event_vel": (np.float64, lambda n, c: (n, 3)),
}
_FRAME_FIELDS = ("t", "pos", "vel", "kinetic", "potential", "tension")
_EVENT_FIELDS = ("event_t", "event_index", "event_pos", "event_vel")


def _npz_path(path: PathLike) -> Path:
    path = Path(path)
    if path.suffix != ".npz":
        raise ValueError(f"{path}: expected a .npz file")
    return path


def _format(path: PathLike) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".npz":
        return "npz"
    if suffix in (".h5", ".hdf5"):
        return "hdf5"
    raise ValueError(f"{path}: unsupported extension; use .npz, .h5 or .hdf5")


def _h5py():
    try:
        import h5py
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ImportError("HDF5 files need h5py: pip install h5py") from e
    return h5py


# ---------------------------------------------------------------------------
# Checkpoints


def save_checkpoint(world: World, path: PathLike) -> None:
    """Writes ``world.checkpoint()`` to a ``.npz`` file."""
    path = _npz_path(path)
    data = world.checkpoint()
    arrays = {k: np.asarray(data.pop(k)) for k in _CHECKPOINT_ARRAYS}
    for k in _OPTIONAL_CHECKPOINT_ARRAYS:
        if k in data:
            arrays[k] = np.asarray(data.pop(k))
    arrays["t"] = np.float64(data.pop("t"))
    meta = json.dumps(data)
    np.savez(path, meta=np.array(meta), **arrays)


def load_checkpoint(
    path: PathLike, custom_forces: Optional[Mapping[Any, Any]] = None
) -> World:
    """Rebuilds the World saved by :func:`save_checkpoint`.

    ``custom_forces`` maps the name (or id) of each ``CustomForce`` in the saved world to
    the force to use now, since Python functions cannot be saved.
    """
    with np.load(_npz_path(path), allow_pickle=False) as f:
        data = json.loads(str(f["meta"]))
        data.update({k: f[k] for k in _CHECKPOINT_ARRAYS})
        data.update({k: f[k] for k in _OPTIONAL_CHECKPOINT_ARRAYS if k in f})
        data["t"] = float(f["t"])
    return World.from_checkpoint(data, dict(custom_forces or {}))


# ---------------------------------------------------------------------------
# Trajectories


class TrajectoryWriter:
    """Appends trajectory chunks to a file; pass it as ``sink`` to :meth:`World.run`.

    Several runs may write to the same writer one after another (same particle count and
    energy setting); a frame repeated at a run boundary is stored twice. The file is
    complete once :meth:`close` is called (or the ``with`` block ends).

    ``.npz`` output is staged in raw files next to the target and packed on close;
    ``.h5`` output is appended in place.
    """

    def __init__(self, path: PathLike, metadata: Optional[Mapping[str, Any]] = None):
        self.path = Path(path)
        self._format = _format(self.path)
        if self._format == "hdf5":
            _h5py()
        self._metadata = dict(metadata) if metadata is not None else None
        self._n: Optional[int] = None
        self._c = 0
        self._energies: Optional[bool] = None
        self._terminated_by: Optional[int] = None
        self._counts = {k: 0 for k in _FIELDS}
        self._closed = False
        self._h5 = None
        self._tmp: Optional[Path] = None
        self._raw: dict[str, Any] = {}

    def __enter__(self) -> "TrajectoryWriter":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    @property
    def n_frames(self) -> int:
        """Frames written so far."""
        return self._counts["t"]

    def __call__(self, chunk: Trajectory) -> None:
        self.write(chunk)

    def write(self, chunk: Trajectory) -> None:
        """Appends the frames and events of ``chunk``."""
        if self._closed:
            raise ValueError("TrajectoryWriter is closed")
        energies = chunk.kinetic is not None
        c = chunk.tension.shape[1]
        if self._n is None:
            self._open(chunk.n_particles, c, energies)
            if self._metadata is None:
                self._metadata = dict(chunk.metadata)
        elif (chunk.n_particles, energies) != (self._n, self._energies):
            raise ValueError(
                f"chunk has {chunk.n_particles} particles (energies={energies}); "
                f"this file has {self._n} (energies={self._energies})"
            )
        elif c != self._c and chunk.n_frames > 0:
            raise ValueError(f"chunk has {c} constraints; this file has {self._c}")
        for name in _FRAME_FIELDS + _EVENT_FIELDS:
            value = getattr(chunk, name)
            if value is not None and name in self._fields():
                self._append(name, np.ascontiguousarray(value, dtype=_FIELDS[name][0]))
        if chunk.terminated_by is not None:
            self._terminated_by = chunk.terminated_by

    def close(self) -> None:
        """Finishes the file. Safe to call more than once."""
        if self._closed:
            return
        self._closed = True
        if self._n is None:
            self._open(0, 0, True)  # an empty but valid file
        if self._format == "hdf5":
            # A dataset rather than an attribute: attributes are limited to 64 KiB.
            self._h5.create_dataset("metadata", data=self._metadata_json())
            self._h5.attrs["terminated_by"] = -1 if self._terminated_by is None else self._terminated_by
            self._h5.close()
        else:
            self._close_npz()

    # -- internals ----------------------------------------------------------

    def _fields(self):
        skip = set() if self._energies else {"kinetic", "potential"}
        if self._c == 0:
            skip.add("tension")  # nothing to store, and HDF5 rejects zero-width chunks
        return [k for k in _FIELDS if k not in skip]

    def _metadata_json(self) -> str:
        return json.dumps(self._metadata or {})

    def _open(self, n: int, c: int, energies: bool) -> None:
        self._n, self._c, self._energies = n, c, energies
        if self._format == "hdf5":
            h5py = _h5py()
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._h5 = h5py.File(self.path, "w")
            for name in self._fields():
                dtype, row = _FIELDS[name]
                shape = row(n, c)
                self._h5.create_dataset(
                    name,
                    shape=(0,) + shape,
                    maxshape=(None,) + shape,
                    dtype=dtype,
                    # ~1 MiB chunks: few enough that the chunk index stays small.
                    chunks=(max(1, (1 << 17) // max(1, int(np.prod(shape)))),) + shape,
                )
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._tmp = Path(tempfile.mkdtemp(prefix=self.path.name + ".", dir=self.path.parent))
            self._raw = {name: open(self._tmp / name, "wb") for name in self._fields()}

    def _append(self, name: str, value: np.ndarray) -> None:
        rows = value.shape[0]
        if rows == 0:
            return
        if self._format == "hdf5":
            ds = self._h5[name]
            ds.resize(ds.shape[0] + rows, axis=0)
            ds[-rows:] = value
        else:
            value.tofile(self._raw[name])
        self._counts[name] += rows

    def _close_npz(self) -> None:
        # Copy each staged file into the archive block by block, so packing a file larger
        # than memory works and memory use stays flat.
        try:
            with zipfile.ZipFile(self.path, "w", zipfile.ZIP_STORED, allowZip64=True) as zf:
                for name, fh in self._raw.items():
                    fh.close()
                    dtype, row = _FIELDS[name]
                    header = {
                        "descr": np.lib.format.dtype_to_descr(np.dtype(dtype)),
                        "fortran_order": False,
                        "shape": (self._counts[name],) + row(self._n, self._c),
                    }
                    with zf.open(name + ".npy", "w", force_zip64=True) as out:
                        np.lib.format.write_array_header_1_0(out, header)
                        with open(self._tmp / name, "rb") as src:
                            shutil.copyfileobj(src, out, 1 << 20)
                small = {
                    "metadata": np.array(self._metadata_json()),
                    "terminated_by": np.int64(
                        -1 if self._terminated_by is None else self._terminated_by
                    ),
                }
                for name, value in small.items():
                    with zf.open(name + ".npy", "w") as out:
                        np.lib.format.write_array(out, value, allow_pickle=False)
        finally:
            shutil.rmtree(self._tmp, ignore_errors=True)


def save_trajectory(
    traj: Trajectory, path: PathLike, metadata: Optional[Mapping[str, Any]] = None
) -> None:
    """Writes a trajectory (frames, energies, events, metadata) to ``.npz`` or ``.h5``.

    ``metadata`` defaults to ``traj.metadata``.
    """
    with TrajectoryWriter(path, metadata if metadata is not None else traj.metadata) as w:
        w.write(traj)


def load_trajectory(path: PathLike) -> Trajectory:
    """Reads a file written by :func:`save_trajectory` or :class:`TrajectoryWriter`.

    The whole file is loaded into memory; for files too large for that, read the ``.h5``
    datasets directly with h5py (``f["pos"][start:stop]``).
    """
    if _format(path) == "hdf5":
        with _h5py().File(path, "r") as f:
            arrays = {k: f[k][()] for k in _FIELDS if k in f}
            metadata = json.loads(f["metadata"][()].decode())
            terminated_by = int(f.attrs["terminated_by"])
    else:
        with np.load(path, allow_pickle=False) as f:
            arrays = {k: f[k] for k in _FIELDS if k in f.files}
            metadata = json.loads(str(f["metadata"]))
            terminated_by = int(f["terminated_by"])
    return Trajectory(
        arrays["t"],
        arrays["pos"],
        arrays["vel"],
        arrays.get("kinetic"),
        arrays.get("potential"),
        tension=arrays.get("tension"),
        event_t=arrays["event_t"],
        event_index=[int(i) for i in arrays["event_index"]],
        event_pos=arrays["event_pos"],
        event_vel=arrays["event_vel"],
        terminated_by=None if terminated_by < 0 else terminated_by,
        metadata=metadata,
    )
