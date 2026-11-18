# KSHELL 2019 release (vendored)

`kshell.tar.gz` is the KSHELL shell-model code as published with N. Shimizu, T. Mizusaki, Y. Utsuno and
Y. Tsunoda, Comput. Phys. Commun. 244, 372 (2019), from Mendeley Data `kgzdz6ryyk`. It is unchanged:
382 339 bytes, SHA-256 `ce0659c3b39dc8c8ed623f6ebe1e10766445fc72372b6e33e952a578c7227526` (pinned in
`scripts/build_kshell.py`, which checks it before every build).

KSHELL is free software under the GNU General Public License, version 3 or later (header of `src/kshell.F90`;
the GPL text is `installer/windows/licenses/GPL-3.0.txt`). It is a separate program from physim (MIT): physim
runs its executables and is not linked with it.

It is kept here for two reasons:

- **Reproducible builds.** Mendeley Data refuses automated downloads, so CI and `scripts/build_kshell.py` build from
  this copy.
- **The source ships with the binaries,** as the GPL asks: the build copies this archive into
  `tools/kshell/source/` with physim's patch (`physim-rename.patch`) and the build script.

physim's only change, applied at build time and never to this file, is the rename patch described in
`scripts/build_kshell.py` (POSIX `rename` behaviour on Windows for Lanczos vectors on disk).
