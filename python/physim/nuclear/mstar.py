"""MSTAR: empirical electronic stopping of ions from Li to Ar, scaled from α particles (ICRU Report 73).

A port of the fitting functions of MSTAR 3.12 by H. Paul and A. Schinner (routines MSPAUL, MSTAR1 and MSID2Z;
freeware, IAEA Nuclear Data Section, https://www-nds.iaea.org/stopping/mstar). References: H. Paul and
A. Schinner, Nucl. Instr. Meth. B 179 (2001) 299 and B 195 (2002) 166; At. Data Nucl. Data Tables 85 (2003) 377.

The stopping of an ion Z₁ at energy E (MeV/u) is

    S(Z₁, E) = S_α(4E MeV) · F(Z₁, E) · Z₁² / 4,

with S_α from ICRU 49 (the ASTAR tables) and F a cumulative Weibull function of log₁₀ E fitted to measured data:
one set of coefficients for solids (mode "c", with a small target correction in mode "d"), one for gases (mode "g",
H₂ and He with their own coefficients) and one per ion for gases other than H₂ and He (mode "h"). :func:`mode`
picks them as MSTAR's mode "b" does.
"""

from __future__ import annotations

import math

import numpy as np

#: ICRU 49 target IDs that MSTAR treats as gases (MSTAR1, GASIDS).
GAS_IDS = frozenset({1, 2, 7, 8, 10, 18, 36, 54, 101, 104, 134, 155, 197, 238, 263, 264, 277})

#: Effective Z of the ICRU 49 compounds and mixtures, weighted from Table 1.2 of the report (MSID2Z, DZTAB).
COMPOUND_Z = {
    99: 5.942618, 101: 5.612910, 103: 5.905021, 104: 7.372747, 106: 10.646255, 111: 8.732800, 119: 9.247112,
    120: 10.628862, 126: 7.313059, 130: 14.646751, 134: 7.454168, 138: 7.131620, 139: 7.463094, 141: 54.023098,
    155: 5.281445, 160: 7.347248, 169: 10.322468, 179: 6.359930, 185: 7.394490, 189: 6.822535, 191: 6.461884,
    197: 4.743470, 200: 6.776222, 201: 7.129601, 202: 7.095221, 203: 6.964475, 204: 7.010650, 209: 5.916677,
    213: 5.256975, 215: 35.808854, 216: 5.575000, 219: 6.100061, 221: 5.281445, 222: 6.456261, 223: 6.236538,
    225: 5.281445, 226: 5.612910, 227: 8.279451, 232: 11.997960, 238: 5.085725, 245: 10.804610, 252: 46.558334,
    255: 5.664495, 263: 6.339387, 264: 6.108394, 266: 5.562450, 276: 7.216742, 277: 7.216742, 906: 6.0,
}

#: Mode "h" coefficients (a, b, c) for particular ions in gases other than H₂ and He.
_H_COEF = {3: (0.59251, -0.71211, 1.91956), 4: (0.35030, -0.28291, 2.92675), 5: (0.32192, -0.13151, 2.97854),
           6: (0.27061, -0.01952, 2.80761), 7: (0.21054, 0.01528, 3.99705), 8: (0.16969, 0.1075, 3.75755),
           9: (0.12192, 0.08955, 3.59459), 10: (0.09358, 0.08692, 3.39734), 11: (0.06353, 0.21303, 4.29082),
           16: (0.05698, 0.34491, 3.01333), 17: (0.05163, 0.37948, 2.81286), 18: (0.048, 0.38887, 2.84076)}

Z1_MIN, Z1_MAX = 3, 18


def target_z(icru_id: int) -> float:
    """The (effective) Z MSTAR uses for an ICRU 49 target ID."""
    return float(icru_id) if icru_id < 99 else COMPOUND_Z[icru_id]


def _d_supported(z2: int) -> bool:
    return 4 < z2 < 10 or 10 < z2 < 19 or 19 < z2 < 36 or 36 < z2 < 62 or 62 < z2 < 93


def mode(z1: int, icru_id: int) -> str:
    """MSTAR's mode "b": the particular fit ("d" for solids, "h" for gases) where there is one, else "c" or "g"."""
    z2 = int(target_z(icru_id) + 0.5)
    if icru_id in GAS_IDS:
        return "h" if z2 >= 3 and z1 in _H_COEF else "g"
    return "d" if _d_supported(z2) else "c"


def factor(z1: int, e_per_u, mode_: str, z2: float = 1.0) -> np.ndarray:
    """F(Z₁, E): the ion's stopping over Z₁²/4 times the α stopping at the same velocity (MSPAUL)."""
    if not Z1_MIN <= z1 <= Z1_MAX:
        raise ValueError(f"MSTAR covers Z1 = {Z1_MIN}..{Z1_MAX}, not {z1}")
    x, lx = float(z1), math.log(z1)
    iz2 = int(z2 + 0.5)
    if mode_ in ("c", "d"):
        a = (152.3281067 - 509.175005 * x**0.5 + 764.0597161 * x - 666.456225 * x**1.5 + 369.6135044 * x**2
             - 134.506826 * x**2.5 + 32.02328376 * x**3 - 4.80247840 * x**3.5 + 0.411498619 * x**4
             - 0.015353889 * x**4.5)
        b = ((-1.60299361 + 2.594689067 * x**0.5 - 1.73972608 * x + 0.601973528 * x**1.5 - 0.10683182 * x**2
              + 0.0077231543 * x**2.5)
             / (1.0 - 1.44605189 * x**0.5 + 0.871321669 * x - 0.24591206 * x**1.5 + 0.027726119 * x**2
                - 0.00042121782 * x**2.5))
        c = 1.0 / (0.15607719 - 0.0023161329 * x**2 * lx + 0.0016853943 * x**2.5)
        d1 = 17.18567 - 0.65734 * x
    elif mode_ == "g" and iz2 > 2:
        a = (38.150347 - 5.0998443 * x - 151.25405 / x + 0.3682134 * x**2 + 301.77474 / x**2
             - 1.3542411e-2 * x**3 - 225.08009 / x**3 + 1.9888596e-4 * x**4)
        b = -0.1983923 + 3.2653672e-2 * x - 12.281832 * math.exp(-x)
        c = math.exp(-1.3058678 - 0.5460869 * lx**2 + 2.3894315 * lx)
        d1 = 9.1963 - 0.29643 * x
    elif mode_ == "g":  # H₂ and He
        a = (82.333704 - 1.7353839e2 * x**0.5 + 1.4994474e2 * x - 67.608335 * x**1.5 + 1.6771605e1 * x**2
             - 2.1732999 * x**2.5 + 1.1514042e-1 * x**3)
        b = 0.74789173 - 1.6701292 / lx
        c = 3.2370042 - 0.049755185 * x
        d1 = 9.1963 - 0.29643 * x
    elif mode_ == "h" and iz2 >= 3 and z1 in _H_COEF:
        a, b, c = _H_COEF[z1]
        d1 = 9.1963 - 0.29643 * x
    else:
        raise ValueError(f"no MSTAR mode {mode_!r} for Z1 = {z1}, Z2 = {z2}")
    lg = np.log10(np.asarray(e_per_u, dtype=float))
    tp = np.maximum((lg + c * math.log(2.0) ** (1 / d1) - b) / c, 0.0)
    f = a + (1.01 - a) * (1.0 - np.exp(-tp**d1))
    if mode_ not in ("c", "d"):
        return f
    f = f * (1.0 - 0.0041466419 - 0.0010098198 * np.exp(-lg))
    if mode_ == "d":
        if 4 < iz2 < 10:  # C and Mylar
            f = f * (1.0 + 0.0012171149 - 0.04021975 / (1.0 + ((lg + 2.0388341) / 0.48193661) ** 2))
        elif 10 < iz2 < 19:  # Al, Si
            f = f * (1.0 + 0.003491232 - 0.11423452 / (1.0 + ((lg + 1.5906109) / 0.41741941) ** 2))
        elif 36 < iz2 < 62:  # Mo, Ag, Sn
            f = f * (1.0 - 0.0030928609 + 0.11488954 / (1.0 + ((lg + 1.4359953) / 0.37177939) ** 2))
        elif 62 < iz2 < 93:  # W, Pt, Au, Pb, U
            a0, b0, c0, d0 = -0.069155044, -0.65644678, 0.18379076, -0.0010104478
            n0 = d0**2 + c0**2
            erf = np.vectorize(math.erf)(d0 * (lg - b0) / (math.sqrt(2.0 * n0) * c0))
            f = f * (1.0 + a0 * c0 * np.exp(-0.5 * (lg - b0) ** 2 / n0) * (1.0 + erf) / math.sqrt(n0))
        elif not 19 < iz2 < 36:
            raise ValueError(f"no MSTAR mode 'd' for Z2 = {z2}")
    return f
