"""SLE-curve intersection -> predicted eutectic (T, x) for a binary A-B system.

Theory (Martins, Pinho & Coutinho 2018): the eutectic point is where the two
components' melting-point-depression (liquidus) curves cross. Each curve is
governed by the SLE equation:

    ln(x_i * gamma_i(x_i, T)) = -(dHfus_i / R) * (1/T - 1/Tfus_i)

gamma_i is supplied by a caller-provided function so this module has no direct
dependency on opencosmorspy/COSMORS -- see make_cosmors_gamma_fn() for the glue.
"""
import numpy as np
from scipy.optimize import brentq

R = 8.314  # J/(mol K)


def liquidus_x(idx, T, dHfus, Tfus, gamma_fn):
    """Mole fraction of component `idx` on its own liquidus curve at T, or
    None if no root exists in (0,1) at this T (T above the curve's range)."""
    rhs = -(dHfus / R) * (1.0 / T - 1.0 / Tfus)

    def resid(x_i):
        g = gamma_fn(idx, x_i, T)
        return np.log(x_i * g) - rhs

    try:
        return brentq(resid, 1e-6, 1 - 1e-6, xtol=1e-6)
    except ValueError:
        return None


def find_eutectic(dHfus_a, Tfus_a, dHfus_b, Tfus_b, gamma_fn, T_lo=250.0, T_hi=None, n_scan=60):
    """Scan T for the crossing of component A's and component B's liquidus
    curves (x_a(T) == 1 - x_b(T)); refine with brentq. Returns (T_eut, x_a_eut)
    or (None, None) if no sign change is found in [T_lo, T_hi]."""
    if T_hi is None:
        T_hi = min(Tfus_a, Tfus_b)

    def diff(T):
        xa = liquidus_x(0, T, dHfus_a, Tfus_a, gamma_fn)
        xb = liquidus_x(1, T, dHfus_b, Tfus_b, gamma_fn)
        if xa is None or xb is None:
            return None
        return xa - (1.0 - xb)

    Ts = np.linspace(T_lo, T_hi, n_scan)
    vals = [diff(T) for T in Ts]
    for i in range(len(Ts) - 1):
        if vals[i] is None or vals[i + 1] is None:
            continue
        if vals[i] * vals[i + 1] < 0:
            T_eut = brentq(lambda T: diff(T), Ts[i], Ts[i + 1], xtol=1e-3)
            x_eut = liquidus_x(0, T_eut, dHfus_a, Tfus_a, gamma_fn)
            return T_eut, x_eut
    return None, None


def make_cosmors_gamma_fn(crs):
    """Wrap a COSMORS instance (2 molecules already added via add_molecule,
    in [A, B] order) into the gamma_fn(idx, x_i, T) signature above."""
    def gamma_fn(idx, x_i, T):
        x = np.array([x_i, 1 - x_i]) if idx == 0 else np.array([1 - x_i, x_i])
        crs.clear_jobs()
        crs.add_job(x, T, refst="pure_component")
        res = crs.calculate()
        return float(np.exp(res["tot"]["lng"][0, idx]))
    return gamma_fn


def _self_check():
    """Ideal-solution (gamma=1) sanity check: no COSMO files needed. If the
    solver is broken (wrong sign, swapped indices, ...) the residuals below
    will not vanish at the returned point."""
    ideal_gamma_fn = lambda idx, x_i, T: 1.0
    dHfus_a, Tfus_a = 18000.0, 300.0  # J/mol, K -- arbitrary but distinct components
    dHfus_b, Tfus_b = 25000.0, 340.0

    T_eut, x_eut = find_eutectic(dHfus_a, Tfus_a, dHfus_b, Tfus_b, ideal_gamma_fn)
    assert T_eut is not None, "no eutectic crossing found for ideal test case"

    resid_a = np.log(x_eut * 1.0) - (-(dHfus_a / R) * (1.0 / T_eut - 1.0 / Tfus_a))
    resid_b = np.log((1 - x_eut) * 1.0) - (-(dHfus_b / R) * (1.0 / T_eut - 1.0 / Tfus_b))
    assert abs(resid_a) < 1e-3, f"component A SLE not satisfied at solution: {resid_a}"
    assert abs(resid_b) < 1e-3, f"component B SLE not satisfied at solution: {resid_b}"
    print(f"self-check OK: T_eut={T_eut:.2f} K, x_a_eut={x_eut:.4f}")


if __name__ == "__main__":
    _self_check()
