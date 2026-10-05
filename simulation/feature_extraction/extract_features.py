"""Batch-extract ML-ready features from .orcacosmo files via opencosmorspy.

Output: one tidy CSV row per (hba, hbd, x_hba, T) combination, containing
  - per-molecule sigma-profile (10-bin) + sigma-moments (computed once per
    molecule and reused across every mixture/composition -- no extra DFT)
  - mixture activity coefficients (tot/enth/comb, both components)
  - H-bond vs. misfit interaction-energy decomposition (pm_E_hb, pm_E_mf)
  - excess Gibbs energy G^E/RT, derived as sum(x_i * ln(gamma_i))
  - predicted eutectic (T, x) and this row's distance from it (see eutectic.py)

Usage:
    python extract_features.py --pairs pairs.csv --out features.csv

pairs.csv columns: hba_name,hba_file,hbd_name,hbd_file[,dHfus_hba,Tfus_hba,dHfus_hbd,Tfus_hbd]
(fusion columns optional -- omit to skip eutectic-distance features for that pair)
"""
import argparse
import time

import numpy as np
import pandas as pd
from opencosmorspy import COSMORS
from opencosmorspy.input_parsers import SigmaProfileParser

from eutectic import find_eutectic, make_cosmors_gamma_fn

SIGMA_BIN_EDGES = np.arange(-0.025, 0.0251, 0.005)  # 11 edges -> 10 bins, Lemaoui 2020 / Mohan 2024 scheme
X_GRID = np.round(np.arange(0.05, 1.0, 0.05), 2)
T_GRID = [298.15, 318.15, 338.15, 358.15, 373.15]

_descriptor_cache = {}


def molecule_descriptors(filepath):
    """Per-molecule sigma-profile/moments, computed once per file and cached."""
    if filepath in _descriptor_cache:
        return _descriptor_cache[filepath]

    sp = SigmaProfileParser(filepath)
    sp.calculate_averaged_sigmas()
    _, sp_areas = sp.cluster_and_create_sigma_profile(
        sigmas="seg_sigma_averaged", sigmas_range=SIGMA_BIN_EDGES)
    sp.calculate_sigma_moments()

    row = {f"sigma_profile_bin{i}": a for i, a in enumerate(sp_areas)}
    row.update({f"sigma_moment{i}": m for i, m in enumerate(sp["sigma_moments"])})
    row["area"] = sp.get("area")
    row["volume"] = sp.get("volume")
    dip = sp.get("dipole_moment")
    row["dipole_moment_mag"] = float(np.linalg.norm(dip)) if dip is not None else None

    _descriptor_cache[filepath] = row
    return row


def mixture_rows(hba_name, hba_file, hbd_name, hbd_file, fusion=None, par="default_orca",
                 is_ionic_hba=False):
    crs = COSMORS(par=par)
    crs.par.calculate_contact_statistics_molecule_properties = True  # needed for pm_*/aim_*
    crs.add_molecule([hba_file])
    crs.add_molecule([hbd_file])

    eutectic_T, eutectic_x = (None, None)
    if fusion is not None:
        gamma_fn = make_cosmors_gamma_fn(crs)
        eutectic_T, eutectic_x = find_eutectic(
            fusion["dHfus_hba"], fusion["Tfus_hba"],
            fusion["dHfus_hbd"], fusion["Tfus_hbd"], gamma_fn)

    hba_desc = {f"hba_{k}": v for k, v in molecule_descriptors(hba_file).items()}
    hbd_desc = {f"hbd_{k}": v for k, v in molecule_descriptors(hbd_file).items()}

    rows = []
    for x_hba in X_GRID:
        for T in T_GRID:
            crs.clear_jobs()
            crs.add_job(np.array([x_hba, 1 - x_hba]), T, refst="pure_component")
            res = crs.calculate()

            lng_tot = res["tot"]["lng"][0]
            lng_enth = res["enth"]["lng"][0]
            lng_comb = res["comb"]["lng"][0]
            g_e_over_rt = float(x_hba * lng_tot[0] + (1 - x_hba) * lng_tot[1])

            row = {
                "hba_name": hba_name, "hbd_name": hbd_name,
                "x_hba": x_hba, "T": T, "is_ionic_hba": is_ionic_hba,
                "lng_hba_tot": lng_tot[0], "lng_hbd_tot": lng_tot[1],
                "lng_hba_enth": lng_enth[0], "lng_hbd_enth": lng_enth[1],
                "lng_hba_comb": lng_comb[0], "lng_hbd_comb": lng_comb[1],
                "gE_over_RT": g_e_over_rt,
            }
            for key in ("pm_E_hb", "pm_E_mf", "pm_A_int"):
                arr = res["enth"].get(key)
                if arr is not None:
                    row[f"hba_{key}"] = arr[0][0]
                    row[f"hbd_{key}"] = arr[0][1]

            row["T_eutectic_pred"] = eutectic_T
            row["x_eutectic_pred"] = eutectic_x
            row["delta_x_from_eutectic"] = (
                abs(x_hba - eutectic_x) if eutectic_x is not None else None)

            row.update(hba_desc)
            row.update(hbd_desc)
            rows.append(row)

    return rows


def mixture_rows_salt(salt_name, cation_name, cation_file, anion_name, anion_file,
                       hbd_name, hbd_file, par="default_orca"):
    """Salt (e.g. NaCl) + molecular HBD, treated as 3 separate COSMO species
    (cation, anion, HBD) with electroneutrality x_cation == x_anion enforced.
    Ionic HBAs cannot be geometry-optimized as one neutral molecule -- see
    Velez & Acevedo (2022) Type I/III DES classification and Lemaoui (2020),
    both of which treat the halide salt's cation/anion as separate species.

    No eutectic-distance feature here yet (binary SLE theory in eutectic.py
    does not extend to a 3-species ionic system) -- skipped, add if needed.
    """
    crs = COSMORS(par=par)
    crs.par.calculate_contact_statistics_molecule_properties = True
    crs.add_molecule([cation_file])  # idx 0
    crs.add_molecule([anion_file])   # idx 1
    crs.add_molecule([hbd_file])     # idx 2

    cation_desc = {f"hba_cation_{k}": v for k, v in molecule_descriptors(cation_file).items()}
    anion_desc = {f"hba_anion_{k}": v for k, v in molecule_descriptors(anion_file).items()}
    hbd_desc = {f"hbd_{k}": v for k, v in molecule_descriptors(hbd_file).items()}

    rows = []
    for x_salt in X_GRID:  # mole fraction of the NaCl formula unit
        x_ion = x_salt / 2.0
        x_hbd = 1 - x_salt
        for T in T_GRID:
            crs.clear_jobs()
            crs.add_job(np.array([x_ion, x_ion, x_hbd]), T, refst="pure_component")
            res = crs.calculate()

            lng_tot = res["tot"]["lng"][0]
            lng_enth = res["enth"]["lng"][0]
            lng_comb = res["comb"]["lng"][0]
            g_e_over_rt = float(x_ion * lng_tot[0] + x_ion * lng_tot[1] + x_hbd * lng_tot[2])

            row = {
                "hba_name": salt_name, "hbd_name": hbd_name,
                "x_hba": x_salt, "T": T, "is_salt_system": True,
                "lng_hba_cation_tot": lng_tot[0], "lng_hba_anion_tot": lng_tot[1], "lng_hbd_tot": lng_tot[2],
                "lng_hba_cation_enth": lng_enth[0], "lng_hba_anion_enth": lng_enth[1], "lng_hbd_enth": lng_enth[2],
                "lng_hba_cation_comb": lng_comb[0], "lng_hba_anion_comb": lng_comb[1], "lng_hbd_comb": lng_comb[2],
                "gE_over_RT": g_e_over_rt,
            }
            for key in ("pm_E_hb", "pm_E_mf", "pm_A_int"):
                arr = res["enth"].get(key)
                if arr is not None:
                    row[f"hba_cation_{key}"] = arr[0][0]
                    row[f"hba_anion_{key}"] = arr[0][1]
                    row[f"hbd_{key}"] = arr[0][2]

            row.update(cation_desc)
            row.update(anion_desc)
            row.update(hbd_desc)
            rows.append(row)

    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", help="CSV: hba_name,hba_file,hbd_name,hbd_file[,dHfus_hba,Tfus_hba,dHfus_hbd,Tfus_hbd]")
    ap.add_argument("--salt-pairs", help="CSV: salt_name,cation_name,cation_file,anion_name,anion_file,hbd_name,hbd_file")
    ap.add_argument("--out", required=True)
    ap.add_argument("--par", default="default_orca")
    args = ap.parse_args()

    n_pairs = (len(pd.read_csv(args.pairs)) if args.pairs else 0) + \
              (len(pd.read_csv(args.salt_pairs)) if args.salt_pairs else 0)
    done = 0
    wrote_header = False
    t_start = time.time()

    def flush(rows):
        nonlocal wrote_header
        if not rows:
            return
        pd.DataFrame(rows).to_csv(args.out, mode="a", header=not wrote_header, index=False)
        wrote_header = True

    def progress(label):
        nonlocal done
        done += 1
        elapsed = time.time() - t_start
        eta = elapsed / done * (n_pairs - done)
        print(f"[{done}/{n_pairs}] {label} done, elapsed={elapsed/60:.1f}min, eta={eta/60:.1f}min", flush=True)

    if args.pairs:
        pairs = pd.read_csv(args.pairs)
        fusion_cols = {"dHfus_hba", "Tfus_hba", "dHfus_hbd", "Tfus_hbd"}
        for _, p in pairs.iterrows():
            fusion = {c: p[c] for c in fusion_cols} if fusion_cols.issubset(pairs.columns) and p[list(fusion_cols)].notna().all() else None
            ionic = bool(p["is_ionic_hba"]) if "is_ionic_hba" in pairs.columns else False
            rows = mixture_rows(p["hba_name"], p["hba_file"], p["hbd_name"], p["hbd_file"], fusion, args.par, ionic)
            flush(rows)
            progress(f"{p['hba_name']}-{p['hbd_name']}")

    if args.salt_pairs:
        salt_pairs = pd.read_csv(args.salt_pairs)
        for _, p in salt_pairs.iterrows():
            rows = mixture_rows_salt(
                p["salt_name"], p["cation_name"], p["cation_file"],
                p["anion_name"], p["anion_file"], p["hbd_name"], p["hbd_file"], args.par)
            flush(rows)
            progress(p["salt_name"])

    print(f"done: {done} pairs written to {args.out}")


if __name__ == "__main__":
    main()
