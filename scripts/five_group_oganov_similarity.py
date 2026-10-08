#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Five-group similarity analysis for generated carbon structures
using the Oganov–Valle fingerprint cosine distance.

Reference:
A. R. Oganov and M. Valle,
How to quantify energy landscapes of solids,
J. Chem. Phys. 130, 104504 (2009).
DOI: 10.1063/1.3079326

Pure-carbon fingerprint:
    F(R) = g(R) - 1

Structural distance:
    D = 1/2 * [1 - (F1 · F2) / (|F1||F2|)]

Interpretation:
    D = 0     -> identical fingerprint
    larger D  -> greater structural difference

For intuitive plotting:
    S = 1 - D
so larger S means more similar.

Input files:
    resultsC0000.csv
    resultsC1000.csv
    ...
    resultsC9000.csv

Grouping:
    Group 1 = resultsC0000 + resultsC1000
    Group 2 = resultsC2000 + resultsC3000
    Group 3 = resultsC4000 + resultsC5000
    Group 4 = resultsC6000 + resultsC7000
    Group 5 = resultsC8000 + resultsC9000

Default parameters from Oganov & Valle examples:
    Rmax  = 15.0 Å
    dR    = 0.05 Å
    sigma = 0.075 Å

Dependencies:
    pip install numpy pandas pymatgen scipy matplotlib

Usage:
    python five_group_oganov_similarity.py .
"""

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from pymatgen.io.vasp.inputs import Poscar

EXPECTED_RUNS = list(range(0, 10000, 1000))
GROUPS = [
    ("Group 1", [0, 1000]),
    ("Group 2", [2000, 3000]),
    ("Group 3", [4000, 5000]),
    ("Group 4", [6000, 7000]),
    ("Group 5", [8000, 9000]),
]

def find_generation_files(base: Path):
    found = {}
    for p in base.glob("resultsC*.csv"):
        m = re.search(r"resultsC\s*(\d{4})", p.name, flags=re.I)
        if m:
            found[int(m.group(1))] = p
    missing = [i for i in EXPECTED_RUNS if i not in found]
    if missing:
        raise FileNotFoundError(
            "Missing generation files: "
            + ", ".join(f"resultsC{i:04d}.csv" for i in missing)
        )
    return found

def detect_poscar_column(df: pd.DataFrame, requested=None):
    df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
    lower_map = {c.lower(): c for c in df.columns}
    if requested:
        key = requested.strip().lower()
        if key in lower_map:
            return lower_map[key]
    preferred = [
        "poscar", "structure", "vasp", "contcar",
        "poscar_text", "generated_poscar", "structure_poscar"
    ]
    for key in preferred:
        if key in lower_map:
            return lower_map[key]
    raise ValueError(
        "Could not detect POSCAR column automatically. "
        f"Available columns: {list(df.columns)}"
    )

def parse_poscar_text(value):
    text = str(value)
    if "\\n" in text:
        text = text.replace("\\n", "\n")
    return Poscar.from_str(text).structure

def is_pure_carbon(structure):
    return {el.symbol for el in structure.composition.elements} == {"C"}

def ordered_pair_distances(structure, rmax: float):
    distances = []
    for site in structure:
        for neigh in structure.get_neighbors(site, rmax):
            d = float(neigh.nn_distance)
            if d > 1e-10:
                distances.append(d)
    return np.asarray(distances, dtype=np.float64)

def oganov_fingerprint(structure, rmax=15.0, dr=0.05, sigma=0.075):
    n_atoms = len(structure)
    volume = float(structure.volume)
    if n_atoms <= 0 or volume <= 0:
        raise ValueError("Invalid structure: non-positive atom count or volume.")

    r = np.arange(0.0, rmax + 0.5 * dr, dr, dtype=np.float64)
    distances = ordered_pair_distances(structure, rmax)
    if len(distances) == 0:
        raise ValueError("No periodic pair distances found within Rmax.")

    pair_density_norm = (n_atoms ** 2) / volume
    g = np.zeros_like(r, dtype=np.float64)

    chunk_size = 2000
    for start in range(0, len(distances), chunk_size):
        d = distances[start:start + chunk_size]
        diff = r[:, None] - d[None, :]
        gaussian = np.exp(
            -0.5 * (diff / sigma) ** 2
        ) / (np.sqrt(2.0 * np.pi) * sigma)
        weights = 1.0 / (
            4.0 * np.pi * (d ** 2) * pair_density_norm
        )
        g += gaussian @ weights

    F = g - 1.0
    norm = np.linalg.norm(F)
    if norm <= 1e-15:
        raise ValueError("Zero-norm Oganov fingerprint.")

    return (F / norm).astype(np.float32)

def load_run_fingerprints(csv_path: Path, rmax, dr, sigma, poscar_column=None):
    df = pd.read_csv(csv_path)
    col = detect_poscar_column(df, poscar_column)

    fingerprints = []
    invalid = []

    for idx, value in df[col].items():
        try:
            structure = parse_poscar_text(value)
            if not is_pure_carbon(structure):
                invalid.append({
                    "source_file": csv_path.name,
                    "row": int(idx),
                    "reason": f"not pure carbon: {structure.composition.reduced_formula}",
                })
                continue
            fp = oganov_fingerprint(
                structure,
                rmax=rmax,
                dr=dr,
                sigma=sigma,
            )
            fingerprints.append(fp)
        except Exception as exc:
            invalid.append({
                "source_file": csv_path.name,
                "row": int(idx),
                "reason": str(exc),
            })

    if not fingerprints:
        raise RuntimeError(f"No valid pure-carbon structures found in {csv_path.name}")

    return np.vstack(fingerprints), invalid

def pairwise_oganov_distance(normalized_fingerprints):
    cosine_matrix = np.clip(
        normalized_fingerprints @ normalized_fingerprints.T,
        -1.0,
        1.0,
    )
    distance_matrix = 0.5 * (1.0 - cosine_matrix)
    iu = np.triu_indices(distance_matrix.shape[0], k=1)
    distances = distance_matrix[iu].astype(np.float64)
    similarities = 1.0 - distances
    return distances, similarities

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", nargs="?", default=".")
    parser.add_argument("--poscar-column", default=None)
    parser.add_argument("--rmax", type=float, default=15.0)
    parser.add_argument("--dr", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.075)
    parser.add_argument("--bin-width", type=float, default=0.025)
    args = parser.parse_args()

    base = Path(args.directory).resolve()
    files = find_generation_files(base)

    edges = np.arange(0.0, 1.0 + args.bin_width, args.bin_width, dtype=np.float64)
    edges[-1] = 1.0
    centers = (edges[:-1] + edges[1:]) / 2.0

    distance_distribution = pd.DataFrame({
        "distance_left": edges[:-1],
        "distance_right": edges[1:],
        "distance_bin_center": centers,
    })

    similarity_distribution = pd.DataFrame({
        "similarity_left": edges[:-1],
        "similarity_right": edges[1:],
        "similarity_bin_center": centers,
    })

    summary_rows = []
    distance_hist_probs = []
    invalid_all = []

    print("=" * 78)
    print("Five-group Oganov–Valle fingerprint analysis")
    print("D = 0.5 * [1 - cos(F1,F2)]")
    print("Smaller D = more similar")
    print("=" * 78)

    for group_name, run_ids in GROUPS:
        parts = []

        for run_id in run_ids:
            path = files[run_id]
            fps, invalid = load_run_fingerprints(
                path,
                rmax=args.rmax,
                dr=args.dr,
                sigma=args.sigma,
                poscar_column=args.poscar_column,
            )
            parts.append(fps)
            invalid_all.extend(invalid)
            print(f"{group_name} | {path.name}: valid pure-C = {fps.shape[0]}")

        group_fps = np.vstack(parts)
        distances, similarities = pairwise_oganov_distance(group_fps)

        d_counts, _ = np.histogram(distances, bins=edges)
        d_pct = d_counts / d_counts.sum() * 100.0
        d_prob = d_counts / d_counts.sum()
        distance_distribution[group_name.replace(" ", "_") + "_pct"] = d_pct
        distance_hist_probs.append(d_prob)

        s_counts, _ = np.histogram(similarities, bins=edges)
        s_pct = s_counts / s_counts.sum() * 100.0
        similarity_distribution[group_name.replace(" ", "_") + "_pct"] = s_pct

        summary_rows.append({
            "group": group_name,
            "run_1": f"resultsC{run_ids[0]:04d}.csv",
            "run_2": f"resultsC{run_ids[1]:04d}.csv",
            "valid_pure_C_structures": int(group_fps.shape[0]),
            "within_pair_count": int(len(distances)),
            "mean_oganov_distance_D": float(np.mean(distances)),
            "median_oganov_distance_D": float(np.median(distances)),
            "std_oganov_distance_D": float(np.std(distances, ddof=1)),
            "min_oganov_distance_D": float(np.min(distances)),
            "max_oganov_distance_D": float(np.max(distances)),
            "mean_similarity_1_minus_D": float(np.mean(similarities)),
            "median_similarity_1_minus_D": float(np.median(similarities)),
            "std_similarity_1_minus_D": float(np.std(similarities, ddof=1)),
        })

        print(
            f"  -> n={group_fps.shape[0]}, pairs={len(distances)}, "
            f"mean D={np.mean(distances):.6f}, "
            f"mean S={np.mean(similarities):.6f}"
        )

    group_cols = [f"Group_{i}_pct" for i in range(1, 6)]

    distance_distribution["five_group_mean_pct"] = (
        distance_distribution[group_cols].mean(axis=1)
    )
    distance_distribution["five_group_sd_pct"] = (
        distance_distribution[group_cols].std(axis=1, ddof=1)
    )

    similarity_distribution["five_group_mean_pct"] = (
        similarity_distribution[group_cols].mean(axis=1)
    )
    similarity_distribution["five_group_sd_pct"] = (
        similarity_distribution[group_cols].std(axis=1, ddof=1)
    )

    summary = pd.DataFrame(summary_rows)

    js_rows = []
    js_values = []

    for i in range(5):
        for j in range(i + 1, 5):
            js_div = float(
                jensenshannon(
                    distance_hist_probs[i],
                    distance_hist_probs[j],
                    base=2.0,
                ) ** 2
            )
            js_values.append(js_div)
            js_rows.append({
                "group_i": f"Group {i+1}",
                "group_j": f"Group {j+1}",
                "JS_divergence_distance_distribution": js_div,
            })

    mean_D = float(summary["mean_oganov_distance_D"].mean())
    sd_D = float(summary["mean_oganov_distance_D"].std(ddof=1))
    cv_D = sd_D / mean_D * 100.0 if mean_D != 0 else np.nan

    mean_S = float(summary["mean_similarity_1_minus_D"].mean())
    sd_S = float(summary["mean_similarity_1_minus_D"].std(ddof=1))
    cv_S = sd_S / mean_S * 100.0 if mean_S != 0 else np.nan

    overall = pd.DataFrame([{
        "mean_of_group_mean_D": mean_D,
        "SD_of_group_mean_D": sd_D,
        "CV_of_group_mean_D_percent": cv_D,
        "mean_of_group_mean_similarity": mean_S,
        "SD_of_group_mean_similarity": sd_S,
        "CV_of_group_mean_similarity_percent": cv_S,
        "mean_pairwise_JS_divergence": float(np.mean(js_values)),
        "max_pairwise_JS_divergence": float(np.max(js_values)),
    }])

    distance_distribution.to_csv(
        base / "01_oganov_distance_distribution_five_groups.csv",
        index=False, encoding="utf-8-sig", float_format="%.8f"
    )
    similarity_distribution.to_csv(
        base / "02_oganov_similarity_distribution_five_groups.csv",
        index=False, encoding="utf-8-sig", float_format="%.8f"
    )
    summary.to_csv(
        base / "03_oganov_five_group_summary.csv",
        index=False, encoding="utf-8-sig", float_format="%.8f"
    )
    pd.DataFrame(js_rows).to_csv(
        base / "04_oganov_five_group_JS_divergence.csv",
        index=False, encoding="utf-8-sig", float_format="%.8f"
    )
    overall.to_csv(
        base / "05_oganov_five_group_overall.csv",
        index=False, encoding="utf-8-sig", float_format="%.8f"
    )

    if invalid_all:
        pd.DataFrame(invalid_all).to_csv(
            base / "06_oganov_invalid_rows.csv",
            index=False, encoding="utf-8-sig"
        )

    fig, ax = plt.subplots(figsize=(8.5, 5.7))
    xD = distance_distribution["distance_bin_center"].to_numpy()
    for i in range(1, 6):
        ax.plot(
            xD,
            distance_distribution[f"Group_{i}_pct"].to_numpy(),
            marker="o", markersize=3, linewidth=1.5,
            label=f"Group {i}"
        )
    ax.set_xlabel("Oganov fingerprint cosine distance, D")
    ax.set_ylabel("Fraction of structure pairs (%)")
    ax.set_xlim(0.0, 1.0)
    ax.legend(frameon=True)
    ax.grid(alpha=0.20)
    fig.tight_layout()
    fig.savefig(base / "07_oganov_distance_distribution_five_groups.png",
                dpi=300, bbox_inches="tight")
    fig.savefig(base / "07_oganov_distance_distribution_five_groups.pdf",
                bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.7))
    xS = similarity_distribution["similarity_bin_center"].to_numpy()
    for i in range(1, 6):
        ax.plot(
            xS,
            similarity_distribution[f"Group_{i}_pct"].to_numpy(),
            marker="o", markersize=3, linewidth=1.5,
            label=f"Group {i}"
        )
    ax.set_xlabel("Structural fingerprint similarity, S = 1 - D")
    ax.set_ylabel("Fraction of structure pairs (%)")
    ax.set_xlim(0.0, 1.0)
    ax.legend(frameon=True)
    ax.grid(alpha=0.20)
    fig.tight_layout()
    fig.savefig(base / "08_oganov_similarity_distribution_five_groups.png",
                dpi=300, bbox_inches="tight")
    fig.savefig(base / "08_oganov_similarity_distribution_five_groups.pdf",
                bbox_inches="tight")
    plt.close(fig)

    method_text = f"""
Oganov–Valle five-group structural similarity analysis

Reference:
A. R. Oganov and M. Valle,
How to quantify energy landscapes of solids,
J. Chem. Phys. 130, 104504 (2009).
DOI: 10.1063/1.3079326

Groups:
Group 1 = resultsC0000 + resultsC1000
Group 2 = resultsC2000 + resultsC3000
Group 3 = resultsC4000 + resultsC5000
Group 4 = resultsC6000 + resultsC7000
Group 5 = resultsC8000 + resultsC9000

Fingerprint:
F(R) = g(R) - 1

Distance:
D = 1/2 * [1 - (F1 dot F2)/(|F1||F2|)]

Interpretation:
D = 0    : identical fingerprint
larger D : more dissimilar

Convenience similarity:
S = 1 - D
larger S : more similar

Parameters:
Rmax  = {args.rmax} Å
dR    = {args.dr} Å
sigma = {args.sigma} Å

Histogram normalization:
fraction (%) = number of within-group structure pairs in bin
               / total within-group pairs * 100
""".strip()

    (base / "09_oganov_method_settings.txt").write_text(
        method_text + "\n", encoding="utf-8"
    )

    print("\nOVERALL")
    print(overall.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    print("\nSaved outputs in:", base)

if __name__ == "__main__":
    main()
