#!/usr/bin/env python3
"""
Identify the G1 peak DNA-content (PE-A) channel for each Parental/Resistant/
Fused flow cytometry replicate.

Background: the colleague's FlowJo workspaces (data_analysis_rep_1.wsp,
data_analysis_rep_2.wsp) only contain a single FSC-A/SSC-A scatter gate
named "cells" (debris exclusion) -- no G1/S/G2 gating or cell-cycle model
was ever applied. So instead of reading a G1 gate statistic out of FlowJo,
this script:
  1. parses each wsp to recover the "cells" polygon gate (dims + vertices)
     for the corresponding sample,
  2. applies that gate to the raw .fcs events (matplotlib Path.contains_points
     -- same polygon-in-polygon test FlowJo itself uses),
  3. finds the G1 and G2/M peaks in the gated PE-A histogram via
     scipy.signal.find_peaks, choosing the peak pair whose channel ratio is
     closest to the expected ~2x (G2/M cells have replicated DNA),
  4. refines the G1 peak location with a local Gaussian fit (more robust
     than a single histogram bin),
  5. reports per-replicate values and the per-condition (replicate-averaged)
     value ready to paste into notebooks/cnv_ploidy_postprocessing.ipynb.

Run with the project's isolated venv (created to avoid touching the base
conda env's numpy version, which fcsparser pins to <2):
    ./.venv-flow/bin/python scripts/analyze_flow_g1_peaks.py
"""

import xml.etree.ElementTree as ET
from pathlib import Path as FsPath
from collections import defaultdict

import numpy as np
import pandas as pd
import fcsparser
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

DATA_DIR = FsPath("data/Analysis_cell_cycle_170925")
WSP_FILES = [
    DATA_DIR / "FlowJo_Analysis" / "data_analysis_rep_1.wsp",
    DATA_DIR / "FlowJo_Analysis" / "data_analysis_rep_2.wsp",
]
FCS_DIR = DATA_DIR / "Raw_flow_cytometry_files"
OUT_DIR = FsPath("analysis/flow_cytometry")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DNA_CHANNEL = "PE-A"

SAMPLE_CONDITION = {
    "Par_1.fcs": "parental",
    "Par2.fcs": "parental",
    "Res1.fcs": "resistant",
    "Res2.fcs": "resistant",
    "Fus3-1b.fcs": "fused",
    "Fus4-1c.fcs": "fused",
}


def local(tag):
    return tag.split("}")[-1] if "}" in tag else tag


def strip_ns_key(k):
    return k.split("}")[-1] if "}" in k else k


def load_cells_gates(wsp_path):
    """Returns {fcs_filename: (dims, vertices)} for the 'cells' PolygonGate."""
    gates = {}
    tree = ET.parse(wsp_path)
    root = tree.getroot()
    for sample in root.iter():
        if local(sample.tag) != "Sample":
            continue
        sample_node = next(c for c in sample if local(c.tag) == "SampleNode")
        fcs_name = sample_node.attrib.get("name")

        for pop in sample.iter():
            if local(pop.tag) != "Population" or pop.attrib.get("name") != "cells":
                continue
            for gate in pop:
                if local(gate.tag) != "Gate":
                    continue
                for poly in gate:
                    if local(poly.tag) != "PolygonGate":
                        continue
                    dims, verts = [], []
                    for child in poly:
                        ctag = local(child.tag)
                        if ctag == "dimension":
                            for gc in child:
                                if local(gc.tag) == "fcs-dimension":
                                    attrib = {strip_ns_key(k): v for k, v in gc.attrib.items()}
                                    dims.append(attrib.get("name"))
                        if ctag == "vertex":
                            coords = []
                            for gc in child:
                                if local(gc.tag) == "coordinate":
                                    attrib = {strip_ns_key(k): v for k, v in gc.attrib.items()}
                                    coords.append(float(attrib.get("value")))
                            verts.append(coords)
                    gates[fcs_name] = (dims, verts)
    return gates


def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))


def refine_peak(pe_a_values, coarse_channel, window_frac=0.15, bins=200):
    window = coarse_channel * window_frac
    mask = (pe_a_values > coarse_channel - window) & (pe_a_values < coarse_channel + window)
    local_vals = pe_a_values[mask]
    if len(local_vals) < 20:
        return coarse_channel, None

    hist, edges = np.histogram(local_vals, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    try:
        popt, _ = curve_fit(
            gaussian, centers, hist,
            p0=[hist.max(), coarse_channel, window / 3],
            maxfev=5000,
        )
        return popt[1], popt  # refined mean
    except RuntimeError:
        return coarse_channel, None


def find_g1_g2_peaks(pe_a_values, min_channel=1000, n_bins=300, prominence=50):
    vals = pe_a_values[pe_a_values > min_channel]
    hist, edges = np.histogram(vals, bins=n_bins)
    centers = (edges[:-1] + edges[1:]) / 2

    peak_idx, props = find_peaks(hist, prominence=prominence, distance=5)
    if len(peak_idx) == 0:
        return None, None, hist, centers

    candidates = sorted(centers[peak_idx].tolist())

    best_pair, best_err = None, np.inf
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            ratio = candidates[j] / candidates[i]
            err = abs(ratio - 2.0)
            if err < best_err:
                best_err = err
                best_pair = (candidates[i], candidates[j])

    if best_pair is None or best_err > 0.5:
        # fall back to the single tallest peak as G1, no confident G2/M match
        tallest = centers[peak_idx[np.argmax(hist[peak_idx])]]
        return tallest, None, hist, centers

    return best_pair[0], best_pair[1], hist, centers


def analyze_sample(fcs_name, gates):
    fcs_path = FCS_DIR / fcs_name
    meta, data = fcsparser.parse(str(fcs_path), reformat_meta=True)

    dims, verts = gates[fcs_name]
    mplpath = MplPath(verts)
    mask = mplpath.contains_points(data[dims].values)
    gated = data[mask]

    pe_a = gated[DNA_CHANNEL].values
    coarse_g1, coarse_g2, hist, centers = find_g1_g2_peaks(pe_a)
    refined_g1, gfit = refine_peak(pe_a, coarse_g1)

    return {
        "fcs_name": fcs_name,
        "n_events_total": len(data),
        "n_events_gated": len(gated),
        "coarse_g1_channel": coarse_g1,
        "coarse_g2_channel": coarse_g2,
        "refined_g1_channel": refined_g1,
        "pe_a_gated": pe_a,
        "hist": hist,
        "centers": centers,
        "gauss_fit": gfit,
    }


def plot_diagnostic(result, condition):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(result["centers"], result["hist"], width=np.diff(result["centers"]).mean(), alpha=0.5)
    ax.axvline(result["coarse_g1_channel"], color="green", linestyle="--", label=f"G1 (coarse) {result['coarse_g1_channel']:.0f}")
    ax.axvline(result["refined_g1_channel"], color="darkgreen", linestyle="-", label=f"G1 (refined) {result['refined_g1_channel']:.0f}")
    if result["coarse_g2_channel"]:
        ax.axvline(result["coarse_g2_channel"], color="crimson", linestyle="--", label=f"G2/M {result['coarse_g2_channel']:.0f}")
    ax.set_title(f"{condition}: {result['fcs_name']}")
    ax.set_xlabel(f"{DNA_CHANNEL} (DNA content)")
    ax.set_ylabel("count")
    ax.legend(fontsize=8)
    plt.tight_layout()
    out_path = OUT_DIR / f"{FsPath(result['fcs_name']).stem}_g1_peak.png"
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    gates = {}
    for wsp in WSP_FILES:
        gates.update(load_cells_gates(wsp))

    rows = []
    for fcs_name, condition in SAMPLE_CONDITION.items():
        if fcs_name not in gates:
            print(f"[!] no 'cells' gate found for {fcs_name} in either wsp -- skipping")
            continue
        result = analyze_sample(fcs_name, gates)
        plot_path = plot_diagnostic(result, condition)

        ratio_ok = result["coarse_g2_channel"] is not None
        rows.append({
            "condition": condition,
            "fcs_name": fcs_name,
            "n_events_gated": result["n_events_gated"],
            "g1_channel_refined": result["refined_g1_channel"],
            "g2_channel": result["coarse_g2_channel"],
            "g2_g1_ratio": (result["coarse_g2_channel"] / result["refined_g1_channel"]) if ratio_ok else None,
            "diagnostic_plot": str(plot_path),
        })
        print(f"{condition:10s} {fcs_name:15s} G1={result['refined_g1_channel']:.0f}  "
              f"G2/M={result['coarse_g2_channel']}  plot={plot_path}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "g1_peak_per_replicate.csv", index=False)

    print("\n--- per-condition average (across replicates) ---")
    flow_g1_peak = {}
    for condition, group in df.groupby("condition"):
        avg = group["g1_channel_refined"].mean()
        flow_g1_peak[condition] = avg
        print(f"  {condition}: {group['g1_channel_refined'].tolist()} -> mean {avg:.1f}")

    print("\nPaste into notebooks/cnv_ploidy_postprocessing.ipynb:")
    print("FLOW_G1_PEAK = {")
    for name in ("parental", "resistant", "fused"):
        print(f'    "{name}": {flow_g1_peak.get(name):.1f},')
    print("}")

    df.to_csv(OUT_DIR / "g1_peak_summary.csv", index=False)


if __name__ == "__main__":
    main()
