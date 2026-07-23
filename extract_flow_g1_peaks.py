#!/usr/bin/env python3
"""
Extract the G1-gate DNA-channel peak value per sample from a FlowJo
workspace (.wsp) + its .fcs files, for filling in FLOW_G1_PEAK in
notebooks/cnv_ploidy_postprocessing.ipynb.

Use this if scripts/inspect_flowjo_workspace.py showed manual gating (no
cell-cycle model). If a Watson/Dean-Jett-Fox model was fit instead, skip
this script -- just read the fitted G1 mean/peak channel directly from
FlowJo's statistics view for each sample, that number is already better
than what this script computes.

Requires: pip install flowkit
(flowkit re-applies the exact gates/compensation/transforms from the .wsp
to the raw .fcs events -- this is the same gating your colleague did in
the FlowJo GUI, just reproduced in Python.)

Fill in GATE_NAME, DNA_CHANNEL, and SAMPLE_NAME_MAP below using the output
of inspect_flowjo_workspace.py, then run:
    python extract_flow_g1_peaks.py /path/to/experiment.wsp /path/to/fcs_dir/
"""

import sys
from pathlib import Path

import flowkit as fk

# --- fill these in using scripts/inspect_flowjo_workspace.py output ---

GATE_NAME = "G1"          # exact gate/population name for the G1 peak in the wsp
DNA_CHANNEL = "PE-A"       # exact DNA-stain channel name (not FSC-A/SSC-A)

# map FlowJo's internal sample names (from the inspector's "Sample:" lines)
# to your notebook's sample keys
SAMPLE_NAME_MAP = {
    "FlowJo_sample_name_for_parental": "parental",
    "FlowJo_sample_name_for_resistant": "resistant",
    "FlowJo_sample_name_for_fused": "fused",
}


def main(wsp_path, fcs_dir):
    session = fk.Session(fcs_samples=fcs_dir)
    session.import_flowjo_workspace(wsp_path)

    g1_peaks = {}
    for sample_id in session.get_sample_ids():
        gated = session.get_gate_events(sample_id, gate_name=GATE_NAME)
        if gated is None or len(gated) == 0:
            print(f"  [!] no events found for sample={sample_id}, gate={GATE_NAME} -- check names")
            continue

        median_val = gated[DNA_CHANNEL].median()
        mapped_name = SAMPLE_NAME_MAP.get(sample_id, sample_id)
        g1_peaks[mapped_name] = median_val
        print(f"{sample_id} -> {mapped_name}: n_events={len(gated)}, "
              f"G1 median {DNA_CHANNEL}={median_val:.2f}")

    print("\nPaste this into the notebook's FLOW_G1_PEAK cell:")
    print("FLOW_G1_PEAK = {")
    for name in ("parental", "resistant", "fused"):
        val = g1_peaks.get(name)
        print(f'    "{name}": {val!r},')
    print("}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_flow_g1_peaks.py /path/to/experiment.wsp /path/to/fcs_dir/")
        sys.exit(1)
    main(sys.argv[1], Path(sys.argv[2]))
