#!/usr/bin/env python3
"""
Peek inside a FlowJo .wsp workspace to see what's actually in it, before
deciding how to pull the G1 peak value out.

No FlowJo license or extra packages needed -- .wsp is XML, this uses only
the standard library. It doesn't apply gates or touch the .fcs data; it just
lists what's there so you know:
  - which sample name maps to which .fcs file (parental/resistant/fused)
  - what the DNA-stain channel/parameter is called (e.g. "PE-A", "DAPI-A")
  - what the gate/population names are (e.g. "G1", "S", "G2/M")
  - whether a cell-cycle model (Watson Pragmatic / Dean-Jett-Fox) was fit,
    which shows up as extra <Statistic> nodes with names like "%G1",
    "Dean-Jett-Fox", "RCS", "G1 Mean", etc.

Usage:
    python inspect_flowjo_workspace.py /path/to/experiment.wsp
"""

import sys
import xml.etree.ElementTree as ET
from collections import defaultdict


def local_tag(elem):
    return elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag


def walk(elem, samples, current_sample=None):
    tag = local_tag(elem)
    attrib = elem.attrib

    if tag == "Sample":
        name = attrib.get("name") or attrib.get("sampleID") or "unnamed_sample"
        current_sample = name
        samples[current_sample]  # touch to create entry

    if tag in ("SampleNode",):
        # SampleNode usually carries the .fcs filename reference
        name = attrib.get("name")
        if current_sample and name:
            samples[current_sample]["fcs_ref"] = name

    if tag in ("Population", "Gate", "PolygonGate", "RectangleGate", "EllipsoidGate"):
        name = attrib.get("name")
        if current_sample and name:
            samples[current_sample]["gates"].append(name)

    if tag in ("Parameter", "fcs-dimension"):
        name = attrib.get("name") or attrib.get("value")
        if current_sample and name:
            samples[current_sample]["channels"].add(name)

    if tag in ("Statistic",):
        name = attrib.get("name") or attrib.get("id") or str(attrib)
        if current_sample:
            samples[current_sample]["statistics"].append(name)

    # generic catch-all: anything with "cellcycle", "watson", "dean" etc in
    # the tag or attributes, in case FlowJo's plugin uses nonstandard tags
    blob = (tag + " " + " ".join(f"{k}={v}" for k, v in attrib.items())).lower()
    for keyword in ("cellcycle", "watson", "dean", "jett", "fox", "%g1", "g1 mean", "g1mean"):
        if keyword in blob and current_sample:
            samples[current_sample]["cell_cycle_hits"].append(f"<{tag}> {attrib}")

    for child in elem:
        walk(child, samples, current_sample)


def main(wsp_path):
    tree = ET.parse(wsp_path)
    root = tree.getroot()

    samples = defaultdict(lambda: {
        "fcs_ref": None,
        "gates": [],
        "channels": set(),
        "statistics": [],
        "cell_cycle_hits": [],
    })

    walk(root, samples)

    for name, info in samples.items():
        print(f"\n=== Sample: {name} ===")
        print(f"  fcs reference : {info['fcs_ref']}")
        print(f"  channels seen : {sorted(info['channels'])}")
        print(f"  gate names    : {info['gates']}")
        print(f"  statistics    : {info['statistics']}")
        if info["cell_cycle_hits"]:
            print("  ** cell-cycle model evidence found: **")
            for hit in info["cell_cycle_hits"]:
                print(f"     {hit}")
        else:
            print("  no cell-cycle model keywords found (Watson/Dean-Jett-Fox) "
                  "-- looks like manual gating only, based on tag names alone")

    print(
        "\nWhat to do with this:\n"
        "  - Match 'fcs reference' to your parental/resistant/fused files.\n"
        "  - The DNA-stain channel is whichever entry in 'channels seen' is your\n"
        "    PI/DAPI detector (not FSC-A/SSC-A).\n"
        "  - If 'cell-cycle model evidence found' printed for a sample, open that\n"
        "    sample's statistics view in FlowJo and read off its fitted G1 mean/\n"
        "    peak channel directly -- that's a better number than a raw gate median.\n"
        "  - If not, the gate named like 'G1'/'G0' in 'gate names' is a manual gate;\n"
        "    use scripts/extract_flow_g1_peaks.py to compute the median DNA-channel\n"
        "    value of events inside it.\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_flowjo_workspace.py /path/to/experiment.wsp")
        sys.exit(1)
    main(sys.argv[1])
