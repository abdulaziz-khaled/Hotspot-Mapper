#!/usr/bin/env python3
"""
ppi_hotspot_analysis.py
=======================
Hotspot identification at protein–protein interfaces from PDBePISA JSON output.

Usage
-----
    python ppi_hotspot_analysis.py --json <path_to_pisa_json> --outdir <output_directory>

    Optional arguments:
        --chain1_label  Display label for chain 1 (default: "Chain 1")
        --chain2_label  Display label for chain 2 (default: "Chain 2")
        --top_n         Number of top hotspot residues to display/plot (default: 10)

Examples
--------
    # Basic usage (TNF-α / TNFR1 case study):
    python ppi_hotspot_analysis.py --json data/interface_7.json --outdir results/

    # With custom chain labels:
    python ppi_hotspot_analysis.py \
        --json data/interface_7.json \
        --outdir results/ \
        --chain1_label "Chain F (TNF-α)" \
        --chain2_label "Chain L (TNFR1)"

Input
-----
    PDBePISA JSON export for a single interface. The JSON must contain three
    top-level keys:
        - "interface"            : global interface metrics (BSA, ΔiG, p-value)
        - "interfacing_residues" : per-residue data for both chains
        - "interfacing_bonds"    : explicit bond inventory (hbonds, salt_bridges)

    To export this JSON:
        1. Go to https://www.ebi.ac.uk/pdbe/pisa/
        2. Enter your PDB ID and click Analyse
        3. Click the Interfaces tab
        4. Identify the interface index for your chain pair of interest
        5. Export the JSON for that specific interface

Output
------
    CSV tables (saved to --outdir):
        - interface_residues_chain1.csv
        - interface_residues_chain2.csv
        - hotspot_scores_chain1.csv
        - hotspot_scores_chain2.csv
        - bond_inventory_hbonds.csv
        - bond_inventory_saltbridges.csv

    Figures (saved to figures/ subfolder of --outdir):
        - figure_hotspot_ranking.svg / .png
        - figure_BSA_vs_energy.svg / .png
        - figure_2D_interaction_map.svg / .png

Portability note
----------------
    This script is designed to work with any PDBePISA interface JSON, not just
    the TNF-α/TNFR1 case study. Before running on a new complex:

        1. Verify JSON keys:
               print(list(data["interfacing_residues"].keys()))
               # Expected: ['chain1', 'chain2']

        2. BSA normalization is data-driven by default (bsa_norm=None in
           compute_hotspot_scores). Do NOT hardcode BSA_NORM to a fixed value
           from a different complex — it will distort the hotspot ranking.

        3. Update --chain1_label and --chain2_label to match your complex.

        4. Ensure you exported the correct interface index from PDBePISA.
           Exporting the wrong interface produces no error but gives results
           for the wrong chain pair.

References
----------
    [1] Berman et al. (2000) Nucleic Acids Res. 28:235-242
        https://doi.org/10.1093/nar/28.1.235
    [2] Krissinel & Henrick (2007) J. Mol. Biol. 372:774-797
        https://doi.org/10.1016/j.jmb.2007.05.022
    [3] Krissinel (2015) Nucleic Acids Res. 43:W314-W319
        https://doi.org/10.1093/nar/gkv314
    [4] McKinney (2010) Proc. SciPy 2010, 51-56
        https://doi.org/10.25080/Majora-92bf1922-00a
    [5] Harris et al. (2020) Nature 585:357-362
        https://doi.org/10.1038/s41586-020-2649-2
"""

import json
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'


# =============================================================================
# Constants
# =============================================================================

# Bond-type bonus weights for the composite hotspot score.
# These reflect the relative energetic importance of each interaction type.
# H-bond only: +1.0 | Salt bridge only: +1.5 | Both: +2.5 | None: +0.0
BOND_BONUS = {
    "H":  1.0,
    "S":  1.5,
    "HS": 2.5,
    None: 0.0,
}

# Human-readable labels for PDBePISA bond type codes.
INTERACTION_LABEL = {
    "H":  "H-bond",
    "S":  "Salt bridge",
    "HS": "H-bond + Salt bridge",
    None: "Hydrophobic/VdW",
}

# Colorblind-friendly palette for interaction types.
COLOR_MAP = {
    "H-bond":                "#4C72B0",
    "Salt bridge":           "#DD8452",
    "H-bond + Salt bridge":  "#55A868",
    "Hydrophobic/VdW":       "#C0C0C0",
}


# =============================================================================
# Parsing functions
# =============================================================================

def load_pisa_json(filepath):
    """
    Load and return a parsed PDBePISA interface JSON file.

    Parameters
    ----------
    filepath : str
        Path to the PDBePISA JSON export file.

    Returns
    -------
    dict
        Parsed JSON as a Python dictionary.
    """
    with open(filepath, "r") as f:
        return json.load(f)


def extract_interface_residues(data, chain_key):
    """
    Extract interface residues (BSA > 0 Å²) for a given chain.

    PDBePISA defines interface residues as those with buried surface area
    (BSA) > 0 upon complex formation. Residues with BSA = 0 are surface-
    exposed but non-interfacing and are excluded.

    Parameters
    ----------
    data : dict
        Parsed PDBePISA JSON (top-level dictionary).
    chain_key : str
        Key identifying the chain in 'interfacing_residues'.
        Typically 'chain1' or 'chain2'. Verify with:
            print(list(data["interfacing_residues"].keys()))

    Returns
    -------
    pd.DataFrame
        DataFrame of interface residues with columns:
            residue, seq_id, bsa, solv_en, bond_type, interaction
        Sorted by BSA descending.
    """
    residues = data["interfacing_residues"][chain_key]
    records = []
    for res in residues:
        bsa = float(res.get("bsa", 0))
        if bsa > 0:
            bt = res.get("bond_type", None)
            records.append({
                "residue":     res["residue_name"],
                "seq_id":      int(res["auth_seq_id"]),
                "bsa":         bsa,
                "solv_en":     float(res.get("solv_en", 0)),
                "bond_type":   bt,
                "interaction": INTERACTION_LABEL.get(bt, "Hydrophobic/VdW"),
            })
    df = pd.DataFrame(records).sort_values("bsa", ascending=False).reset_index(drop=True)
    return df


def parse_bonds(data, bond_category):
    """
    Parse the explicit bond inventory from PDBePISA JSON output.

    PDBePISA reports each cross-chain bond in both directions (chain1→chain2
    and chain2→chain1). Use get_unique_bonds() to deduplicate.

    Parameters
    ----------
    data : dict
        Parsed PDBePISA JSON.
    bond_category : str
        'hbonds' or 'salt_bridges'. Verify available keys with:
            print(list(data["interfacing_bonds"].keys()))

    Returns
    -------
    pd.DataFrame
        Bond table with columns:
            chain1_res, chain1_seqid, chain1_atom,
            chain2_res, chain2_seqid, chain2_atom, distance
        Sorted by distance ascending.
    """
    bonds_raw = data["interfacing_bonds"].get(bond_category, [])
    records = []
    for bond in bonds_raw:
        records.append({
            "chain1_res":   bond["chain1_residue_name"],
            "chain1_seqid": int(bond["chain1_auth_seq_id"]),
            "chain1_atom":  bond["chain1_atom"],
            "chain2_res":   bond["chain2_residue_name"],
            "chain2_seqid": int(bond["chain2_auth_seq_id"]),
            "chain2_atom":  bond["chain2_atom"],
            "distance":     float(bond["dist"]),
        })
    df = pd.DataFrame(records).sort_values("distance").reset_index(drop=True)
    return df


def get_unique_bonds(df_bonds):
    """
    Deduplicate a PDBePISA bond table.

    PDBePISA lists each cross-chain bond bidirectionally. This function
    retains only one canonical entry per unique bond pair (identified by
    the frozenset of both atom endpoints).

    Parameters
    ----------
    df_bonds : pd.DataFrame
        Raw bond DataFrame from parse_bonds().

    Returns
    -------
    pd.DataFrame
        Deduplicated bond DataFrame.
    """
    seen = set()
    unique = []
    for _, row in df_bonds.iterrows():
        key = frozenset([
            (row["chain1_seqid"], row["chain1_atom"]),
            (row["chain2_seqid"], row["chain2_atom"]),
        ])
        if key not in seen:
            seen.add(key)
            unique.append(row)
    return pd.DataFrame(unique).reset_index(drop=True)


# =============================================================================
# Hotspot scoring
# =============================================================================

def compute_hotspot_scores(df, bsa_norm=None):
    """
    Compute composite hotspot scores for interface residues.

    The composite score integrates three complementary descriptors:

        Score = (BSA / bsa_norm) + |ΔG_solv| + bond_bonus

    Components:
        BSA / bsa_norm  : Normalized buried surface area. Captures the
                          geometric contribution — residues that bury more
                          surface area contribute more to shape complementarity.
                          bsa_norm defaults to the maximum BSA in the dataset,
                          making this component data-driven and portable across
                          different complexes.

        |ΔG_solv|       : Absolute per-residue solvation energy contribution
                          (kcal/mol). Captures the thermodynamic contribution.
                          Both favorable (negative ΔG, hydrophobic burial) and
                          unfavorable (positive ΔG, polar desolvation) extremes
                          indicate energetically important residues.

        bond_bonus      : Rewards residues engaged in specific polar contacts.
                          H-bond: +1.0 | Salt bridge: +1.5 | Both: +2.5
                          Purely hydrophobic residues receive no bonus (+0.0).

    IMPORTANT — BSA normalization:
        bsa_norm defaults to None, which triggers automatic computation as
        df["bsa"].max(). Do NOT hardcode this to a value from a different
        complex — it will distort the ranking for your dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Interface residue DataFrame from extract_interface_residues().
    bsa_norm : float or None
        BSA normalization factor (Å²). If None, computed as df["bsa"].max().

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added columns: bond_bonus, hotspot_score, rank.
        Sorted by hotspot_score descending.
    """
    df = df.copy()
    if bsa_norm is None:
        bsa_norm = df["bsa"].max()   # Data-driven: adapts to any complex
    df["bond_bonus"]    = df["bond_type"].map(BOND_BONUS).fillna(0.0)
    df["hotspot_score"] = (
        df["bsa"] / bsa_norm
        + df["solv_en"].abs()
        + df["bond_bonus"]
    )
    df = df.sort_values("hotspot_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


# =============================================================================
# Visualization
# =============================================================================

def plot_hotspot_ranking(df_1, df_2, fig_dir, chain1_label, chain2_label, top_n=10):
    """
    Horizontal bar chart of top hotspot residues for both chains.

    Parameters
    ----------
    df_1, df_2 : pd.DataFrame
        Hotspot-scored DataFrames for chain 1 and chain 2.
    fig_dir : str
        Output directory for figures.
    chain1_label, chain2_label : str
        Display labels for the two chains.
    top_n : int
        Number of top residues to display.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, df_scored, label in zip(
        axes, [df_1, df_2], [chain1_label, chain2_label]
    ):
        top = df_scored.head(top_n)
        labels = [f"{r['residue']}{r['seq_id']}" for _, r in top.iterrows()]
        colors = [COLOR_MAP[i] for i in top["interaction"]]
        ax.barh(labels[::-1], top["hotspot_score"].values[::-1],
                color=colors[::-1], edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Composite Hotspot Score", fontsize=11)
        ax.set_title(label, fontsize=12, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=9)

    patches = [mpatches.Patch(color=v, label=k) for k, v in COLOR_MAP.items()]
    fig.legend(handles=patches, loc="lower center", ncol=4,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.04))
    plt.suptitle("Interface Residue Hotspot Ranking",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    for ext in ("svg", "png"):
        plt.savefig(os.path.join(fig_dir, f"figure_hotspot_ranking.{ext}"),
                    bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close()
    print("  Saved: figure_hotspot_ranking.svg / .png")


def plot_bsa_vs_energy(df_1, df_2, fig_dir, chain1_label, chain2_label, top_label=5):
    """
    Scatter plot of BSA vs. solvation energy for both chains.

    Parameters
    ----------
    df_1, df_2 : pd.DataFrame
        Hotspot-scored DataFrames for chain 1 and chain 2.
    fig_dir : str
        Output directory for figures.
    chain1_label, chain2_label : str
        Display labels for the two chains.
    top_label : int
        Number of top hotspot residues to annotate on the plot.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, df_scored, label in zip(
        axes, [df_1, df_2], [chain1_label, chain2_label]
    ):
        for interaction, group in df_scored.groupby("interaction"):
            ax.scatter(group["bsa"], group["solv_en"],
                       label=interaction, color=COLOR_MAP[interaction],
                       s=60, alpha=0.85, edgecolors="white", linewidths=0.5)
        for _, row in df_scored.head(top_label).iterrows():
            ax.annotate(f"{row['residue']}{row['seq_id']}",
                        (row["bsa"], row["solv_en"]),
                        fontsize=8, ha="left", va="bottom",
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("Buried Surface Area (Å²)", fontsize=11)
        ax.set_ylabel("Solvation Energy (kcal/mol)", fontsize=11)
        ax.set_title(label, fontsize=12, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.tick_params(axis="both", labelsize=9)

    axes[0].legend(fontsize=8, frameon=False)
    plt.suptitle("BSA vs. Solvation Energy at the PPI Interface",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    for ext in ("svg", "png"):
        plt.savefig(os.path.join(fig_dir, f"figure_BSA_vs_energy.{ext}"),
                    bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close()
    print("  Saved: figure_BSA_vs_energy.svg / .png")


def plot_2d_interaction_map(df_hbonds_unique, df_sb_unique, fig_dir,
                             chain1_label, chain2_label):
    """
    2D schematic interaction map of cross-chain bonds.

    Chain 1 residues are drawn on the left; chain 2 residues on the right.
    Hydrogen bonds are shown as dashed blue lines; salt bridges as solid red lines.
    Residue nodes are colored by physicochemical property.

    Parameters
    ----------
    df_hbonds_unique : pd.DataFrame
        Deduplicated hydrogen bond DataFrame.
    df_sb_unique : pd.DataFrame
        Deduplicated salt bridge DataFrame.
    fig_dir : str
        Output directory for figures.
    chain1_label, chain2_label : str
        Display labels for the two chains.
    """
    # Collect unique residues from both bond tables
    c1_residues, c2_residues = [], []
    for _, row in pd.concat([df_hbonds_unique, df_sb_unique]).iterrows():
        lbl1 = f"{row['chain1_res']}{row['chain1_seqid']}"
        lbl2 = f"{row['chain2_res']}{row['chain2_seqid']}"
        if lbl1 not in c1_residues:
            c1_residues.append(lbl1)
        if lbl2 not in c2_residues:
            c2_residues.append(lbl2)

    n1 = len(c1_residues)
    n2 = len(c2_residues)
    height = max(n1, n2) * 1.2 + 2

    fig, ax = plt.subplots(figsize=(9, height))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, height)
    ax.axis("off")

    def node_color(res_label):
        aa = res_label[:3]
        if aa in {"LYS", "ARG", "HIS"}: return "#4C72B0"
        if aa in {"ASP", "GLU"}:        return "#C44E52"
        if aa in {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}: return "#55A868"
        if aa == "GLY":                 return "#AAAAAA"
        return "#DDDDDD"

    y1 = {res: height - 1.5 - i * (height - 2) / max(n1 - 1, 1)
          for i, res in enumerate(c1_residues)}
    y2 = {res: height - 1.5 - i * (height - 2) / max(n2 - 1, 1)
          for i, res in enumerate(c2_residues)}

    # Draw residue nodes
    for res, y in y1.items():
        ax.add_patch(plt.Circle((2, y), 0.35, color=node_color(res), zorder=3))
        ax.text(1.5, y, res, ha="right", va="center", fontsize=9, fontweight="bold")
    for res, y in y2.items():
        ax.add_patch(plt.Circle((8, y), 0.35, color=node_color(res), zorder=3))
        ax.text(8.5, y, res, ha="left", va="center", fontsize=9, fontweight="bold")

    # Draw hydrogen bonds (dashed blue)
    for _, row in df_hbonds_unique.iterrows():
        l1 = f"{row['chain1_res']}{row['chain1_seqid']}"
        l2 = f"{row['chain2_res']}{row['chain2_seqid']}"
        if l1 in y1 and l2 in y2:
            ax.plot([2.35, 7.65], [y1[l1], y2[l2]],
                    color="#4C72B0", linewidth=1.5, linestyle="--", alpha=0.8, zorder=2)
            mx, my = 5.0, (y1[l1] + y2[l2]) / 2
            ax.text(mx, my + 0.15, f"{row['distance']:.3f} Å",
                    ha="center", va="bottom", fontsize=7.5, color="#4C72B0")

    # Draw salt bridges (solid red)
    for _, row in df_sb_unique.iterrows():
        l1 = f"{row['chain1_res']}{row['chain1_seqid']}"
        l2 = f"{row['chain2_res']}{row['chain2_seqid']}"
        if l1 in y1 and l2 in y2:
            ax.plot([2.35, 7.65], [y1[l1], y2[l2]],
                    color="#C44E52", linewidth=2.0, linestyle="-", alpha=0.85, zorder=2)
            mx, my = 5.0, (y1[l1] + y2[l2]) / 2
            ax.text(mx, my - 0.15, f"{row['distance']:.3f} Å",
                    ha="center", va="top", fontsize=7.5, color="#C44E52")

    # Column headers
    ax.text(2, height - 0.5, chain1_label,
            ha="center", va="center", fontsize=10, fontweight="bold")
    ax.text(8, height - 0.5, chain2_label,
            ha="center", va="center", fontsize=10, fontweight="bold")

    # Legend
    legend_elements = [
        mpatches.Patch(color="#4C72B0", label="Positively charged (K/R/H)"),
        mpatches.Patch(color="#C44E52", label="Negatively charged (D/E)"),
        mpatches.Patch(color="#55A868", label="Polar (S/T/N/Q/Y/C)"),
        mpatches.Patch(color="#AAAAAA", label="Glycine"),
        plt.Line2D([0], [0], color="#4C72B0", lw=1.5, linestyle="--", label="H-bond"),
        plt.Line2D([0], [0], color="#C44E52", lw=2.0, linestyle="-",  label="Salt bridge"),
    ]
    ax.legend(handles=legend_elements, loc="lower center",
              bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=8, frameon=False)
    ax.set_title("2D Interaction Map — PPI Interface",
                 fontsize=12, fontweight="bold", pad=10)

    plt.tight_layout()
    for ext in ("svg", "png"):
        plt.savefig(os.path.join(fig_dir, f"figure_2D_interaction_map.{ext}"),
                    bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close()
    print("  Saved: figure_2D_interaction_map.svg / .png")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="PPI hotspot identification from PDBePISA JSON output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--json", required=True,
        help="Path to PDBePISA interface JSON file."
    )
    parser.add_argument(
        "--outdir", default="results/",
        help="Output directory for CSV tables and figures (default: results/)."
    )
    parser.add_argument(
        "--chain1_label", default="Chain 1",
        help='Display label for chain 1 (default: "Chain 1"). '
             'Example: "Chain F (TNF-alpha)"'
    )
    parser.add_argument(
        "--chain2_label", default="Chain 2",
        help='Display label for chain 2 (default: "Chain 2"). '
             'Example: "Chain L (TNFR1)"'
    )
    parser.add_argument(
        "--top_n", type=int, default=10,
        help="Number of top hotspot residues to display and plot (default: 10)."
    )
    args = parser.parse_args()

    # Set up output directories
    os.makedirs(args.outdir, exist_ok=True)
    fig_dir = os.path.join(args.outdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    sep = "=" * 60
    print(f"\n{sep}")
    print("  PPI Hotspot Identification — PDBePISA JSON Analysis")
    print(f"{sep}\n")

    # ------------------------------------------------------------------
    # Step 1: Load JSON
    # ------------------------------------------------------------------
    print(f"[1] Loading PDBePISA JSON: {args.json}")
    data = load_pisa_json(args.json)

    # Verify expected keys
    top_keys = list(data.keys())
    res_keys  = list(data.get("interfacing_residues", {}).keys())
    bond_keys = list(data.get("interfacing_bonds", {}).keys())
    print(f"    Top-level keys:        {top_keys}")
    print(f"    Interfacing res keys:  {res_keys}")
    print(f"    Interfacing bond keys: {bond_keys}")

    iface = data["interface"]
    print(f"\n    Interface summary:")
    print(f"      Chain 1 ID : {iface.get('chain1_id', 'N/A')}")
    print(f"      Chain 2 ID : {iface.get('chain2_id', 'N/A')}")
    print(f"      Interface ASA : {iface.get('int_area', 'N/A')} Å²")
    print(f"      ΔiG           : {iface.get('int_solv_en', 'N/A')} kcal/mol")
    print(f"      p-value       : {iface.get('pvalue', 'N/A')}")

    # ------------------------------------------------------------------
    # Step 2: Extract interface residues
    # ------------------------------------------------------------------
    print(f"\n[2] Extracting interface residues (BSA > 0 Å²)...")
    df_1 = extract_interface_residues(data, "chain1")
    df_2 = extract_interface_residues(data, "chain2")
    print(f"    {args.chain1_label}: {len(df_1)} residues  "
          f"(total BSA = {df_1['bsa'].sum():.1f} Å²)")
    print(f"    {args.chain2_label}: {len(df_2)} residues  "
          f"(total BSA = {df_2['bsa'].sum():.1f} Å²)")

    # ------------------------------------------------------------------
    # Step 3: Parse bond inventory
    # ------------------------------------------------------------------
    print(f"\n[3] Parsing bond inventory...")
    df_hb   = parse_bonds(data, "hbonds")
    df_sb   = parse_bonds(data, "salt_bridges")
    df_hb_u = get_unique_bonds(df_hb)
    df_sb_u = get_unique_bonds(df_sb)
    print(f"    Hydrogen bonds : {len(df_hb)} entries → {len(df_hb_u)} unique pairs")
    print(f"    Salt bridges   : {len(df_sb)} entries → {len(df_sb_u)} unique pairs")

    if len(df_hb_u) > 0:
        print(f"\n    Hydrogen bonds (unique, sorted by distance):")
        for _, row in df_hb_u.iterrows():
            print(f"      {row['chain1_res']}{row['chain1_seqid']} {row['chain1_atom']}"
                  f"  ↔  {row['chain2_res']}{row['chain2_seqid']} {row['chain2_atom']}"
                  f"  {row['distance']:.3f} Å")

    if len(df_sb_u) > 0:
        print(f"\n    Salt bridges (unique, sorted by distance):")
        for _, row in df_sb_u.iterrows():
            print(f"      {row['chain1_res']}{row['chain1_seqid']} {row['chain1_atom']}"
                  f"  ↔  {row['chain2_res']}{row['chain2_seqid']} {row['chain2_atom']}"
                  f"  {row['distance']:.3f} Å")

    # ------------------------------------------------------------------
    # Step 4: Compute hotspot scores
    # ------------------------------------------------------------------
    print(f"\n[4] Computing composite hotspot scores (BSA_NORM = data-driven)...")
    df_1_scored = compute_hotspot_scores(df_1)
    df_2_scored = compute_hotspot_scores(df_2)

    print(f"\n    Top {args.top_n} hotspot residues — {args.chain1_label}:")
    cols = ["rank", "residue", "seq_id", "bsa", "solv_en", "interaction", "hotspot_score"]
    for _, row in df_1_scored.head(args.top_n).iterrows():
        print(f"      #{int(row['rank']):>2}: {row['residue']}{row['seq_id']:<6}  "
              f"BSA={row['bsa']:>6.1f} Å²  ΔG={row['solv_en']:>+7.3f}  "
              f"Score={row['hotspot_score']:.2f}  [{row['interaction']}]")

    # ------------------------------------------------------------------
    # Step 5: Save CSV outputs
    # ------------------------------------------------------------------
    print(f"\n[5] Saving CSV outputs to: {args.outdir}")
    df_1.to_csv(os.path.join(args.outdir, "interface_residues_chain1.csv"), index=False)
    df_2.to_csv(os.path.join(args.outdir, "interface_residues_chain2.csv"), index=False)
    df_1_scored.to_csv(os.path.join(args.outdir, "hotspot_scores_chain1.csv"), index=False)
    df_2_scored.to_csv(os.path.join(args.outdir, "hotspot_scores_chain2.csv"), index=False)
    df_hb.to_csv(os.path.join(args.outdir, "bond_inventory_hbonds.csv"), index=False)
    df_sb.to_csv(os.path.join(args.outdir, "bond_inventory_saltbridges.csv"), index=False)
    print("    Saved 6 CSV files.")

    # ------------------------------------------------------------------
    # Step 6: Generate figures
    # ------------------------------------------------------------------
    print(f"\n[6] Generating figures to: {fig_dir}")
    plot_hotspot_ranking(
        df_1_scored, df_2_scored, fig_dir,
        args.chain1_label, args.chain2_label, top_n=args.top_n
    )
    plot_bsa_vs_energy(
        df_1_scored, df_2_scored, fig_dir,
        args.chain1_label, args.chain2_label
    )
    plot_2d_interaction_map(
        df_hb_u, df_sb_u, fig_dir,
        args.chain1_label, args.chain2_label
    )

    print(f"\n{sep}")
    print("  Analysis complete.")
    print(f"  Tables : {args.outdir}")
    print(f"  Figures: {fig_dir}")
    print(f"{sep}\n")


if __name__ == "__main__":
    main()
