# PPI Hotspot Identification via PDBePISA
### A Computational Pipeline for Protein–Protein Interface Analysis and Binding Determinant Characterization

> **Case study:** TNF-α / TNFR1 complex (PDB: [8ZUI](https://www.rcsb.org/structure/8ZUI))

---

## Table of Contents

1. [Overview](#overview)
2. [Scientific Background](#scientific-background)
3. [Repository Structure](#repository-structure)
4. [Dependencies](#dependencies)
5. [Workflow Overview](#workflow-overview)
6. [Step-by-Step Pipeline](#step-by-step-pipeline)
   - [Step 1 — Retrieve the Structure from PDB](#step-1--retrieve-the-structure-from-pdb)
   - [Step 2 — Run PDBePISA Interface Analysis](#step-2--run-pdbepisa-interface-analysis)
   - [Step 3 — Export PDBePISA Output as JSON](#step-3--export-pdbepisa-output-as-json)
   - [Step 4 — Parse the JSON Output with Python](#step-4--parse-the-json-output-with-python)
   - [Step 5 — Extract Interface Residues](#step-5--extract-interface-residues)
   - [Step 6 — Parse the Bond Inventory](#step-6--parse-the-bond-inventory)
   - [Step 7 — Compute Composite Hotspot Scores](#step-7--compute-composite-hotspot-scores)
   - [Step 8 — Visualize Results](#step-8--visualize-results)
7. [Complete Python Script](#complete-python-script)
8. [Results: TNF-α / TNFR1 Case Study](#results-tnf--tnfr1-case-study)
   - [Interface Metrics](#interface-metrics)
   - [Hydrogen Bond Network](#hydrogen-bond-network)
   - [Salt Bridge Cluster](#salt-bridge-cluster)
   - [Hotspot Residue Rankings](#hotspot-residue-rankings)
   - [Biological Interpretation](#biological-interpretation)
9. [Output Files](#output-files)
10. [PyMOL Visualization Scripts](#pymol-visualization-scripts)
11. [References](#references)

---

## Overview

This repository provides a complete, reproducible pipeline for identifying **hotspot residues** at protein–protein interfaces (PPIs) using thermodynamic data from [PDBePISA](https://www.ebi.ac.uk/pdbe/pisa/) (Proteins, Interfaces, Structures and Assemblies; EMBL-EBI). The pipeline is demonstrated on the **TNF-α / TNFR1** complex (PDB: 8ZUI), a clinically important cytokine–receptor interaction and validated drug target in inflammatory disease.

The pipeline:
- Retrieves and validates a PPI structure from the RCSB Protein Data Bank
- Runs PDBePISA to compute buried surface areas (BSA), solvation energies (ΔiG), and non-covalent bond inventories
- Parses the PDBePISA JSON output programmatically using Python
- Ranks interface residues by a composite hotspot score integrating BSA, solvation energy, and interaction type
- Generates publication-quality figures

This workflow is directly applicable to any protein–protein complex deposited in the PDB and is designed to support **structure-based drug discovery** targeting PPIs.

---

## Scientific Background

### What is a Hotspot Residue?

Protein–protein interfaces are not energetically uniform. A small subset of interface residues — typically 10–20% — contribute disproportionately to binding free energy. These are termed **hotspot residues** and are operationally defined as residues whose alanine substitution (Ala-scan mutagenesis) reduces binding affinity by ≥2 kcal/mol. Hotspot residues are the primary pharmacological targets for PPI inhibitor design, as small molecules or peptides that mimic or disrupt these contacts can block the interaction with high selectivity.

### What is PDBePISA?

**PDBePISA** (Proteins, Interfaces, Structures and Assemblies) is a web server and API developed by the European Bioinformatics Institute (EMBL-EBI) that performs thermodynamic analysis of macromolecular interfaces from crystal structures. For each interface, PDBePISA computes:

| Parameter | Description |
|---|---|
| **BSA** | Buried Surface Area (Å²) — solvent-accessible surface area lost upon complexation |
| **ΔiG** | Solvation energy gain upon complexation (kcal/mol) — negative values indicate favorable burial |
| **H-bonds** | Hydrogen bonds with donor/acceptor atoms and distances (Å) |
| **Salt bridges** | Electrostatic interactions between oppositely charged residues |
| **p-value** | Statistical significance of the solvation energy gain |

PDBePISA output can be exported as JSON, which enables systematic, reproducible downstream analysis in Python.

### Why TNF-α / TNFR1?

Tumor Necrosis Factor-alpha (TNF-α) is a homotrimeric cytokine that signals through TNFR1 (TNF Receptor 1) to activate NF-κB and apoptotic pathways. The TNF–TNFR1 interaction is a validated therapeutic target in rheumatoid arthritis, Crohn's disease, and other inflammatory conditions. Blocking this PPI with biologics (e.g., etanercept, adalimumab) is clinically effective, but small-molecule PPI inhibitors remain an active area of research. This interface therefore serves as an ideal case study for hotspot identification and structure-based drug design.

---

## Repository Structure

```
ppi-hotspot-identification/
├── data/
│   └── interface_7.json          # PDBePISA JSON output for PDB 8ZUI (Chain F vs. Chain L)
├── scripts/
│   └── ppi_hotspot_analysis.py   # Main analysis script (see below)
├── figures/
│   ├── figure_hotspot_ranking.svg
│   ├── figure_2D_interaction_map.svg
│   └── figure_BSA_vs_energy.svg
├── results/
│   ├── interface_residues_chainF.csv
│   ├── interface_residues_chainL.csv
│   ├── hotspot_scores.csv
│   └── bond_inventory.csv
└── README.md
```

---

## Dependencies

All dependencies are standard scientific Python packages:

```bash
pip install pandas numpy matplotlib
```

| Package | Version | Purpose |
|---|---|---|
| `pandas` | ≥1.3 | DataFrame construction, filtering, and tabular output |
| `numpy` | ≥1.21 | Numerical operations and normalization |
| `matplotlib` | ≥3.4 | Figure generation (bar charts, scatter plots, interaction maps) |
| `json` | stdlib | Parsing PDBePISA JSON output |

Python ≥ 3.8 is required.

---

## Workflow Overview

```
PDB Structure (8ZUI)
        │
        ▼
  PDBePISA Web Server
  (Interface Analysis)
        │
        ▼
  JSON Export (interface_7.json)
        │
        ▼
  Python Parsing (pandas + numpy)
  ┌─────────────────────────────┐
  │  1. Parse top-level keys    │
  │  2. Extract interface data  │
  │  3. Filter BSA > 0          │
  │  4. Assign interaction types│
  │  5. Parse bond inventory    │
  │  6. Compute hotspot scores  │
  └─────────────────────────────┘
        │
        ▼
  Ranked Hotspot Residues
  + Bond Tables
  + Publication Figures
```

---

## Step-by-Step Pipeline

### Step 1 — Retrieve the Structure from PDB

Navigate to the [RCSB Protein Data Bank](https://www.rcsb.org/) and search for your structure of interest. For this case study, the TNF-α / TNFR1 complex is deposited under PDB ID **8ZUI**.

Download the structure in PDB or mmCIF format:

```bash
# Using wget
wget https://files.rcsb.org/download/8ZUI.pdb -O 8ZUI.pdb

# Or using the RCSB REST API
curl -o 8ZUI.pdb "https://files.rcsb.org/download/8ZUI.pdb"
```

Inspect the SEQRES and ATOM records to identify the chain identifiers for your protein pair of interest. For 8ZUI:
- **Chain F** = TNF-α (ligand)
- **Chain L** = TNFR1 (receptor)

![image](https://github.com/abdulaziz-khaled/Hotspot-Mapper/blob/main/screen%20chain%20F%20and%20L.png)

---

### Step 2 — Run PDBePISA Interface Analysis

1. Navigate to [https://www.ebi.ac.uk/pdbe/pisa/](https://www.ebi.ac.uk/pdbe/pisa/)
2. Enter your PDB ID (e.g., `8ZUI`) in the search box and click **Analyse**
3. In the results page, click on the **Interfaces** tab
4. Identify the interface corresponding to your chain pair of interest (Chain F vs. Chain L)
5. Note the interface index number (e.g., Interface 7 for the F–L pair in 8ZUI)

PDBePISA will display:
- Interface area (BSA per chain)
- Solvation energy (ΔiG per chain)
- p-value
- Interfacing residues table
- Non-covalent bond inventory (H-bonds, salt bridges)

**Validation criteria for a biologically relevant interface:**
- BSA per chain > 500 Å² (typically 800–2,000 Å² for genuine PPIs)
- ΔiG < 0 kcal/mol (favorable solvation energy gain)
- Symmetric residue counts and near-equal BSA on both sides
- Presence of specific polar contacts (H-bonds, salt bridges)

---

### Step 3 — Export PDBePISA Output as JSON

PDBePISA provides a programmatic API that returns interface data in JSON format. This is the recommended approach for reproducible, automated analysis.

**API endpoint format:**
```
https://www.ebi.ac.uk/pdbe/api/pisa/interfacelist/{pdb_id}
https://www.ebi.ac.uk/pdbe/api/pisa/interfacedetail/{pdb_id}/{interface_id}
```

Alternatively, the JSON can be exported directly from the PDBePISA web interface by selecting the interface of interest and using the download/export option.

The exported JSON file (`interface_7.json`) contains three top-level keys:

```json
{
  "interface":            { ... },   // Global interface metrics (BSA, ΔiG, p-value, atom counts)
  "interfacing_residues": { ... },   // Per-residue data for both chains
  "interfacing_bonds":    { ... }    // Explicit bond inventory (H-bonds, salt bridges)
}
```

---

### Step 4 — Parse the JSON Output with Python

Load the JSON file and inspect its structure:

```python
import json
import pandas as pd
import numpy as np

# Load PDBePISA JSON output
with open("data/interface_7.json", "r") as f:
    data = json.load(f)

# Inspect top-level keys
print("Top-level keys:", list(data.keys()))
# Output: ['interface', 'interfacing_residues', 'interfacing_bonds']

# Inspect interface-level metrics
iface = data["interface"]
print(f"Chain 1: {iface['chain1_id']}, Chain 2: {iface['chain2_id']}")
print(f"Interface ASA: {iface['int_area']} Å²")
print(f"ΔiG: {iface['int_solv_en']} kcal/mol")
print(f"p-value: {iface['pvalue']}")
```

**Understanding the JSON structure:**

The `interfacing_residues` key contains a nested dictionary with entries for each chain. Each residue entry includes:

| JSON field | Description |
|---|---|
| `auth_seq_id` | Author sequence number (matches PDB ATOM records) |
| `residue_name` | Three-letter amino acid code |
| `bsa` | Buried surface area (Å²) upon complexation |
| `solv_en` | Per-residue solvation energy contribution (kcal/mol) |
| `bond_type` | Interaction type: `"H"` (H-bond), `"S"` (salt bridge), `"HS"` (both), or `null` |

The `interfacing_bonds` key contains separate sub-dictionaries for `hbonds` and `salt_bridges`, each listing atom-level contacts with distances.

---

### Step 5 — Extract Interface Residues

Filter residues with BSA > 0 Å² (the PDBePISA definition of an interface residue) and assign interaction categories:

```python
def extract_interface_residues(data, chain_key):
    """
    Extract interface residues for a given chain from PDBePISA JSON output.
    
    Parameters
    ----------
    data : dict
        Parsed PDBePISA JSON (top-level dictionary).
    chain_key : str
        Key identifying the chain in 'interfacing_residues' (e.g., 'chain1', 'chain2').
    
    Returns
    -------
    pd.DataFrame
        DataFrame of interface residues with BSA > 0, sorted by BSA descending.
    """
    residues = data["interfacing_residues"][chain_key]
    records = []
    for res in residues:
        bsa = float(res.get("bsa", 0))
        if bsa > 0:
            bond_type = res.get("bond_type", None)
            # Map PDBePISA bond_type codes to human-readable labels
            if bond_type == "H":
                interaction = "H-bond"
            elif bond_type == "S":
                interaction = "Salt bridge"
            elif bond_type == "HS":
                interaction = "H-bond + Salt bridge"
            else:
                interaction = "Hydrophobic/VdW"
            records.append({
                "residue":      res["residue_name"],
                "seq_id":       int(res["auth_seq_id"]),
                "bsa":          bsa,
                "solv_en":      float(res.get("solv_en", 0)),
                "bond_type":    bond_type,
                "interaction":  interaction,
            })
    df = pd.DataFrame(records).sort_values("bsa", ascending=False).reset_index(drop=True)
    return df

# Extract for both chains
df_F = extract_interface_residues(data, "chain1")   # Chain F = TNF-α
df_L = extract_interface_residues(data, "chain2")   # Chain L = TNFR1

print(f"Chain F interface residues: {len(df_F)}")
print(f"Chain L interface residues: {len(df_L)}")
print(f"\nChain F total BSA: {df_F['bsa'].sum():.1f} Å²")
print(f"Chain L total BSA: {df_L['bsa'].sum():.1f} Å²")
```

---

### Step 6 — Parse the Bond Inventory

Extract explicit hydrogen bonds and salt bridges with atom-level detail and distances:

```python
def parse_bonds(data, bond_category):
    """
    Parse the explicit bond inventory from PDBePISA JSON output.
    
    Parameters
    ----------
    data : dict
        Parsed PDBePISA JSON.
    bond_category : str
        Either 'hbonds' or 'salt_bridges'.
    
    Returns
    -------
    pd.DataFrame
        DataFrame of bonds with donor/acceptor residues, atoms, and distances.
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

df_hbonds  = parse_bonds(data, "hbonds")
df_sbridges = parse_bonds(data, "salt_bridges")

print("Hydrogen bonds:")
print(df_hbonds.to_string(index=False))
print("\nSalt bridges:")
print(df_sbridges.to_string(index=False))
```

> **Note on PDBePISA bond directionality:** PDBePISA reports each cross-chain bond in **both directions** (Chain F → Chain L and Chain L → Chain F). When reporting unique bond pairs, deduplicate by keeping only the canonical direction (lower chain ID as donor) or by filtering on `chain1_seqid < chain2_seqid`.

---

### Step 7 — Compute Composite Hotspot Scores

Rank interface residues by a composite score that integrates three complementary descriptors:

```python
def compute_hotspot_scores(df, bsa_norm=140.0):
    """
    Compute composite hotspot scores for interface residues.
    
    The composite score integrates:
      (1) Normalized BSA: captures the geometric contribution (surface burial)
      (2) |ΔG_solv|: captures the thermodynamic contribution (solvation energy)
      (3) Bond-type bonus: rewards residues engaged in specific polar contacts
    
    Score = (BSA / bsa_norm) + |ΔG_solv| + bond_bonus
    
    Bond bonuses:
      - Hydrogen bond only:          +1.0
      - Salt bridge only:            +1.5
      - Hydrogen bond + salt bridge: +2.5
      - Hydrophobic/VdW only:        +0.0
    
    Parameters
    ----------
    df : pd.DataFrame
        Interface residue DataFrame from extract_interface_residues().
    bsa_norm : float
        Normalization factor for BSA (default: maximum BSA in dataset, ~140 Å²).
    
    Returns
    -------
    pd.DataFrame
        Input DataFrame with added 'hotspot_score' column, sorted descending.
    """
    bond_bonus_map = {
        "H":    1.0,   # Hydrogen bond
        "S":    1.5,   # Salt bridge
        "HS":   2.5,   # Hydrogen bond + salt bridge
        None:   0.0,   # Hydrophobic / van der Waals
    }
    df = df.copy()
    df["bond_bonus"]     = df["bond_type"].map(bond_bonus_map).fillna(0.0)
    df["hotspot_score"]  = (
        df["bsa"] / bsa_norm
        + df["solv_en"].abs()
        + df["bond_bonus"]
    )
    df = df.sort_values("hotspot_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df

# Compute hotspot scores for both chains
df_F_scored = compute_hotspot_scores(df_F)
df_L_scored = compute_hotspot_scores(df_L)

# Display top 10 hotspot residues
cols = ["rank", "residue", "seq_id", "bsa", "solv_en", "interaction", "hotspot_score"]
print("Top 10 hotspot residues — Chain F (TNF-α):")
print(df_F_scored[cols].head(10).to_string(index=False))
```

**Rationale for the composite score:**

| Component | Rationale |
|---|---|
| `BSA / 140` | Residues that bury more surface area upon binding contribute more to shape complementarity and are harder to replace. Normalized to the maximum observed BSA in the dataset. |
| `\|ΔG_solv\|` | The absolute solvation energy contribution captures both favorable (negative ΔG, hydrophobic burial) and unfavorable (positive ΔG, polar desolvation) contributions. Both extremes indicate energetically important residues. |
| Bond bonus | Residues forming specific polar contacts (H-bonds, salt bridges) are known to be disproportionately important for binding specificity and affinity. The bonus weights these residues above purely hydrophobic contacts. |

---

### Step 8 — Visualize Results

Generate three publication-quality figures:

**Figure 1 — Hotspot Ranking Bar Chart:**
```python
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'

color_map = {
    "H-bond":                "#4C72B0",
    "Salt bridge":           "#DD8452",
    "H-bond + Salt bridge":  "#55A868",
    "Hydrophobic/VdW":       "#C0C0C0",
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, df_scored, chain_label in zip(
    axes,
    [df_F_scored, df_L_scored],
    ["Chain F (TNF-α)", "Chain L (TNFR1)"]
):
    top = df_scored.head(10)
    labels = [f"{r['residue']}{r['seq_id']}" for _, r in top.iterrows()]
    colors = [color_map[i] for i in top["interaction"]]
    ax.barh(labels[::-1], top["hotspot_score"].values[::-1], color=colors[::-1])
    ax.set_xlabel("Composite Hotspot Score", fontsize=11)
    ax.set_title(chain_label, fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig("figures/figure_hotspot_ranking.svg", bbox_inches="tight")
plt.savefig("figures/figure_hotspot_ranking.png", dpi=300, bbox_inches="tight")
```

**Figure 2 — BSA vs. Solvation Energy Scatter:**
```python
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, df_scored, chain_label in zip(
    axes,
    [df_F_scored, df_L_scored],
    ["Chain F (TNF-α)", "Chain L (TNFR1)"]
):
    for interaction, group in df_scored.groupby("interaction"):
        ax.scatter(group["bsa"], group["solv_en"],
                   label=interaction, color=color_map[interaction],
                   s=60, alpha=0.85, edgecolors="white", linewidths=0.5)
    # Label top 5 hotspots
    for _, row in df_scored.head(5).iterrows():
        ax.annotate(f"{row['residue']}{row['seq_id']}",
                    (row["bsa"], row["solv_en"]),
                    fontsize=8, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Buried Surface Area (Å²)", fontsize=11)
    ax.set_ylabel("Solvation Energy (kcal/mol)", fontsize=11)
    ax.set_title(chain_label, fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)

axes[0].legend(fontsize=8, frameon=False)
plt.tight_layout()
plt.savefig("figures/figure_BSA_vs_energy.svg", bbox_inches="tight")
```

---

## Complete Python Script

The full self-contained analysis script is provided below. Copy it to `scripts/ppi_hotspot_analysis.py` and run with:

```bash
python scripts/ppi_hotspot_analysis.py --json data/interface_7.json --outdir results/
```

```python
#!/usr/bin/env python3
"""
ppi_hotspot_analysis.py
=======================
Hotspot identification at protein–protein interfaces from PDBePISA JSON output.

Usage
-----
    python ppi_hotspot_analysis.py --json <path_to_pisa_json> --outdir <output_directory>

Input
-----
    PDBePISA JSON export for a single interface (keys: interface,
    interfacing_residues, interfacing_bonds).

Output
------
    - interface_residues_chain1.csv / interface_residues_chain2.csv
    - hotspot_scores_chain1.csv / hotspot_scores_chain2.csv
    - bond_inventory_hbonds.csv / bond_inventory_saltbridges.csv
    - figure_hotspot_ranking.svg / .png
    - figure_BSA_vs_energy.svg / .png
    - figure_2D_interaction_map.svg / .png

References
----------
    Krissinel & Henrick (2007) J. Mol. Biol. 372:774-797
    Krissinel (2015) Nucleic Acids Res. 43:W314-W319
    Berman et al. (2000) Nucleic Acids Res. 28:235-242
    McKinney (2010) Proc. SciPy 2010, 51-56
    Harris et al. (2020) Nature 585:357-362
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

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

BOND_BONUS = {
    "H":  1.0,   # Hydrogen bond only
    "S":  1.5,   # Salt bridge only
    "HS": 2.5,   # Hydrogen bond + salt bridge
    None: 0.0,   # Hydrophobic / van der Waals
}

INTERACTION_LABEL = {
    "H":  "H-bond",
    "S":  "Salt bridge",
    "HS": "H-bond + Salt bridge",
    None: "Hydrophobic/VdW",
}

COLOR_MAP = {
    "H-bond":                "#4C72B0",
    "Salt bridge":           "#DD8452",
    "H-bond + Salt bridge":  "#55A868",
    "Hydrophobic/VdW":       "#C0C0C0",
}

BSA_NORM = 140.0   # Normalization factor (Å²) — maximum BSA in this dataset


# ─────────────────────────────────────────────────────────────────────────────
# Parsing functions
# ─────────────────────────────────────────────────────────────────────────────

def load_pisa_json(filepath):
    """Load and return parsed PDBePISA JSON."""
    with open(filepath, "r") as f:
        return json.load(f)


def extract_interface_residues(data, chain_key):
    """
    Extract interface residues (BSA > 0) for a given chain.

    Parameters
    ----------
    data : dict
        Parsed PDBePISA JSON.
    chain_key : str
        Key in 'interfacing_residues' dict (e.g., 'chain1', 'chain2').

    Returns
    -------
    pd.DataFrame
        Interface residues sorted by BSA descending.
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
    Parse explicit bond inventory (hbonds or salt_bridges).

    Parameters
    ----------
    data : dict
        Parsed PDBePISA JSON.
    bond_category : str
        'hbonds' or 'salt_bridges'.

    Returns
    -------
    pd.DataFrame
        Bond table sorted by distance ascending.
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
    Deduplicate PDBePISA bond table (which lists each bond bidirectionally).
    Keeps the entry with the smaller chain1_seqid as the canonical direction.
    """
    seen = set()
    unique = []
    for _, row in df_bonds.iterrows():
        key = tuple(sorted([(row["chain1_seqid"], row["chain1_atom"]),
                             (row["chain2_seqid"], row["chain2_atom"])]))
        if key not in seen:
            seen.add(key)
            unique.append(row)
    return pd.DataFrame(unique).reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# Hotspot scoring
# ─────────────────────────────────────────────────────────────────────────────

def compute_hotspot_scores(df, bsa_norm=BSA_NORM):
    """
    Compute composite hotspot scores.

    Score = (BSA / bsa_norm) + |ΔG_solv| + bond_bonus

    Parameters
    ----------
    df : pd.DataFrame
        Interface residue DataFrame.
    bsa_norm : float
        BSA normalization factor (Å²).

    Returns
    -------
    pd.DataFrame
        Input DataFrame with 'bond_bonus', 'hotspot_score', and 'rank' columns.
    """
    df = df.copy()
    df["bond_bonus"]    = df["bond_type"].map(BOND_BONUS).fillna(0.0)
    df["hotspot_score"] = (
        df["bsa"] / bsa_norm
        + df["solv_en"].abs()
        + df["bond_bonus"]
    )
    df = df.sort_values("hotspot_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

def plot_hotspot_ranking(df_F, df_L, outdir, top_n=10):
    """Horizontal bar chart of top hotspot residues for both chains."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, df_scored, chain_label in zip(
        axes,
        [df_F, df_L],
        ["Chain F (TNF-α)", "Chain L (TNFR1)"]
    ):
        top = df_scored.head(top_n)
        labels = [f"{r['residue']}{r['seq_id']}" for _, r in top.iterrows()]
        colors = [COLOR_MAP[i] for i in top["interaction"]]
        ax.barh(labels[::-1], top["hotspot_score"].values[::-1],
                color=colors[::-1], edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Composite Hotspot Score", fontsize=11)
        ax.set_title(chain_label, fontsize=12, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=9)

    # Shared legend
    patches = [mpatches.Patch(color=v, label=k) for k, v in COLOR_MAP.items()]
    fig.legend(handles=patches, loc="lower center", ncol=4,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.04))
    plt.suptitle("Interface Residue Hotspot Ranking", fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    for ext in ("svg", "png"):
        path = os.path.join(outdir, f"figure_hotspot_ranking.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close()
    print(f"  Saved: figure_hotspot_ranking.svg / .png")


def plot_bsa_vs_energy(df_F, df_L, outdir, top_label=5):
    """Scatter plot of BSA vs. solvation energy for both chains."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, df_scored, chain_label in zip(
        axes,
        [df_F, df_L],
        ["Chain F (TNF-α)", "Chain L (TNFR1)"]
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
        ax.set_title(chain_label, fontsize=12, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.tick_params(axis="both", labelsize=9)

    axes[0].legend(fontsize=8, frameon=False)
    plt.suptitle("BSA vs. Solvation Energy at the PPI Interface",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    for ext in ("svg", "png"):
        path = os.path.join(outdir, f"figure_BSA_vs_energy.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close()
    print(f"  Saved: figure_BSA_vs_energy.svg / .png")


def plot_2d_interaction_map(df_hbonds_unique, df_sb_unique, outdir):
    """
    2D schematic interaction map of cross-chain bonds.
    Chain F residues on the left, Chain L residues on the right.
    H-bonds: dashed blue lines. Salt bridges: solid red lines.
    """
    # Collect unique residues from both bond tables
    f_residues, l_residues = [], []
    for _, row in pd.concat([df_hbonds_unique, df_sb_unique]).iterrows():
        label_f = f"{row['chain1_res']}{row['chain1_seqid']}"
        label_l = f"{row['chain2_res']}{row['chain2_seqid']}"
        if label_f not in f_residues:
            f_residues.append(label_f)
        if label_l not in l_residues:
            l_residues.append(label_l)

    n_f = len(f_residues)
    n_l = len(l_residues)
    height = max(n_f, n_l) * 1.2 + 2

    fig, ax = plt.subplots(figsize=(9, height))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, height)
    ax.axis("off")

    # Residue node colors by physicochemical property
    def node_color(res_label):
        aa = res_label[:3]
        positive = {"LYS", "ARG", "HIS"}
        negative = {"ASP", "GLU"}
        polar    = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
        glycine  = {"GLY"}
        if aa in positive: return "#4C72B0"
        if aa in negative: return "#C44E52"
        if aa in polar:    return "#55A868"
        if aa in glycine:  return "#AAAAAA"
        return "#DDDDDD"

    y_f = {res: height - 1.5 - i * (height - 2) / max(n_f - 1, 1)
           for i, res in enumerate(f_residues)}
    y_l = {res: height - 1.5 - i * (height - 2) / max(n_l - 1, 1)
           for i, res in enumerate(l_residues)}

    # Draw nodes
    for res, y in y_f.items():
        ax.add_patch(plt.Circle((2, y), 0.35, color=node_color(res), zorder=3))
        ax.text(1.5, y, res, ha="right", va="center", fontsize=9, fontweight="bold")
    for res, y in y_l.items():
        ax.add_patch(plt.Circle((8, y), 0.35, color=node_color(res), zorder=3))
        ax.text(8.5, y, res, ha="left", va="center", fontsize=9, fontweight="bold")

    # Draw bonds
    for _, row in df_hbonds_unique.iterrows():
        lf = f"{row['chain1_res']}{row['chain1_seqid']}"
        ll = f"{row['chain2_res']}{row['chain2_seqid']}"
        if lf in y_f and ll in y_l:
            ax.plot([2.35, 7.65], [y_f[lf], y_l[ll]],
                    color="#4C72B0", linewidth=1.5, linestyle="--",
                    alpha=0.8, zorder=2)
            mx, my = 5.0, (y_f[lf] + y_l[ll]) / 2
            ax.text(mx, my + 0.15, f"{row['distance']:.3f} Å",
                    ha="center", va="bottom", fontsize=7.5, color="#4C72B0")

    for _, row in df_sb_unique.iterrows():
        lf = f"{row['chain1_res']}{row['chain1_seqid']}"
        ll = f"{row['chain2_res']}{row['chain2_seqid']}"
        if lf in y_f and ll in y_l:
            ax.plot([2.35, 7.65], [y_f[lf], y_l[ll]],
                    color="#C44E52", linewidth=2.0, linestyle="-",
                    alpha=0.85, zorder=2)
            mx, my = 5.0, (y_f[lf] + y_l[ll]) / 2
            ax.text(mx, my - 0.15, f"{row['distance']:.3f} Å",
                    ha="center", va="top", fontsize=7.5, color="#C44E52")

    # Column headers
    ax.text(2, height - 0.5, "Chain F\n(TNF-α)",
            ha="center", va="center", fontsize=10, fontweight="bold")
    ax.text(8, height - 0.5, "Chain L\n(TNFR1)",
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

    ax.set_title("2D Interaction Map — TNF-α / TNFR1 Interface",
                 fontsize=12, fontweight="bold", pad=10)
    plt.tight_layout()
    for ext in ("svg", "png"):
        path = os.path.join(outdir, f"figure_2D_interaction_map.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close()
    print(f"  Saved: figure_2D_interaction_map.svg / .png")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="PPI hotspot identification from PDBePISA JSON output."
    )
    parser.add_argument("--json",   required=True, help="Path to PDBePISA JSON file")
    parser.add_argument("--outdir", default="results/", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fig_dir = os.path.join(args.outdir, "..", "figures")
    os.makedirs(fig_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print("  PPI Hotspot Identification — PDBePISA JSON Analysis")
    print(f"{'='*60}\n")

    # 1. Load JSON
    print(f"[1] Loading PDBePISA JSON: {args.json}")
    data = load_pisa_json(args.json)
    iface = data["interface"]
    print(f"    Chain 1: {iface.get('chain1_id', 'N/A')}  |  "
          f"Chain 2: {iface.get('chain2_id', 'N/A')}")
    print(f"    Interface ASA: {iface.get('int_area', 'N/A')} Å²")
    print(f"    ΔiG: {iface.get('int_solv_en', 'N/A')} kcal/mol")
    print(f"    p-value: {iface.get('pvalue', 'N/A')}")

    # 2. Extract interface residues
    print("\n[2] Extracting interface residues (BSA > 0)...")
    df_F = extract_interface_residues(data, "chain1")
    df_L = extract_interface_residues(data, "chain2")
    print(f"    Chain 1: {len(df_F)} residues  |  Chain 2: {len(df_L)} residues")

    # 3. Parse bonds
    print("\n[3] Parsing bond inventory...")
    df_hb  = parse_bonds(data, "hbonds")
    df_sb  = parse_bonds(data, "salt_bridges")
    df_hb_u = get_unique_bonds(df_hb)
    df_sb_u = get_unique_bonds(df_sb)
    print(f"    Hydrogen bonds: {len(df_hb)} entries ({len(df_hb_u)} unique pairs)")
    print(f"    Salt bridges:   {len(df_sb)} entries ({len(df_sb_u)} unique pairs)")

    # 4. Compute hotspot scores
    print("\n[4] Computing composite hotspot scores...")
    df_F_scored = compute_hotspot_scores(df_F)
    df_L_scored = compute_hotspot_scores(df_L)
    print("    Top 5 hotspot residues (Chain 1):")
    for _, row in df_F_scored.head(5).iterrows():
        print(f"      #{int(row['rank'])}: {row['residue']}{row['seq_id']}  "
              f"BSA={row['bsa']:.1f} Å²  ΔG={row['solv_en']:.3f}  "
              f"Score={row['hotspot_score']:.2f}  [{row['interaction']}]")

    # 5. Save CSV outputs
    print("\n[5] Saving CSV outputs...")
    df_F.to_csv(os.path.join(args.outdir, "interface_residues_chain1.csv"), index=False)
    df_L.to_csv(os.path.join(args.outdir, "interface_residues_chain2.csv"), index=False)
    df_F_scored.to_csv(os.path.join(args.outdir, "hotspot_scores_chain1.csv"), index=False)
    df_L_scored.to_csv(os.path.join(args.outdir, "hotspot_scores_chain2.csv"), index=False)
    df_hb.to_csv(os.path.join(args.outdir, "bond_inventory_hbonds.csv"), index=False)
    df_sb.to_csv(os.path.join(args.outdir, "bond_inventory_saltbridges.csv"), index=False)
    print("    Saved 6 CSV files.")

    # 6. Generate figures
    print("\n[6] Generating figures...")
    plot_hotspot_ranking(df_F_scored, df_L_scored, fig_dir)
    plot_bsa_vs_energy(df_F_scored, df_L_scored, fig_dir)
    plot_2d_interaction_map(df_hb_u, df_sb_u, fig_dir)

    print(f"\n{'='*60}")
    print("  Analysis complete.")
    print(f"  Results: {args.outdir}")
    print(f"  Figures: {fig_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
```

---

## Results: TNF-α / TNFR1 Case Study

### Interface Metrics

PDBePISA analysis of PDB 8ZUI (Chain F vs. Chain L) identified a well-defined, reciprocal protein–protein interface:

| Parameter | Chain F (TNF-α) | Chain L (TNFR1) |
|---|---|---|
| Interface atoms | 102 | 102 |
| Interface residues | 30 | 30 |
| Interface ASA (Å²) | 958.8 | 954.3 |
| Total ASA (Å²) | 9,916.6 | 9,910.8 |
| % ASA buried | 9.7% | 9.6% |
| ΔiG (kcal/mol) | −6.3 | −4.5 |
| p-value | 0.616 | 0.615 |
| **Total interface area** | **~1,913 Å²** | |
| **Total ΔiG** | **−10.8 kcal/mol** | |

The combined interface area of ~1,913 Å² and total ΔiG of −10.8 kcal/mol are consistent with a genuine, biologically relevant protein–protein interface. The symmetric residue counts (30 per chain) and near-identical ASA values confirm a well-defined, reciprocal interaction surface.

---

### Hydrogen Bond Network

Three unique cross-chain hydrogen bonds stabilize the interface:

| Chain F Residue | Atom | Chain L Residue | Atom | Distance (Å) | Type |
|---|---|---|---|---|---|
| **HIS63** | NE2 (sidechain) | **LYS61** | O (backbone) | **2.542** | Sidechain–backbone |
| **GLN46** | NE2 (sidechain) | **ASP120** | OD2 (sidechain) | 2.868 | Sidechain–sidechain |
| **ASP78** | N (backbone) | **GLY76** | O (backbone) | 2.922 | Backbone–backbone |

> PDBePISA reports each bond bidirectionally; three unique cross-chain pairs are listed here in canonical direction (Chain F → Chain L).

---

### Salt Bridge Cluster

Two unique cross-chain salt bridges define the electrostatic core:

| Chain F Residue | Atom | Chain L Residue | Atom | Distance (Å) |
|---|---|---|---|---|
| **LYS64** | NZ | **GLU93** | OE2 | **2.750** |
| **HIS63** | ND1 | **GLU93** | OE1 | 2.867 |

**Key observation:** GLU93 (Chain L) acts as a **bifurcated salt bridge acceptor**, simultaneously engaging both HIS63 and LYS64 from Chain F. This creates a charge-complementary cluster (HIS63–LYS64 dyad on TNF-α ↔ GLU93 on TNFR1) that constitutes the electrostatic core of the interface.

---

### Hotspot Residue Rankings

| Rank | Residue | BSA (Å²) | ΔG (kcal/mol) | Interaction Type | Hotspot Score |
|---|---|---|---|---|---|
| 1 | **HIS63** | 135.9 / 136.2 | +0.888 / +0.889 | H-bond + Salt bridge | **4.36** |
| 2 | **LYS64** | 57.4 / 57.1 | −1.285 / −1.297 | Salt bridge | **3.20** |
| 3 | **GLN46** | 118.9 / 118.5 | −0.820 / −0.821 | H-bond | **2.67** |
| 4 | **GLU93** | 52.0 / 51.6 | −0.118 / −0.117 | Salt bridge | **1.99** |
| 5 | **ASP78** | 66.7 / 65.8 | +0.319 / +0.317 | H-bond | **1.79** |
| 6 | **GLY76** | 48.0 / 48.2 | −0.083 / −0.080 | H-bond | **1.43** |
| 7 | **ASP120** | 17.8 / 17.5 | −0.043 / −0.039 | H-bond | **1.17** |
| 8 | **LYS61** | 7.0 / 6.1 | −0.025 / −0.019 | H-bond | **1.07** |

*Values reported as Chain F / Chain L. Hotspot Score = (BSA/140) + |ΔG| + bond-type bonus.*

The top 8 hotspot residues are **identical between Chain F and Chain L**, reflecting the symmetric nature of the interface — consistent with the TNF homotrimer architecture where each protomer presents an equivalent binding surface.

---

### Biological Interpretation

The **HIS63–LYS64–GLU93 triad** constitutes the energetic and geometric core of the TNF-α/TNFR1 interface and represents the primary pharmacological target for structure-based inhibitor design:

- **HIS63** (top hotspot, score 4.36): Largest BSA of any interface residue (~136 Å²); dual role in both H-bonding (NE2 → LYS61 backbone, 2.542 Å) and salt bridge formation (ND1 → GLU93 OE1, 2.867 Å). A single residue bridging both polar interaction networks.
- **LYS64** (second hotspot, score 3.20): Strongest salt bridge at the interface (NZ → GLU93 OE2, 2.750 Å); most favorable per-residue solvation energy (−1.285/−1.297 kcal/mol).
- **GLN46** (third hotspot, score 2.67): Largest-BSA hydrogen bond residue (~119 Å²); sidechain–sidechain contact with ASP120 (NE2 → OD2, 2.868 Å).
- **GLU93** (fourth hotspot, score 1.99): Bifurcated acceptor simultaneously engaging HIS63 and LYS64; central node of the electrostatic network.

---

## Output Files

After running the pipeline, the following files are generated:

| File | Description |
|---|---|
| `results/interface_residues_chain1.csv` | All interface residues for Chain F (BSA, ΔG, interaction type) |
| `results/interface_residues_chain2.csv` | All interface residues for Chain L |
| `results/hotspot_scores_chain1.csv` | Hotspot-ranked residues for Chain F |
| `results/hotspot_scores_chain2.csv` | Hotspot-ranked residues for Chain L |
| `results/bond_inventory_hbonds.csv` | Full hydrogen bond inventory with atom-level detail |
| `results/bond_inventory_saltbridges.csv` | Full salt bridge inventory |
| `figures/figure_hotspot_ranking.svg` | Horizontal bar chart of top hotspot residues (both chains) |
| `figures/figure_BSA_vs_energy.svg` | Scatter plot of BSA vs. solvation energy |
| `figures/figure_2D_interaction_map.svg` | 2D schematic of cross-chain H-bonds and salt bridges |

---

## PyMOL Visualization Scripts

Use these commands in PyMOL to generate publication-quality structural figures:

**Figure A — Overall complex (cartoon):**
```python
fetch 8ZUI
hide everything
show cartoon
color marine, chain F
color orange, chain L
color gray80, not (chain F or chain L)
set cartoon_transparency, 0.5, not (chain F or chain L)
orient chain F or chain L
ray 1200, 900
png figure_complex_overview.png, dpi=300
```

**Figure B — Interface residues with bonds (sticks):**
```python
select iface_F, chain F and resi 45+46+47+48+60+61+62+63+64+65+66+76+77+78+79+81+83+92+93+119+120+155+156+159+162+165+166+167+174+177
select iface_L, chain L and resi 45+46+47+48+60+61+62+63+64+65+66+76+77+78+79+81+83+92+93+119+120+155+156+159+162+165+166+167+174+177
show sticks, iface_F or iface_L
color yellow, iface_F and name C*
color salmon, iface_L and name C*

# Hydrogen bonds (dashed yellow)
distance hb1, /8ZUI//F/GLN`46/NE2,  /8ZUI//L/ASP`120/OD2
distance hb2, /8ZUI//F/HIS`63/NE2,  /8ZUI//L/LYS`61/O
distance hb3, /8ZUI//F/ASP`78/N,    /8ZUI//L/GLY`76/O
set dash_color, yellow
set dash_width, 2.5

# Salt bridges (dashed magenta)
distance sb1, /8ZUI//F/HIS`63/ND1,  /8ZUI//L/GLU`93/OE1
distance sb2, /8ZUI//F/LYS`64/NZ,   /8ZUI//L/GLU`93/OE2
color magenta, sb1
color magenta, sb2

label iface_F and name CA, "%s%s" % (resn, resi)
ray 1200, 900
png figure_interface_sticks.png, dpi=300
```

---

## References

1. Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., Shindyalov, I.N., and Bourne, P.E. (2000). The Protein Data Bank. *Nucleic Acids Research*, **28**(1), 235–242. https://doi.org/10.1093/nar/28.1.235

2. Krissinel, E., and Henrick, K. (2007). Inference of macromolecular assemblies from crystalline state. *Journal of Molecular Biology*, **372**(3), 774–797. https://doi.org/10.1016/j.jmb.2007.05.022

3. Krissinel, E. (2015). Stock-based detection of protein oligomeric states in jsPISA. *Nucleic Acids Research*, **43**(W1), W314–W319. https://doi.org/10.1093/nar/gkv314

4. McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference (SciPy 2010)*, 51–56. https://doi.org/10.25080/Majora-92bf1922-00a

5. Harris, C.R., Millman, K.J., van der Walt, S.J., et al. (2020). Array programming with NumPy. *Nature*, **585**, 357–362. https://doi.org/10.1038/s41586-020-2649-2

---

*Pipeline developed for the structural analysis of the TNF-α/TNFR1 complex (PDB: 8ZUI). All data derived from PDBePISA thermodynamic analysis and processed with Python/pandas/numpy. Figures generated with matplotlib.*

---

## Adapting This Pipeline to Any PPI Complex

The pipeline is designed to be transferable. If your new PDB complex meets the same structural conditions as 8ZUI, you can swap the JSON file and run the script with minimal changes. This section documents exactly what is fixed, what requires verification, and what must be updated.

---

### What Transfers Directly (No Changes Needed)

If your new complex produces a standard PDBePISA JSON export, the following components work unchanged:

- The JSON parsing logic (`extract_interface_residues`, `parse_bonds`) reads PDBePISA's standardized output format, which is identical regardless of which complex you analyze
- The hotspot scoring formula contains no complex-specific constants (after the fix below)
- All figure generation code is fully data-driven
- The CLI (`--json` / `--outdir`) accepts any input file

---

### What You Must Verify Before Running

#### 1. Identify the correct interface index in PDBePISA

PDBePISA enumerates **all interfaces** in a structure sequentially. For 8ZUI, the TNF-α/TNFR1 interface was Interface 7. For a different complex, the biologically relevant interface may be Interface 1, 3, or 12. **You must identify the correct index** from the PDBePISA web interface before exporting the JSON.

> Exporting the wrong interface produces no error — it silently gives you results for the wrong chain pair.

**How to identify the correct interface:**
1. Go to [https://www.ebi.ac.uk/pdbe/pisa/](https://www.ebi.ac.uk/pdbe/pisa/) and enter your PDB ID
2. Click the **Interfaces** tab
3. Find the row where **Chain 1** and **Chain 2** match your biological pair of interest
4. Note the interface index number and export that specific interface as JSON

#### 2. Confirm the chain key names in the JSON

The script uses:
```python
df_chain1 = extract_interface_residues(data, "chain1")
df_chain2 = extract_interface_residues(data, "chain2")
```

PDBePISA consistently uses `"chain1"` and `"chain2"` as keys inside `interfacing_residues`, but always verify this before running:

```python
import json
with open("your_interface.json") as f:
    data = json.load(f)

print("Top-level keys:       ", list(data.keys()))
print("Interfacing res keys: ", list(data["interfacing_residues"].keys()))
print("Interfacing bond keys:", list(data["interfacing_bonds"].keys()))
```

Expected output:
```
Top-level keys:        ['interface', 'interfacing_residues', 'interfacing_bonds']
Interfacing res keys:  ['chain1', 'chain2']
Interfacing bond keys: ['hbonds', 'salt_bridges']
```

#### 3. Update the chain labels in the plotting functions

The chain labels `"Chain F (TNF-α)"` and `"Chain L (TNFR1)"` are hardcoded strings in the plotting functions. The script will still run and produce correct data, but the figures will display the wrong biological names. Update these two lines at the top of `main()` to match your complex:

```python
# In main(), update these labels to match your complex:
CHAIN1_LABEL = "Chain A (Protein X)"   # e.g., "Chain A (IL-6)"
CHAIN2_LABEL = "Chain B (Receptor Y)"  # e.g., "Chain B (IL-6R)"
```

Then pass them to the plotting functions instead of the hardcoded strings.

---

### The One Code Change Required: BSA Normalization

This is the **most important fix** for portability. The script currently hardcodes:

```python
BSA_NORM = 140.0
```

This value (140 Å²) was the maximum BSA observed among all interface residues in the 8ZUI dataset. For a different complex, the maximum BSA may be higher or lower. If you leave it at 140 and your new complex has a residue with BSA > 140 Å², the normalized score exceeds 1.0 and the hotspot ranking is distorted.

**Fix — make BSA normalization data-driven:**

In `compute_hotspot_scores()`, replace the hardcoded constant with a runtime-computed value:

```python
def compute_hotspot_scores(df, bsa_norm=None):
    """
    Compute composite hotspot scores.
    If bsa_norm is None, it is computed automatically as the maximum BSA
    in the dataset, making the normalization fully data-driven.
    """
    df = df.copy()
    if bsa_norm is None:
        bsa_norm = df["bsa"].max()   # Data-driven normalization
    df["bond_bonus"]    = df["bond_type"].map(BOND_BONUS).fillna(0.0)
    df["hotspot_score"] = (
        df["bsa"] / bsa_norm
        + df["solv_en"].abs()
        + df["bond_bonus"]
    )
    df = df.sort_values("hotspot_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df
```

And remove the global constant from the top of the script:

```python
# Remove this line:
# BSA_NORM = 140.0
```

With this single change, the BSA normalization automatically adapts to any complex.

---

### Portability Summary

| Scenario | Can you run unchanged? | Action required |
|---|---|---|
| Same PDB, different interface index | No | Re-export the correct interface JSON from PDBePISA |
| Different PDB, same chain key structure | Yes* | Fix `BSA_NORM` (see above) + update chain labels |
| Homo-oligomeric complex (symmetric chains) | Yes* | Symmetric results are expected and meaningful |
| Hetero-complex with >2 chains at the interface | Partial | Run on each pairwise interface JSON separately |
| Complex where max BSA > 140 Å² | No | Fix `BSA_NORM` (see above) — required |

*After applying the `BSA_NORM` fix and updating chain labels.

---

### Quick-Start Checklist for a New Complex

Before running the pipeline on a new PDB complex, work through this checklist:

- [ ] Retrieved the PDB structure and identified the biologically relevant chain pair
- [ ] Opened PDBePISA, located the correct interface index for that chain pair
- [ ] Exported the JSON for that specific interface (not the full structure JSON)
- [ ] Confirmed JSON keys: `interfacing_residues` → `chain1` / `chain2`; `interfacing_bonds` → `hbonds` / `salt_bridges`
- [ ] Removed the hardcoded `BSA_NORM = 140.0` constant and switched to `bsa_norm=None` (data-driven)
- [ ] Updated chain label strings in the plotting functions to match the new complex
- [ ] Run the script and verify that the number of interface residues and total BSA are biologically reasonable (BSA per chain > 500 Å² for a genuine PPI)
