# this script is about calculating the significance of shared rearrangements in drosera regia and drosera capensis , we will use the results of our genespace run in wd

# imports
import os
import csv
from collections import defaultdict
from __future__ import annotations
import argparse
from dataclasses import dataclass
from typing import Dict, List, Set, FrozenSet, Tuple
import pandas as pd
import gzip
# files
wd="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/results"
n0_tsv=os.path.join(wd, "N0.tsv")
comb_bed_txt=os.path.join(wd, "combBed.txt")
pangenes_dir = os.path.normpath(os.path.join(wd,"..","pangenes"))
sequence_ids_txt=os.path.join(wd, "SequenceIDs.txt")
pangenes_dir = os.path.normpath(os.path.join(wd, "..", "pangenes"))
species_list = ["A_abbreviatus", "D_capensis", "D_regia", "Dio_muscipula", "N_gra_dom", "T_peltatum"]



def load_pangenes_df(pangenes_gz_path: str) -> pd.DataFrame:
    with gzip.open(pangenes_gz_path, "rt") as f:
        df = pd.read_csv(f, sep="\t", dtype={"genome": str, "chr": str})
    return df

def build_pg_adjacencies_native(df: pd.DataFrame) -> set:
    """
    Build adjacency set in pgID space using the species' native coordinates (chr, ord).
    """
    # numeric ord + pgID
    df = df.copy()
    df["ord"] = pd.to_numeric(df["ord"], errors="coerce")
    df["pgID"] = pd.to_numeric(df["pgID"], errors="coerce")
    df = df.dropna(subset=["chr", "ord", "pgID"])
    df["ord"] = df["ord"].astype(int)
    df["pgID"] = df["pgID"].astype(int)

    # sort by native order
    df = df.sort_values(["chr", "ord"], kind="mergesort")

    adj = set()
    for chr_, sub in df.groupby("chr", sort=False):
        pg = sub["pgID"].tolist()
        for i in range(len(pg) - 1):
            a, b = pg[i], pg[i + 1]
            if a == b:
                continue
            adj.add(frozenset((a, b)))
    return adj


def load_species_pg_adjacencies_from_pangenes(pangenes_dir: str, species: str) -> set:
    """
    Each *species*_pangenes.txt.gz file contains rows for all genomes.
    We load the file and filter to genome==species.
    """
    path = os.path.join(pangenes_dir, f"{species}_pangenes.txt.gz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing: {path}")

    df = load_pangenes_df(path)

    # FILTER TO THIS SPECIES ONLY  <-- the crucial fix
    df = df[df["genome"] == species]

    # optional quality filter
    if "flag" in df.columns:
        df = df[df["flag"] == "PASS"]

    return build_pg_adjacencies_native(df)

def compute_distances_relative_to_reference(ref_adj: set, species_adj: dict) -> pd.DataFrame:
    rows = []
    for sp, adj in species_adj.items():
        preserved = len(ref_adj & adj)
        broken = len(ref_adj - adj)
        rows.append({
            "species": sp,
            "preserved": preserved,
            "broken": broken,
            "prop_broken": broken / len(ref_adj),
            "ref_adj_total": len(ref_adj),
            "species_adj_total": len(adj),
        })
    df = pd.DataFrame(rows).sort_values(["broken", "species"]).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df





# ---------------------------
# Step 1: build adjacency sets in integrated pgID space
# ---------------------------

species_adj_pg = {sp: load_species_pg_adjacencies_from_pangenes(pangenes_dir, sp) for sp in species_list}

# reference
ref_species = "N_gra_dom"
ref_adj = species_adj_pg[ref_species]

print("Reference adj size:", len(ref_adj))
for sp in species_list:
    print(sp, "adj size:", len(species_adj_pg[sp]))

dist_df = compute_distances_relative_to_reference(ref_adj, species_adj_pg)
print(dist_df.to_string(index=False))

ref_species = "N_gra_dom"
ref_adj = species_adj_pg[ref_species]

dist_df = compute_distances_relative_to_reference(ref_adj, species_adj_pg)
print(dist_df.to_string(index=False))



# okay that didnt work













# choose ancestral-like reference (as you stated): N_gra_dom
ref_species = "N_gra_dom"
if ref_species not in species_adj_pg:
    raise ValueError(f"Reference {ref_species} not found. Available: {sorted(species_adj_pg.keys())}")

ref_adj = species_adj_pg[ref_species]

print("Reference species:", ref_species)
print("Reference adjacencies:", len(ref_adj))
print("Species adjacency counts (pgID-space):")
for sp in sorted(species_adj_pg.keys()):
    print(" ", sp, len(species_adj_pg[sp]))

# ---------------------------
# Step 2: compute distances + rank species
# ---------------------------

dist_df = compute_distances_relative_to_reference(ref_adj, species_adj_pg)
print("\nRanked distances relative to", ref_species)
print(dist_df.to_string(index=False))

out_tsv = os.path.join(wd, f"distance_to_{ref_species}_by_pgID.tsv")
dist_df.to_csv(out_tsv, sep="\t", index=False)
print("\nWrote:", out_tsv)

# (Later) you'll use species_adj_pg + ref_adj to define breaks and shared breaks for D_regia vs D_capensis.

