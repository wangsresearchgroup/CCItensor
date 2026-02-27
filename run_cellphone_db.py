#!/usr/bin/env python3

import os
import sys
import gc
import traceback

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import liana as li
import matplotlib.pyplot as plt
import seaborn as sns
import mygene

from DMWT import mDWT


# =====================================================
# Notebook-faithful preprocessing
# =====================================================

def filter_cells(adata, min_genes, max_genes):
    """
    EXACT match to Jupyter notebook behavior.
    Operates in-place.
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)


def log_norm(adata):
    """
    EXACT match to Jupyter notebook behavior.
    """
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
        exclude_highly_expressed=True
    )
    sc.pp.log1p(adata)


def preprocess(adata, min_genes=200, max_genes=2500):
    """
    EXACT notebook preprocessing pipeline (in-place).
    """
    filter_cells(adata, min_genes, max_genes)
    log_norm(adata)


# =====================================================
# Gene ID conversion (CRITICAL for CellPhoneDB)
# =====================================================

def convert_ensembl_to_hgnc(adata):
    """
    Convert Ensembl gene IDs → HGNC symbols (human).
    Drops unmapped + duplicated genes.
    """
    print("  • Converting Ensembl → HGNC", flush=True)

    mg = mygene.MyGeneInfo()
    q = mg.querymany(
        adata.var_names.tolist(),
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=True,
        verbose=False
    )

    q = q[~q.index.duplicated(keep="first")]
    q = q[q["symbol"].notna()]

    adata = adata[:, adata.var_names.isin(q.index)].copy()
    adata.var["hgnc"] = q.reindex(adata.var_names)["symbol"].values
    adata.var_names = adata.var["hgnc"]

    adata = adata[:, ~adata.var_names.duplicated()].copy()

    print(
        f"    Genes after HGNC mapping: {adata.n_vars}",
        flush=True
    )

    return adata


# =====================================================
# Utilities
# =====================================================

def sparse_relu(X):
    X.data[X.data < 0] = 0.0
    X.eliminate_zeros()
    return X


def generate_plots(cpdb_df, donor_id, layer_name, output_path):
    # ---- Heatmap ----
    heat = (
        cpdb_df
        .groupby(["source", "target"])["lr_means"]
        .mean()
        .unstack(fill_value=0)
    )

    plt.figure(figsize=(10, 8))
    sns.heatmap(heat, cmap="viridis", square=True)
    plt.title(f"Interactions: {donor_id} ({layer_name})")
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f"{layer_name}_heatmap.png"))
    plt.close()

    # ---- Dot plot (Top 20 LR) ----
    top_lr = (
        cpdb_df
        .groupby(["ligand", "receptor"])["lr_means"]
        .mean()
        .sort_values(ascending=False)
        .head(20)
        .index
    )

    if len(top_lr) > 0:
        lr_df = (
            cpdb_df
            .set_index(["ligand", "receptor"])
            .loc[top_lr]
            .reset_index()
        )

        lr_df["lr_pair"] = lr_df["ligand"] + " → " + lr_df["receptor"]

        plt.figure(figsize=(10, 10))
        sns.scatterplot(
            data=lr_df,
            x="source",
            y="lr_pair",
            size="lr_means",
            hue="lr_means",
            sizes=(20, 200),
            palette="viridis"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"{layer_name}_dotplot.png"))
        plt.close()


# =====================================================
# Core processing
# =====================================================

def process_single_h5ad(h5ad_path, results_base):
    if "/LUAD/" in h5ad_path:
        cancer = "LUAD"
    elif "/LUSC/" in h5ad_path:
        cancer = "LUSC"
    else:
        cancer = "UNKNOWN"

    donor_id = os.path.basename(h5ad_path).replace(".h5ad", "")
    donor_out = os.path.join(results_base, cancer, donor_id)
    os.makedirs(donor_out, exist_ok=True)

    

    print(f"▶ Processing donor: {donor_id}", flush=True)

    # ---- Load ----
    adata = ad.read_h5ad(h5ad_path)

    # ---- Notebook preprocessing (EXACT) ----
    preprocess(
        adata,
        min_genes=200,
        max_genes=2500
    )

    # ---- Ensembl → HGNC ----
    adata = convert_ensembl_to_hgnc(adata)

    # ---- DMWT ----
    print("  • Running DMWT (m=4)", flush=True)
    dmwt_results = mDWT(adata.X, 4)

    for name, comp in dmwt_results.items():
        adata.layers[name] = sparse_relu(comp)

    # ---- Explicit CellPhoneDB resource ----
    li.resource.select_resource("CellPhoneDB")

    # ---- CellPhoneDB per layer ----
    layers_to_analyze = [None, "low", "high1", "high2", "high3"]

    for layer in layers_to_analyze:
        label = layer if layer is not None else "standard"
        print(f"  • CellPhoneDB on layer: {label}", flush=True)

        cpdb_res = li.mt.cellphonedb(
            adata,
            layer=layer,
            groupby="cell_type_predicted",
            use_raw=False,
            inplace=False
        )

        cpdb_res.to_csv(
            os.path.join(donor_out, f"{label}_results.csv"),
            index=False
        )

        generate_plots(cpdb_res, donor_id, label, donor_out)

        del cpdb_res
        gc.collect()

    # ---- Save processed AnnData ----
    adata.write(os.path.join(donor_out, f"{donor_id}_processed.h5ad"))

    del adata
    gc.collect()

    print(f"✔ Finished donor: {donor_id}", flush=True)


# =====================================================
# Entry point
# =====================================================

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print(
            "Usage:\n"
            "  python run_dmwt_cellphonedb.py <input.h5ad> <results_base_dir>",
            file=sys.stderr
        )
        sys.exit(1)

    try:
        process_single_h5ad(sys.argv[1], sys.argv[2])
    except Exception:
        print("❌ Fatal error during processing", file=sys.stderr)
        traceback.print_exc()
        sys.exit(2)
