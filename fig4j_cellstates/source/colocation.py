#!/usr/bin/env python3
import os
import numpy as np
import scanpy as sc
from scipy import sparse
import torch


def compute_total_gene_counts(adata):
    X = adata.X
    if sparse.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel()
    return X.sum(axis=1)


def ensure_spatial_has_hires(adata):
    """
    cell2location.run_colocation internally calls scanpy sc.pl.spatial, which often defaults to img_key='hires'.
    If your merged h5ad was created with only lowres images, plotting will fail with KeyError: 'hires'.

    This function makes plotting robust by:
      - Creating images['hires'] from images['lowres'] if missing
      - Setting tissue_hires_scalef from tissue_lowres_scalef if missing
    """
    if "spatial" not in adata.uns or len(adata.uns["spatial"]) == 0:
        raise KeyError("adata.uns['spatial'] missing or empty (no Visium images attached).")

    for lib_key in list(adata.uns["spatial"].keys()):
        sd = adata.uns["spatial"][lib_key]

        # images dict
        images = sd.get("images", {})
        if "hires" not in images:
            if "lowres" in images:
                images["hires"] = images["lowres"]
            else:
                raise KeyError(f"Library '{lib_key}' has no 'lowres' or 'hires' in uns['spatial'][...]['images'].")

        # scalefactors dict
        scalefactors = sd.get("scalefactors", {})
        if "tissue_hires_scalef" not in scalefactors:
            if "tissue_lowres_scalef" in scalefactors:
                scalefactors["tissue_hires_scalef"] = scalefactors["tissue_lowres_scalef"]
            else:
                # If neither exists, scanpy plotting will not scale properly
                raise KeyError(
                    f"Library '{lib_key}' missing scalefactors tissue_hires_scalef and tissue_lowres_scalef."
                )

        sd["images"] = images
        sd["scalefactors"] = scalefactors
        adata.uns["spatial"][lib_key] = sd


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Merged h5ad with images (e.g. sp_all_samples_merged_with_images.h5ad)")
    parser.add_argument("-o", "--outdir", default="coloc", help="Base output directory (default: coloc)")
    parser.add_argument("-e", "--environment", default="Tumour", help="Label for outputs (default: Tumour)")
    parser.add_argument("--min_counts", type=int, default=800, help="Filter spots by total_gene_counts >= min_counts (set 0 to disable)")
    parser.add_argument("--n_fact_min", type=int, default=7, help="Min number of NMF factors (inclusive)")
    parser.add_argument("--n_fact_max", type=int, default=10, help="Max number of NMF factors (inclusive)")
    parser.add_argument("--n_restarts", type=int, default=5, help="Number of training restarts")
    args = parser.parse_args()

    # Reproducibility
    seed = 42
    torch.manual_seed(seed)

    print(f"Loading: {args.input}")
    adata_vis = sc.read_h5ad(args.input)

    # Critical: prevent KeyError 'hires' during plotting inside run_colocation
    ensure_spatial_has_hires(adata_vis)

    # Optional filter
    if args.min_counts and args.min_counts > 0:
        print("Computing total_gene_counts...")
        adata_vis.obs["total_gene_counts"] = compute_total_gene_counts(adata_vis)
        n0 = adata_vis.n_obs
        adata_vis = adata_vis[adata_vis.obs["total_gene_counts"] >= args.min_counts].copy()
        n1 = adata_vis.n_obs
        print(f"Filtered spots on total_gene_counts ≥ {args.min_counts}: {n0} → {n1}")
    else:
        print("No filtering applied (min_counts=0).")

    # Run colocation
    from cell2location import run_colocation

    n_fact = np.arange(args.n_fact_min, args.n_fact_max + 1)
    export_path = os.path.join(args.outdir, args.environment, "colocation_model")
    os.makedirs(export_path, exist_ok=True)

    print(f"Running run_colocation with n_fact={list(n_fact)}")
    print(f"Export path: {export_path}")

    res_dict, adata_out = run_colocation(
        adata_vis,
        model_name="CoLocatedGroupsSklearnNMF",
        train_args={
            "n_fact": n_fact,
            "sample_name_col": "sample",
            "n_restarts": args.n_restarts,
        },
        export_args={"path": export_path},
    )

    print("Done.")
    return res_dict, adata_out


if __name__ == "__main__":
    main()
