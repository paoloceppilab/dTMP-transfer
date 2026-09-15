#!/usr/bin/env python3
"""
Automatically aligns and merges ALL samples found in a Cell2Location posterior file
with their corresponding Visium histology images, producing a SINGLE merged .h5ad
that can be plotted with scanpy sc.pl.spatial without library_id guessing.

Each slide's `uns['spatial']` key is set to the sample name (e.g. `P10_T1`),
so later plotting can use `library_id=sample_name` deterministically.

Usage:
  python3 source/save_all_samples_merged.py \
    -i cell2location/spatial_model/Tumour/sp.h5ad \
    --spatial_dir SpatialDatasets_prepared \
    -o . \
    -f sp_all_samples_merged_with_images.h5ad

Optional:
  --min_counts 800        # apply total_gene_counts filter
  --keep_img lowres       # keep only lowres image to reduce file size
"""

import os
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from visium_qc_and_visualisation import read_and_qc


# --------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------

def compute_total_gene_counts(adata):
    """Return total counts per spot from adata.X (robust to sparse)."""
    X = adata.X
    if sparse.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel()
    else:
        return X.sum(axis=1)


def force_spatial_key_to_sample(slide_vis, sample_name):
    """
    Ensure slide_vis.uns['spatial'] has exactly one key and that key == sample_name.
    This prevents library_id mismatches after concatenation.
    """
    if "spatial" not in slide_vis.uns or len(slide_vis.uns["spatial"]) == 0:
        raise KeyError("slide_vis.uns['spatial'] missing or empty.")

    old_key = list(slide_vis.uns["spatial"].keys())[0]
    if old_key != sample_name:
        slide_vis.uns["spatial"] = {sample_name: slide_vis.uns["spatial"][old_key]}
    return slide_vis


def optionally_keep_one_image(slide_vis, keep_img="lowres"):
    """
    Reduce file size by keeping only one image resolution in uns['spatial'].
    keep_img: 'lowres', 'hires', or 'both'
    """
    if keep_img == "both":
        return slide_vis

    lib_key = list(slide_vis.uns["spatial"].keys())[0]
    sd = slide_vis.uns["spatial"][lib_key]

    images = sd.get("images", {})
    if keep_img == "lowres":
        images.pop("hires", None)
    elif keep_img == "hires":
        images.pop("lowres", None)
    else:
        raise ValueError("--keep_img must be one of: lowres, hires, both")

    sd["images"] = images
    slide_vis.uns["spatial"][lib_key] = sd
    return slide_vis


def align_visium_and_posterior(adata_post, sample_name, spatial_dir, min_counts=0, keep_img="lowres"):
    """
    Align posterior data with Visium data for a single sample.

    Returns:
      slide_post: posterior restricted to this sample and common barcodes,
                  with spatial info attached and spatial key renamed to sample_name.
    """
    # 1) Load Visium (images + scalefactors + coords + counts)
    slide_vis = read_and_qc(sample_name, path=spatial_dir)

    # 2) Set the spatial key to the sample name.
    slide_vis = force_spatial_key_to_sample(slide_vis, sample_name)

    # 3) Optionally reduce size by keeping only lowres/hires image
    slide_vis = optionally_keep_one_image(slide_vis, keep_img=keep_img)

    # 4) Filter posterior to just this sample
    if "sample" in adata_post.obs:
        mask = adata_post.obs["sample"].astype(str) == str(sample_name)
    else:
        mask = adata_post.obs_names.astype(str).str.contains(str(sample_name))

    slide_post = adata_post[mask].copy()

    # 5) Find common barcodes (exact overlap)
    common_barcodes = slide_post.obs_names.intersection(slide_vis.obs_names)
    if len(common_barcodes) == 0:
        raise ValueError(f"No overlapping barcodes found for {sample_name}")

    # subset both
    slide_post = slide_post[common_barcodes].copy()
    slide_vis = slide_vis[common_barcodes].copy()

    # 6) Optional QC filter using Visium counts
    #    (this matches your ">=800 total_gene_counts" pattern)
    if min_counts > 0:
        slide_vis.obs["total_gene_counts"] = compute_total_gene_counts(slide_vis)
        keep = slide_vis.obs["total_gene_counts"].values >= min_counts

        n0 = slide_vis.n_obs
        slide_vis = slide_vis[keep].copy()
        slide_post = slide_post[keep].copy()
        n1 = slide_vis.n_obs
        print(f"    QC filter total_gene_counts >= {min_counts}: {n0} -> {n1}")

    # 7) Transfer spatial metadata to posterior slice
    slide_post.obsm["spatial"] = slide_vis.obsm["spatial"].copy()
    slide_post.uns["spatial"] = slide_vis.uns["spatial"].copy()

    # 8) If Visium has in_tissue annotation and posterior doesn't, transfer it
    for col in ["in_tissue", "array_row", "array_col"]:
        if col in slide_vis.obs.columns and col not in slide_post.obs.columns:
            slide_post.obs[col] = slide_vis.obs[col].values

    # ensure sample column is correct and string
    slide_post.obs["sample"] = str(sample_name)

    return slide_post


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Merge ALL samples from a posterior file with their spatial images (keyed by sample name)."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input posterior .h5ad file (e.g. sp.h5ad)"
    )
    parser.add_argument(
        "--spatial_dir",
        default="SpatialDatasets_prepared",
        help="Folder containing per-sample Visium folders (default: SpatialDatasets_prepared)."
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Directory to save the output file."
    )
    parser.add_argument(
        "-f", "--filename",
        default="sp_all_samples_merged_with_images.h5ad",
        help="Output filename."
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        default=0,
        help="If >0, filter spots by total_gene_counts >= this value before merging (recommended: 800)."
    )
    parser.add_argument(
        "--keep_img",
        choices=["lowres", "hires", "both"],
        default="lowres",
        help="Keep only lowres/hires/both images to control output size (default: lowres)."
    )

    args = parser.parse_args()

    # 1) Load posterior
    print(f"Loading input file: {args.input}")
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    adata_post = sc.read_h5ad(args.input)

    # 2) Detect samples
    if "sample" not in adata_post.obs:
        raise KeyError("Input file must have a 'sample' column in .obs to identify slides.")

    all_samples = [str(s) for s in adata_post.obs["sample"].unique().tolist()]
    all_samples = sorted(all_samples)

    print(f"Detected {len(all_samples)} samples: {all_samples}")
    print(f"Using spatial_dir: {args.spatial_dir}")
    print(f"Using keep_img: {args.keep_img}")
    if args.min_counts > 0:
        print(f"Applying QC filter: total_gene_counts >= {args.min_counts}")

    # 3) Align each sample
    aligned_slides = []
    skipped_samples = []

    print("\n--- Starting Alignment ---")
    for s in all_samples:
        try:
            print(f"Processing: {s}")
            slide = align_visium_and_posterior(
                adata_post,
                s,
                spatial_dir=args.spatial_dir,
                min_counts=args.min_counts,
                keep_img=args.keep_img,
            )
            aligned_slides.append(slide)

        except Exception as e:
            print(f"  [WARNING] Skipping {s}. Reason: {e}")
            skipped_samples.append(s)

    if not aligned_slides:
        print("No samples were successfully processed. Exiting.")
        return

    # 4) Concatenate (keep original spot barcodes; do not append extra suffixes)
    print(f"\n--- Concatenating {len(aligned_slides)} samples ---")
    merged_adata = sc.concat(
        aligned_slides,
        join="outer",
        merge="same",
        label=None,
        index_unique=None,
    )

    # 5) Rebuild merged uns['spatial'] with keys == sample names
    merged_adata.uns["spatial"] = {}
    for slide in aligned_slides:
        merged_adata.uns["spatial"].update(slide.uns["spatial"])

    # 6) Preserve cell2location metadata if present
    if "mod" in adata_post.uns:
        merged_adata.uns["mod"] = adata_post.uns["mod"]

    # 7) Save
    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, args.filename)

    print(f"Saving to: {out_path}")
    merged_adata.write(out_path)

    print("\nDone!")
    if skipped_samples:
        print(f"Note: The following samples were skipped due to errors: {skipped_samples}")


if __name__ == "__main__":
    main()
