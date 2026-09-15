import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc

data_type = 'float32'

# This line forces theano to use the GPU and should go before importing cell2location
#os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=' + data_type + ',force_device=True'
# If using the CPU uncomment this:
os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,openmp=True,force_device=True'

import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.backends.backend_pdf as mpdf

# Silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')

inDir = "/nfs/research1/gerstung/jp27/AC_LNG/3_spaceranger_output/"
outDir = "/nfs/research1/gerstung/nelson/outputs/lung/"

sample_names = [
        'spaceranger110_count_36209_OTAR_LNGsp9476038_GRCh38-2020-A',
        'spaceranger110_count_36209_OTAR_LNGsp9476040_GRCh38-2020-A',
        'spaceranger110_count_36210_OTAR_LNGsp9476042_GRCh38-2020-A',
        'spaceranger110_count_36210_OTAR_LNGsp9476044_GRCh38-2020-A',
        'spaceranger110_count_36209_OTAR_LNGsp9476039_GRCh38-2020-A',
        'spaceranger110_count_36209_OTAR_LNGsp9476041_GRCh38-2020-A',
        'spaceranger110_count_36210_OTAR_LNGsp9476043_GRCh38-2020-A',
        'spaceranger110_count_36210_OTAR_LNGsp9476045_GRCh38-2020-A'
]

def read_and_qc(sample_name, path=inDir):
    """
    Read one 10X Visium experiment into an AnnData object and compute QC metrics.

    Supports two layouts:
    1) Original AC_LNG layout:
         path/<sample_name>/filtered_feature_bc_matrix.h5
    2) Local flat layout (your current setup):
         path/<sample_name>-filtered_feature_bc_matrix.h5
    """
    # Ensure we use os.path.join instead of string concat
    folder_dir   = os.path.join(path, sample_name)  # e.g. SpatialDatasets/P10_T1
    h5_in_folder = os.path.join(folder_dir, "filtered_feature_bc_matrix.h5")
    h5_flat      = os.path.join(path, f"{sample_name}-filtered_feature_bc_matrix.h5")

    # ------------------------------------------------------------------
    # 1) Try original visium folder layout (has images, positions, etc.)
    # ------------------------------------------------------------------
    if os.path.exists(h5_in_folder):
        adata = sc.read_visium(
            folder_dir,
            count_file="filtered_feature_bc_matrix.h5",
            load_images=True
        )

    # ------------------------------------------------------------------
    # 2) Try local flat HDF5 layout (only counts; no images)
    # ------------------------------------------------------------------
    elif os.path.exists(h5_flat):
        # Read raw matrix from 10x h5
        adata = sc.read_10x_h5(h5_flat)

        # For compatibility with the rest of the pipeline,
        # create a minimal .uns["spatial"] entry so that
        # concatenation and cell2location work without images.
        adata.uns["spatial"] = {sample_name: {}}

    else:
        raise FileNotFoundError(
            f"Could not find Visium data for sample '{sample_name}'. "
            f"Tried:\n  {h5_in_folder}\n  {h5_flat}"
        )

    # ----------------------------------------------------------
    # Harmonise var / gene metadata (same as original script)
    # ----------------------------------------------------------
    adata.obs["sample"] = sample_name  # add sample name

    # Store gene symbols
    adata.var["SYMBOL"] = adata.var_names

    # If gene_ids present, use as ENSEMBL IDs and move to var_names
    if "gene_ids" in adata.var.columns:
        adata.var.rename(columns={"gene_ids": "ENSEMBL"}, inplace=True)
        adata.var_names = adata.var["ENSEMBL"]
        adata.var.drop(columns="ENSEMBL", inplace=True)

    # ----------------------------------------------------------
    # QC metrics
    # ----------------------------------------------------------
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Mito genes (make it robust to case)
    adata.var["mt"] = [str(g).lower().startswith("mt-") for g in adata.var["SYMBOL"]]
    adata.obs["mt_frac"] = (
        adata[:, adata.var["mt"].values].X.sum(1).A.squeeze() /
        adata.obs["total_counts"]
    )

    # Add sample name to obs index (as in original code)
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs_names = adata.obs["sample"] + "_" + adata.obs_names
    adata.obs.index.name = "spot_id"

    return adata


def select_slide(adata, s, s_col='sample'):
    """ This function selects the data for one slide from the spatial anndata object.
        :param adata: Anndata object with multiple spatial experiments
        :param s: name of selected experiment
        :param s_col: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}

    return slide

def main():
    pdf = mpdf.PdfPages("Visium_Lung_QC.pdf")
    # Pick up the scRNA-seq gene expression model
    reg_mod_name = 'RegressionGeneBackgroundCoverageTorch_40covariates_54574cells_1470genes'
    reg_path = f'{outDir}regression_model_sub_types/{reg_mod_name}/'
    # scRNA-seq reference (raw counts)
    #adata_sc_model = sc.read(f'{reg_path}sc.h5ad')

    # Pick up the visium data from SpaceRanger
    # Read the data into anndata objects
    slides = [read_and_qc(sample_name, path=inDir) for sample_name in sample_names]

    # Combine anndata objects together
    adata = slides[0].concatenate(slides[1:], batch_key="sample", uns_merge="unique", batch_categories=sample_names, index_unique=None)
    # Mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
    adata = adata[:, ~adata.var['mt'].values]

    fig, axs = plt.subplots(len(slides), 4, figsize=(15, 4*len(slides)-4))
    for i, s in enumerate(adata.obs['sample'].unique()):
        slide = select_slide(adata, s)
        sns.distplot(slide.obs['total_counts'], kde=False, ax = axs[i, 0])
        axs[i, 0].set_xlim(0, adata.obs['total_counts'].max())
        axs[i, 0].set_xlabel(f'total_counts | {s}')

        sns.distplot(slide.obs['total_counts'][slide.obs['total_counts']<20000], kde=False, bins=40, ax = axs[i, 1])
        axs[i, 1].set_xlim(0, 20000)
        axs[i, 1].set_xlabel(f'total_counts | {s}')

        sns.distplot(slide.obs['n_genes_by_counts'], kde=False, bins=60, ax = axs[i, 2])
        axs[i, 2].set_xlim(0, adata.obs['n_genes_by_counts'].max())
        axs[i, 2].set_xlabel(f'n_genes_by_counts | {s}')

        sns.distplot(slide.obs['n_genes_by_counts'][slide.obs['n_genes_by_counts']<6000], kde=False, bins=60, ax = axs[i, 3])
        axs[i, 3].set_xlim(0, 6000)
        axs[i, 3].set_xlabel(f'n_genes_by_counts | {s}')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    pdf.close()

    for sample in adata.obs['sample'].unique():
        slide = select_slide(adata, sample)
        pdf_vis = mpdf.PdfPages("Visium_{0}.pdf".format(sample))
        with mpl.rc_context({'figure.figsize': [6,7], 'axes.facecolor': 'white'}):
            sc.pl.spatial(slide,
                          img_key = "hires",
                          cmap='magma',
                          library_id=list(slide.uns['spatial'].keys())[0],
                          color=['total_counts', 'n_genes_by_counts'],
                          size=1,
                          gene_symbols='SYMBOL',
                          show=False,
                          return_fig=True)
            plt.tight_layout()
            pdf_vis.savefig()
            plt.close()

        with mpl.rc_context({'figure.figsize': [6,7],'axes.facecolor': 'black'}):
            sc.pl.spatial(slide,
                          color=["IL10", "FOXP3", "CTLA4"], # Change these gene names
                          img_key=None,
                          size=1,
                          vmin=0,
                          cmap='magma',
                          vmax='p99.0',
                          gene_symbols='SYMBOL',
                          show=False,
                          return_fig=True)
            plt.tight_layout()
            pdf_vis.savefig()
            plt.close()

        pdf_vis.close()

if __name__=="__main__":
    main()
