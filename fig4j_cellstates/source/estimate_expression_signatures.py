import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os

# Using cell2location
import cell2location
import scvi
import torch

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.backends.backend_pdf as mpdf

# Silence scanpy since it prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')
data_type = 'float32'

# Global settings for CUDA
seed = 42
torch.manual_seed(seed)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

from config import lung_config
inDir = "scRNAseq/"
outDir = "cell2location/"
def main(environment="Tumour"):
    environment_name = environment if environment=="Tumour" else "Healthy_Background"
    infile = "10X_Lung_{0}_Annotated_v2_TYMSTK1_classification.h5ad".format(environment_name)
    annotations = lung_config[environment]["name"]
    to_remove = lung_config[environment]["remove"]
    adata_sc_raw = sc.read_h5ad(inDir+infile)
    adata_sc_raw = adata_sc_raw[~adata_sc_raw.obs[annotations].isin(to_remove),:]
    if environment != "Tumour":
        adata_sc_raw = adata_sc_raw[adata_sc_raw.obs["environment"]==environment,:] # Select on the Healthy/Background tissue
    adata_sc_raw.obs["Cell types"] = adata_sc_raw.obs[annotations]
    adata_sc_raw = adata_sc_raw[~adata_sc_raw.obs["Cell types"].isna(), :]
    adata_sc_raw.var["Gene"] = adata_sc_raw.var.index # Genes as var name columns

    # Recommended NB regression configuration with scvi-tools and cell2location external models
    cell2location.models.RegressionModel.setup_anndata(adata=adata_sc_raw,
                            batch_key="patient_sample",
                            labels_key="Cell types",
                            categorical_covariate_keys=None
                            )
    from cell2location.models import RegressionModel
    mod = RegressionModel(adata_sc_raw)
    epochs = 500 if environment=="Healthy" else 100
    mod.train(max_epochs=epochs, batch_size=2500, train_size=1, lr=0.002) # Note that the validation set is not implemented for this model (so setting train_size=1)
    mod.save(outDir+"regression_model/{0}".format(environment), overwrite=True)
    adata_sc_raw.write(outDir+"regression_model/{0}/sc.h5ad".format(environment))

    train_performance = mod.history["elbo_train"]
    train_performance.plot(logy=True)
    plt.xlabel("Training epochs")
    plt.ylabel("-ELBO loss")
    plt.title(f"Lung ({environment})")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(outDir+"regression_model/{0}/ELBO_Lung_{0}_Regression.png".format(environment))
    plt.close()

if __name__=="__main__":
    # Adding an argument parser
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-e", "--environment", type=str, help="Tumour, Background, Healthy",default="Tumour")

    options = parser.parse_args()

    # Defining dictionary to be passed to the main function
    option_dict = dict( (k, v) for k, v in vars(options).items() if v is not None)
    print(option_dict)
    main(**option_dict)
