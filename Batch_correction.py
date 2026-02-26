#!/usr/bin/env python 

from pathlib import Path
import numpy as np
import pandas as pd

import scvi
import torch

import pickle
import argparse
import sys
import os

import scanpy as sc
import anndata

import logging

logger = logging.getLogger('scvi')
logger.setLevel('INFO')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh = logging.StreamHandler()
sh.setFormatter(formatter)
logger.addHandler(sh)

## settings 
save_dir = '/faststorage/project/TCRseq_AZ/Analysis/Nanna/SingleCell/Jupyter/h5ad/PARSE_10X/10X/'
p = Path(save_dir)
p.mkdir(parents=True, exist_ok=True)


## load adata
adata = sc.read('/faststorage/project/TCRseq_AZ/Analysis/Nanna/SingleCell/Jupyter/h5ad/10X_postQC_Tcells.h5ad')
adata.layers['counts'] = adata.X.copy()

## run scvi
logger.info(f'adata: {adata.shape}')

scvi.model.SCVI.setup_anndata(adata, layer="counts",  batch_key="sample")

## this model definition is based off HLCA atlas integration

model = scvi.model.SCVI(adata,
      n_latent=int(10), # 30 # 64
      n_hidden=int(128), # 32 # 128
      n_layers=int(1), # 3 # 5
      dropout_rate=0.1,
      dispersion='gene')
      #gene_likelihood='nb')

vae_epochs = 500

early_stopping_kwargs = {
    'early_stopping': True, 
    'early_stopping_monitor': 'elbo_validation',
    'early_stopping_patience': 10,
    'early_stopping_min_delta': 0.0,
}

plan_kwargs = {
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

scvi.settings.dl_num_workers = 7

# Train scVI model
model.train(
    max_epochs=vae_epochs,
    plan_kwargs = plan_kwargs,
    **early_stopping_kwargs,
    check_val_every_n_epoch=1
)

latent = model.get_latent_representation(adata)

logger.info(f'Saving to {save_dir}')
model.save(f'{save_dir}/snapshot')
np.save(f'{save_dir}/latent_space.npy', latent)
barcodes = np.array(adata.obs_names)
np.save(f'{save_dir}/barcodes.npy', barcodes)