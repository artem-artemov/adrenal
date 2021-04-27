# Human adrenal gland paper

## Notebooks

- Whole dataset
  * [Whole human dataset - 9 samples](notebooks/human_whole_dataset.ipynb)

- Adrenal medulla
  * [Human adrenal medulla](notebooks/human_adrenal_medulla.ipynb)
  * [RNA velocity](notebooks/human_adrenal_medulla_velocity.ipynb)
  * [Genes specific to bridge population - deviation from linear mix](notebooks/deviation_from_linear_mix.ipynb)


- Other fates
  * [Adrenal cortex, endothelium, blood](notebooks/human_cortex_endothelium_blood.ipynb)

- Mouse trunk dataset
  * [Mouse: whole dataset, medulla, trajectories](notebooks/adrenal_mouse.ipynb)
  * [RNA velocity](notebooks/mouse_adrenal_medulla_velocity.ipynb)

## Download data

```sh
mkdir data/Seurat
wget -P data/Seurat/ http://pklab.med.harvard.edu/artem/adrenal/data/Seurat/adrenal.human.seurat.scrublet.rds
wget -P data/Seurat/ http://pklab.med.harvard.edu/artem/adrenal/data/Seurat/adrenal.mouse.seurat.rds
```

Alternatively, all preprocessed data including h5 files (read counts from cellranger) can be downloaded from GEO: 
- Human: GSE147821 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147821
- Mouse: GSE166262 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166262

