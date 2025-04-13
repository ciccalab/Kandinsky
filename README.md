# Kandinsky <img src="man/figures/logo.png" align="right" height="138" alt="" />
Spatial analysis toolkit expanding Seurat functionalities for Spatial transcriptomic and proteomic data.

To speed up package installation in R, user can first create a conda environment starting from the environment.yaml file.
This will make available most of Kandinsky dependencies before downloading it.

For conda users:
```
conda env create -f environment.yaml
```

For mamba users:
```
mamba env create -f environment.yaml
```

After setting up the conda environment, you can activate the environment and start a new R session to download Kandinsky:

```
conda activate kandinsky
R
```

Install development version from GitHub:
```r
install.packages('devtools')
devtools::install_github('ciccalab/Kandinsky')
```
