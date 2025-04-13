# Kandinsky <img src="man/figures/logo.png" align="right" height="138" alt="" />

##Overview
Kandinsky is a spatial analysis toolkit expanding Seurat functionalities for cell and spot neighbourhood analysis using spatial transcriptomic and proteomic data.

## Installation
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

Kandinsky can be downloaded from Github and installed using `R` package `devtools`:
```r
install.packages('devtools')
devtools::install_github('ciccalab/Kandinsky')
```

Alternatively, install Kandinsky by downloading the repository on your local computer using `git`:

```r
cd path/to/local/folder
git clone https://github.com/ciccalab/Kandinsky.git
#Alternative code for SSH connection: git clone git@github.com:ciccalab/Kandinsky.git
```

After downloading the repository, change directory to `Kandinsky/` and use `devtools` to install the package:

```r
cd Kandinsky/
R
devtools::install()
```
