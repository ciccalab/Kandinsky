# Kandinsky

## Overview <img src="man/figures/logo.png" align="right" height="138" alt="" />
Kandinsky is a spatial analysis toolkit designed to provide a compendium of methods to study cell and spot spatial neighbourhoods using spatial transcriptomic and proteomic data.

<img src="man/figures/Kandinsky_Overview.png" align="center" height="500" alt="" />


As an input, Kandinsky uses gene or protein expression values and cell or spot coordinates deriving from any spatial transcriptomic or proteomic platform and implements helper functions to automate their loading and formatting into a Seurat object.\
Starting from their coordinates, Kandinsky groups cells or spots into neighbourhoods (c/s-NBs) according to their reciprocal spatial relationships inferred with five methods:\
	- K-nearest neighbours (KNN);\
	- Cell/spot centroid distance;\
	- Delaunay triangulation;\
	- Queen contiguity;\
	- cell membrane distance. 

KNN, centroid distance, and Delaunay triangulation are applicable to both spot and cell data, while Queen contiguity is limited to spot data and membrane distances can be measured only from single cell segmentation data. 

Once defined, c/s-NBs can be used to derive clusters with similar composition, measure spatial co-localisation or dispersion and infer hot and cold gene expression areas.

## Setting up Kandinsky environment
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


## Installation
After setting up the conda environment, you can activate the environment and start a new R session to install Kandinsky:

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

```
cd path/to/local/folder
git clone https://github.com/ciccalab/Kandinsky.git
#Alternative code for SSH connection: git clone git@github.com:ciccalab/Kandinsky.git
```

After downloading the repository, change directory to `Kandinsky/` and use `devtools` to install the package:

```
cd Kandinsky/
R
```
```r
devtools::install()
```
