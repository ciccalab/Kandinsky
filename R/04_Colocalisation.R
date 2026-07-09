#' @title same-color (BB) and different-color (BW) join count test
#' @name multi_jc
#' @description
#' Readapted version of joincount.multi() function from the spdep package.
#'
#' @param fx a factor of the same length as the neighbours and weights objects in listw
#' @param listw a listw object
#' @param wc.type input type to use to extract summary weight constants used to calculate join count statistics. Must be one of the following: "listw", "mat"
#' @return data.frame containing join count coefficients for each cell type pair tested
#' @export
multi_jc = function(fx, listw, wc.type = c('mat','listw')){
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!is.factor(fx))
    stop(paste(deparse(substitute(fx)), "is not a factor"))
  if (any(is.na(fx)))
    stop("NA in factor")
  n <- length(listw$neighbours)
  if (n != length(fx))
    stop("objects of different length")
  cards <- spdep::card(listw$neighbours)
  ifx <- as.integer(fx)
  minf = min(as.numeric(table(fx)))
  if(minf<=2){
	  stop('each cell type must contains more than 2 cells')
  }
  k <- length(levels(fx))
  if (k < 2)
    stop("must be at least two levels in factor")
  sn <- spdep::listw2sn(listw)

  ###THIS SECTION CAN BE DONE ONLY ONCE###
  tab <- table(fx)
  ntab <- as.numeric(as.vector(tab))
  if(length(wc.type)>1){
    wc.type = NULL
  }
  wc.type = wc.type %||% 'mat'
  if(wc.type == 'mat'){
    wc <- spweights_mat(listw)
  }
  if(wc.type == 'listw'){
    wc <- spdep::spweights.constants(listw, zero.policy = T,adjust.n = T)
  }
  N <- wc$n
  S02 <- wc$S0 * wc$S0
  Ejc <- (wc$S0 * (ntab * (ntab - 1)))/(2 * N * wc$n1)
  Vjc <- (wc$S1 * (ntab * (ntab - 1)))/(N * wc$n1)
  Vjc <- Vjc + (((wc$S2 - 2 * wc$S1) * ntab * (ntab - 1) *
                   (ntab - 2))/(N * wc$n1 * wc$n2))
  Vjc <- Vjc + (((S02 + wc$S1 - wc$S2) * ntab * (ntab - 1) *
                   (ntab - 2) * (ntab - 3))/(N * wc$n1 * wc$n2 * wc$n3))
  Vjc <- (0.25 * Vjc) - Ejc^2
  nrns <- function(x, op = "*") {
    k <- length(x)
    res <- numeric(((k^2) - k)/2)
    ii <- 1
    for (i in 2:k) {
      for (j in 1:(i - 1)) {
        if (is.character(op) && op == "*") {
          res[ii] <- x[i] * x[j]
        }
        else if (is.character(op) && op == "+") {
          res[ii] <- x[i] + x[j]
        }
        ii <- ii + 1
      }
    }
    res
  }
  Exp <- (wc$S0 * (nrns(ntab, op = "*")))/(N * wc$n1)
  Var <- (2 * wc$S1 * nrns(ntab, op = "*"))/(N * wc$n1)
  Var <- Var + (((wc$S2 - 2 * wc$S1) * nrns(ntab, op = "*") *
                   (nrns(ntab, op = "+") - 2))/(N * wc$n1 * wc$n2))
  Var <- Var + ((4 * (S02 + wc$S1 - wc$S2) * nrns((ntab *
                                                     (ntab - 1)), op = "*"))/(N * wc$n1 * wc$n2 * wc$n3))
  Var <- (0.25 * Var) - Exp^2


  JtotExp <- sum(Exp)
  Jvar <- ((wc$S2/(N * wc$n1)) - ((4 * (S02 + wc$S1 - wc$S2) *
                                     wc$n1)/(N * wc$n1 * wc$n2 * wc$n3))) * sum(nrns(ntab,
                                                                                     op = "*"))
  Jvar <- Jvar + 4 * (((wc$S1 - wc$S2)/(N * wc$n1 * wc$n2 *
                                          wc$n3)) + ((2 * S02 * (2 * n - 3))/((N * wc$n1) * (N *
                                                                                               wc$n1 * wc$n2 * wc$n3)))) * sum(nrns(ntab^2, op = "*"))
  if (k > 2) {
    ntnsnr <- as.numeric(0)
    for (r in 1:(k - 2)) {
      for (s in (r + 1):(k - 1)) {
        for (t in (s + 1):(k)) {
          ntnsnr <- ntnsnr + ntab[r] * ntab[s] * ntab[t]
        }
      }
    }
    Jvar <- Jvar + (((2 * wc$S1 - 5 * wc$S2)/(N * wc$n1 *
                                                wc$n2)) + ((12 * (S02 + wc$S1 - wc$S2))/(N * wc$n1 *
                                                                                           wc$n2 * wc$n3)) + ((8 * S02)/((N * wc$n1 * wc$n2) *
                                                                                                                           wc$n1))) * ntnsnr
  }
  if (k > 3) {
    nuntnsnr <- as.numeric(0)
    for (r in 1:(k - 3)) {
      for (s in (r + 1):(k - 2)) {
        for (t in (s + 1):(k - 1)) {
          for (u in (t + 1):(k)) {
            nuntnsnr <- nuntnsnr + ntab[r] * ntab[s] *
              ntab[t] * ntab[u]
          }
        }
      }
    }
    Jvar <- Jvar - 8 * (((wc$S1 - wc$S2)/(N * wc$n1 * wc$n2 *
                                            wc$n3)) + ((2 * S02 * (2 * N - 3))/((N * wc$n1) *
                                                                                  (N * wc$n1 * wc$n2 * wc$n3)))) * nuntnsnr
  }
  Jvar <- (0.25 * Jvar)

  ##POTENTIAL PERMUTATION SECTION
  y <- factor(paste(ifx[sn[, 1]], ifx[sn[, 2]], sep = ":"),
              levels = as.vector(outer(1:k, 1:k, FUN = function(X,
                                                                Y) paste(X, Y, sep = ":"))))
  res <- matrix(tapply(sn[, 3], y, sum), ncol = k)/2
  res[is.na(res)] <- 0
  rownames(res) <- colnames(res) <- levels(fx)
  ldiag <- numeric(((k^2) - k)/2)
  diffcolnames <- character(((k^2) - k)/2)
  ii <- 1
  for (i in 2:k) {
    for (j in 1:(i - 1)) {
      ldiag[ii] <- res[i, j] + res[j, i]
      diffcolnames[ii] <- paste(levels(fx)[i], levels(fx)[j],
                                sep = ":")
      ii <- ii + 1
    }
  }
  Jtot <- sum(ldiag)
  statistic <- (c(diag(res), ldiag, Jtot) - c(Ejc, Exp, JtotExp))/sqrt(c(Vjc,Var, Jvar))
  lres <- cbind(c(diag(res), ldiag, Jtot), c(Ejc, Exp, JtotExp),
                c(Vjc, Var, Jvar), statistic)
  colnames(lres) <- c("Joincount", "Expected", "Variance","z-value")
  rownames(lres) <- c(paste(levels(fx), ":", levels(fx), sep = ""),
                      diffcolnames, "Jtot")
  return(as.data.frame(lres))
}


#' @title JC colocalization test and plot using Seurat/Kandinsky object
#' @name jc_coloc
#' @description
#' Call multi_jc() function to perform colocalization test, and plot results in the form of a heatmap
#'
#' @param data a Seurat object containing Kandinsky data slot
#' @param label name of the metadata column variable containing cell type annotation to be used for join count test
#' @param max.cap numeric, quantile cap to apply to colocalisation coefficients for colour scale creation. Default is 0.9, will set the coefficient value limit to the 90% quantile of all coefficients
#' @param return.mat boolean, whether returning only colocalization plot (FALSE) or also a dataframe containing colocalization test results
#' @param wc.type input type to use to extract summary weight constants used to calculate join count statistics. Must be one of the following: "listw", "mat"
#' @return a colocalization plot (when return.mat is FALSE) or a list containing colocalization plot and results in a tabular format
#' @seealso [squidpy_coloc()], [cellcharter_coloc()], [hoodscanR_coloc()], [giotto_coloc()] for alternative colocalisation methods
#' @export
jc_coloc = function(data=NULL,label=NULL,max.cap=0.9,return.mat=F,
                    wc.type=c('mat','listw')){
  if(!is.factor(data@meta.data[[label]])){
    data@meta.data[[label]] = as.factor(data@meta.data[[label]])
  }
  if(any(methods::slotNames(KanData(data))=='listw')){
    jc_mat = multi_jc(data@meta.data[[label]],
                      KanData(data,'listw'),wc.type=wc.type)
  }else{
    jc_mat = multi_jc(data@meta.data[[label]],
                      KanData(data,'nb'),wc.type=wc.type)
  }
  #jc_mat = as.data.frame(jc_mat)
  jc_mat$V1 = sapply(rownames(jc_mat),function(x){strsplit(x,split="\\:")[[1]][[1]]})
  jc_mat = jc_mat[jc_mat$V1 != 'Jtot',]
  jc_mat$V2 = sapply(rownames(jc_mat),function(x){strsplit(x,split="\\:")[[1]][[2]]})

  if(is.factor(data@meta.data[[label]])){
    levels= levels(data@meta.data[[label]])
    jc_mat[["V2"]] = factor(jc_mat[["V2"]],levels=levels)
    jc_mat[["V1"]] = factor(jc_mat[["V1"]],levels=levels)
  }
  jc_mat[["oe_log2ratio"]] = log2(jc_mat$Joincount/jc_mat$Expected)
  jc_mat[!is.finite(jc_mat$oe_log2ratio),]$oe_log2ratio = NA
    #g=ggplot(jc_mat,aes(y=.data[["V2"]],x=.data[["V1"]],fill=.data[["z-value"]]))+
    #  theme_minimal()+geom_tile(color='gray20')+
    #  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,
    #                       limits=c(quantile(jc_mat[['z-value']],0.9)*-1,
    #                                quantile(jc_mat[['z-value']],0.9)),
    #                       oob = scales::squish)+
    #  theme(axis.text.x=element_text(angle=45,hjust=1))+
    #  theme(panel.grid.major.x = element_blank(),axis.text = element_text(color='black'))+
    #  labs(x='Cell Type',y='Cell Type',fill='Coefficient')
  breaks = tryCatch(seq(round(min(abs(jc_mat$oe_log2ratio),na.rm=T)),round(max(abs(jc_mat$oe_log2ratio),na.rm=T)),1),error=function(e){e})
  if(length(breaks) >=5){ breaks = breaks[c(1,round(median(length(breaks))),length(breaks))] }
  if(inherits(breaks, "error")){breaks=seq(1,4,1)}
  
  g=ggplot(jc_mat, aes(y = .data[["V2"]], x = .data[["V1"]],fill = .data[["z-value"]],size=abs(.data[["oe_log2ratio"]]))) +
    theme_minimal() + geom_point(shape=21,color = "black") +
    scale_size_continuous(range=c(0,9),
                          breaks=breaks,
                          labels = 2^breaks)+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         limits = c(quantile(jc_mat[["z-value"]],max.cap,na.rm=T) * -1,
                                    quantile(jc_mat[["z-value"]], max.cap,na.rm=T)),
                         oob = scales::squish) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    theme(panel.grid.major.x = element_blank(),
          axis.text = element_text(color = "black")) +
    labs(x = '',y = '', fill = "z-value",size='|Obs/Exp| ratio')+
    theme(axis.title.x  = element_blank(),
          axis.title.y = element_blank())
  if(return.mat ==T){
    coloc_mat=jc_mat
    colnames(coloc_mat)[colnames(coloc_mat) == 'V1'] = 'Type1'
    colnames(coloc_mat)[colnames(coloc_mat) == 'V2'] = 'Type2'
    coloc_mat$oe_ratio = 2^(abs(coloc_mat$oe_log2ratio))
  }
  if(return.mat==F){
    g
  }else{
    return(list(coloc_mat=coloc_mat,plot=g))
  }
}



#nnMat = KanData(seurat,'nnMat')$nnMat
#nnMat[,'tot_nn'] = NULL
#label_mat = reshape2::dcast(nnMat %>% tibble::rownames_to_column(var='ID'),
#                formula=ID~cluster,value.var = 'cluster',
#                fun.aggregate = length) %>%
#  tibble::column_to_rownames(var='ID') %>%
#  as.matrix() %>%
#  as(.,'CsparseMatrix')
#label_var = nnMat[,'cluster']
#nnMat = nnMat %>% dplyr::select_if(is.numeric)
#nnMat = as(as.matrix(nnMat),'CsparseMatrix') + label_mat[rownames(nnMat),]
#make binary
#nnMat[nnMat >1] = 1



#' @title estimate weight summary constansts from spatial weights matrix
#' @name spweights_mat
#' @description
#' Faster extraction of neighbour network summary weights (taken from spweights.constants function in spdep )
#'
#' @details
#' For more information see spweights.constants() function documentation
#'
#' @param listw a weighted network in list (listw) or matrix (CsparseMatrix) format
#' @returns a list obejct containing summary constants
#'
#' @export
spweights_mat = function(listw){
  if(inherits(listw, "listw")){listw=as(listw,'CsparseMatrix')}
  n = ncol(listw)
  n1 <- n - 1
  n2 <- n - 2
  n3 <- n - 3
  nn <- n^2
  c1 = rowSums(listw)
  S0 <- sum(c1)
  S1 <- sum((listw * listw) + (listw * Matrix::t(listw)))
  S2 <- sum((Matrix::rowSums(listw) + Matrix::colSums(listw))^2)

  list(n = n, n1 = n1, n2 = n2, n3 = n3, nn = nn, S0 = S0, S1 = S1, S2 = S2)
}




#' @title Squidpy colocalisation analysis
#' @name squidpy_coloc
#' @description
#' Squidpy is called thrugh reticulate and interacts with Seurat/Kandinsky object after anndata conversion with anndataR
#' 
#' @details
#' Squidpy must be installed within a python environment together with its dependencies
#' 
#' @param data a Seurat object containing Kandinsky data slot
#' @param label name of the metadata column variable containing cell type annotation to be used for join count test
#' @param seed numeric, random seed
#' @param nperms numeric, number of random permutations to perform for colocalisation test
#' @param python_path path of python environment with MESA and its dependencies installed
#' @param max.cap numeric, quantile cap to apply to colocalisation coefficients for colour scale creation. Default is 0.9, will set the coefficient value limit to the 90% quantile of all coefficients
#' @param return.mat boolean, whether returning only colocalization plot (FALSE) or also a dataframe containing colocalization test results
#' @return a colocalization plot (when return.mat is FALSE) or a list containing colocalization plot and results in a tabular format
#' @export
squidpy_coloc = function(data=NULL,
			 label = NULL,
			 seed=347548,nperms=100,python_path=NULL,
			 max.cap=0.9,
			 return.mat=T
			 ){
	message('Initializing Squidpy with reticulate')
	python_path = python_path %||% Sys.getenv('RETICULATE_PYTHON')
	if(is.null(python_path) & Sys.getenv('RETICULATE_PYTHON') == "") stop('Please specificy a path for the python environment containing Squidpy with the parameter "python_path".\nExample: /path/to/my/version/of/python')
	if(Sys.getenv('RETICULATE_PYTHON') == ""){Sys.setenv(RETICULATE_PYTHON = python_path)}
	if(!requireNamespace("reticulate",quietly=T)) stop('You need to install reticulate to run Squidpy.\n Run install.packages("reticulate")')
	if(!requireNamespace("anndataR",quietly=T)) stop('You need to install anndataR to run Squidpy.\n Run BiocManager::install("anndataR") or conda install bioconda::anndataR if you are working within a conda environment')
	message('Loading Squidpy...')
	np = reticulate::import('numpy')
	sq = reticulate::import('squidpy')
	message('Converting Seurat object to anndata...')
	coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
	colnames(coords) = c('X','Y')
	adata = reticulate::r_to_py(anndataR::as_AnnData(data,assay_name=DefaultAssay(data),output_class='ReticulateAnnData'))
	adata$obsm['spatial']= coords
	adata$obsp['spatial_connectivities'] = reticulate::r_to_py(as(KanData(data,'nb'),'CsparseMatrix'))
	adata$obsp['spatial_distances'] = adata$obsp['spatial_connectivities']
	sq$gr$nhood_enrichment(adata, cluster_key=label, n_perms=as.integer(nperms))
	name = paste0(label,'_nhood_enrichment')
	zmat_squidpy = reticulate::py_to_r(adata$uns[name]['zscore'])
	dimnames(zmat_squidpy) = list(levels(data@meta.data[,label]),levels(data@meta.data[,label]))
	res = reshape2::melt(zmat_squidpy)
	colnames(res)=c('central_ct','nb_ct','coef')
	g=ggplot(res, aes(y = .data[["central_ct"]], x = .data[["nb_ct"]],fill = .data[["coef"]])) +
	theme_minimal() + geom_point(shape=21,color = "black") +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red",
			     midpoint = 0,
			     limits = c(quantile(res$coef,max.cap) * -1,
					quantile(res$coef, max.cap)),
			     oob = scales::squish) +
	theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
	theme(panel.grid.major.x = element_blank(),
	      axis.text = element_text(color = "black")) +
	labs(x = '',y = '', fill = "coefficient")+
	theme(axis.title.x  = element_blank(),
	      axis.title.y = element_blank())
	if(return.mat==F){
		g
	}else{
		return(list(coloc_mat=res,plot=g))
	}
}

#' @title CellCharter colocalisation analysis
#' @name cellcharter_coloc
#' @description
#' CellCharter is called thrugh reticulate and interacts with Seurat/Kandinsky object after anndata conversion with anndataR
#' 
#' @details
#' CellCharter must be installed within a python environment together with its dependencies (squidpy)
#' 
#' @param data a Seurat object containing Kandinsky data slot
#' @param label name of the metadata column variable containing cell type annotation to be used for join count test
#' @param seed numeric, random seed
#' @param log_fold_change boolean, whether to return default CellCharter colocalisation scores (observed - expected) or their log-fold change version (log2(observed/expected))
#' @param python_path path of python environment with MESA and its dependencies installed
#' @param max.cap numeric, quantile cap to apply to colocalisation coefficients for colour scale creation. Default is 0.9, will set the coefficient value limit to the 90% quantile of all coefficients
#' @param return.mat boolean, whether returning only colocalization plot (FALSE) or also a dataframe containing colocalization test results
#' @return a colocalization plot (when return.mat is FALSE) or a list containing colocalization plot and results in a tabular format
#' @export
cellcharter_coloc = function(data=NULL,
			     label = NULL,
			     seed=347548,
			     log_fold_change=F,
			     python_path=NULL,
			     max.cap=0.9,
			     return.mat=T){
	message('Initializing CellCharter with reticulate')
	python_path = python_path %||% Sys.getenv('RETICULATE_PYTHON')
	if(is.null(python_path) & Sys.getenv('RETICULATE_PYTHON') == "") stop('Please specificy a path for the python environment containing CellCharter with the parameter "python_path".\nExample: /path/to/my/version/of/python')
	if(Sys.getenv('RETICULATE_PYTHON') == ""){Sys.setenv(RETICULATE_PYTHON = python_path)}
	if(!requireNamespace("reticulate",quietly=T)) stop('You need to install reticulate to run CellCharter.\n Run install.packages("reticulate")')
	if(!requireNamespace("anndataR",quietly=T)) stop('You need to install anndataR to run CellCharter.\n Run BiocManager::install("anndataR") or conda install bioconda::anndataR if you are working within a conda environment')
       message('Loading CellCharter and its dependencies...')
	np = reticulate::import('numpy')
	sq = reticulate::import('squidpy')
	cc = reticulate::import('cellcharter')
	message('Converting Seurat object to anndata...')
	coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
	colnames(coords) = c('X','Y')
	adata = reticulate::r_to_py(anndataR::as_AnnData(data,assay_name=DefaultAssay(data),output_class='ReticulateAnnData'))
	adata$obsm['spatial']= coords
	adata$obsp['spatial_connectivities'] = reticulate::r_to_py(as(KanData(data,'nb'),'CsparseMatrix'))
	adata$obsp['spatial_distances'] = adata$obsp['spatial_connectivities']
	zmat_cellcharter =cc$gr$nhood_enrichment(adata,cluster_key=label,copy=T,only_inter=F,log_fold_change=log_fold_change,
						 observed_expected=F,pvalues=F)$enrichment
	res = reshape2::melt(as.matrix(zmat_cellcharter))
	colnames(res)=c('central_ct','nb_ct','coef')
	g=ggplot(res, aes(y = .data[["central_ct"]], x = .data[["nb_ct"]],fill = .data[["coef"]])) +
	theme_minimal() + geom_point(shape=21,color = "black") +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red",
			     midpoint = 0,
			     limits = c(quantile(res$coef,max.cap) * -1,
					quantile(res$coef, max.cap)),
			     oob = scales::squish) +
	theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
	theme(panel.grid.major.x = element_blank(),
	      axis.text = element_text(color = "black")) +
	labs(x = '',y = '', fill = "coefficient")+
	theme(axis.title.x  = element_blank(),
	      axis.title.y = element_blank())
	if(return.mat==F){
		g
	}else{
		return(list(coloc_mat=res,plot=g))
	}
}

	
	
#' @title hoodscanR colocalisation
#' @name hoodscanR_coloc
#' @description
#' hoodscanR colocalisation performed by calling the R package hoodscanR through Kandinsky
#' 
#' @details
#' hoodscanR must be installed within the same environmnent of Kandinsky..
#' 
#' @param data Seurat object containing Kandinsky data 
#' @param k numeric, number of nearest neighbours to consider for hoodscanR clustering
#' @param label character string specifying the variable name to be used to defne cell annotation groups
#' @param sample_key character string specifying a variable stored in the Seurat object to use as sample/batch annotation. If not NULL, neighbour networks will be defined separately for each sample/batch. Default is NULL
#' @param max.cap numeric, quantile cap to apply to colocalisation coefficients for colour scale creation. Default is 1, will set the coefficient value limit to the highest value of all coefficients
#' @param return.mat boolean, whether returning only colocalization plot (FALSE) or also a dataframe containing colocalization test results
#' @return a colocalization plot (when return.mat is FALSE) or a list containing colocalization plot and results in a tabular format
#' @export
hoodscanR_coloc = function(data=NULL,k=100,label=NULL,sample_key=NULL,max.cap=1,return.mat=T){
	if((!requireNamespace("hoodscanR",quietly = T))){
		stop('Please install hoodscanR first by running BiocManager::install("hoodscanR")')
	}
	if(is.null(label)){
		stop('Please specify "label" argument.')
	}
	if(!is.null(sample_key)){
		samples = unique(data@meta.data[,sample_key])
	}
	message('Preparing data for hoodscanR...')
	spe = SpatialExperiment::SpatialExperiment(assays=list(counts=LayerData(data,layer='counts')),
						   rowData=data.frame(row.names=rownames(data)),
						   colData=data@meta.data,metadata=list(),reducedDims=list(),
						   altExps=list(),sample_id = 'my_spe',
						   spatialCoords=sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
						   )
      message('Running hoodscanR...')
      if(!is.null(sample_key)){
	      nns = list()
	      for(s in samples){
		      nns[[s]] = hoodscanR::findNearCells(spe[,spe@colData[,sample_key] ==s], k = k,anno_col=label)
	      }
	      cells = lapply(nns,function(x){x[[1]]})
	      dists = lapply(nns,function(x){x[[2]]})
	      fnc = list(cells=purrr::reduce(cells,rbind),distance=purrr::reduce(dists,rbind))
      }else{
	      fnc <- hoodscanR::findNearCells(spe, k = k,anno_col=label)
      }
      pm <-  hoodscanR::scanHoods(fnc$distance)
      hoods <-  hoodscanR::mergeByGroup(pm, fnc$cells)
      zmat_hoodscanr = cor(hoods)
      res = reshape2::melt(as.matrix(zmat_hoodscanr))
      colnames(res)=c('central_ct','nb_ct','coef')
      g=ggplot(res, aes(y = .data[["central_ct"]], x = .data[["nb_ct"]],fill = .data[["coef"]])) +
      theme_minimal() + geom_point(shape=21,color = "black") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
			   midpoint = 0,
			   limits = c(quantile(res$coef,max.cap) * -1,
				      quantile(res$coef, max.cap)),
			   oob = scales::squish) +
	theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
	theme(panel.grid.major.x = element_blank(),
	      axis.text = element_text(color = "black")) +
	labs(x = '',y = '', fill = "coefficient")+
	theme(axis.title.x  = element_blank(),
	      axis.title.y = element_blank())
	if(return.mat==F){
		g
	}else{
		return(list(coloc_mat=res,plot=g))
	}
}


#' @title Giotto colocalisation analysis
#' @name giotto_coloc
#' @description
#' Giotto colocalisation performed by calling the R package Giotto through Kandinsky
#' 
#' @details
#' Giotto must be installed within the same environmnent of Kandinsky
#' 
#' @param data a Seurat object containing Kandinsky data slot
#' @param label name of the metadata column variable containing cell type annotation to be used for join count test
#' @param nperms numeric, number of random permutations to perform for colocalisation test
#' @param python_path path of python environment with MESA and its dependencies installed
#' @param max.cap numeric, quantile cap to apply to colocalisation coefficients for colour scale creation. Default is 0.9, will set the coefficient value limit to the 90% quantile of all coefficients
#' @param return.mat boolean, whether returning only colocalization plot (FALSE) or also a dataframe containing colocalization test results
#' @return a colocalization plot (when return.mat is FALSE) or a list containing colocalization plot and results in a tabular format
#' @export
giotto_coloc = function(data=NULL,label=NULL,nperms=1000,python_path=NULL,max.cap=0.9,return.mat=T){
	if((!requireNamespace("Giotto", quietly = TRUE))){
		stop('Please install Giotto first by running pak::pkg_install("giotto-suite/Giotto")')
	}
Giotto::createGiottoInstructions(python_path = python_path)
genv_exists <- Giotto::checkGiottoEnvironment()
if(!genv_exists){
    # The following command need only be run once to install the Giotto environment
	warning('python env not found. Giotto will create/search its default environment')
	Giotto::installGiottoEnvironment()
}
message('Preparing data for Giotto...')
data$cell_ID = colnames(data)
  
#Create Giotto object
  spat_coord = sf::st_coordinates(sf::st_centroid(KanData(data,'sf')$geometry))
  spat_coord = as.data.frame(spat_coord)
  rownames(spat_coord) = colnames(data)
  spat_coord$cell_ID <- rownames(spat_coord)
  colnames(spat_coord) <- c("sdimx", "sdimy", "cell_ID")
  adj = as(KanData(data,'nb'),'CsparseMatrix')
  s <- Matrix::summary(adj)
  
  networkDT <- data.table::data.table(
    from = rownames(adj)[s[,1]],
    to   = colnames(adj)[s[,2]],
    weight = s[,3]
  )
  giotto <- Giotto::createGiottoObject(
    expression = list(normalized=LayerData(data,'counts')),
    spatial_locs = spat_coord,cell_metadata=data.table::as.data.table(data@meta.data)
  )   
		  
  giotto <- Giotto::createSpatialNetwork(giotto,
					 method = "Delaunay", 
					 name = "spatial_network"
					 )
  snet <- Giotto::getSpatialNetwork(
				    giotto,
				    name = "spatial_network"
				    )
		    
  slot(snet, "networkDT") <- networkDT	    
  gobject <- Giotto::setSpatialNetwork(
				       giotto,
				       x = snet,
				       name = "spatial_network"
				       ) 
  giotto = Giotto::identifyTMAcores(giotto,spatial_network_name='spatial_network') # useful to split networks from independent samples
  res = Giotto::cellProximityEnrichment(
					giotto,
					cluster_column        = label,
					spatial_network_name  = "spatial_network",
					number_of_simulations = nperms
					)
  res = as.data.frame(res)
  res$central_ct = sapply(res$unified_int,function(x){strsplit(x,split='--')[[1]][[1]]})
  res$nb_ct = sapply(res$unified_int,function(x){strsplit(x,split='--')[[1]][[2]]})
  res$coef = res$PI_value
  res=res[,c('central_ct','nb_ct','coef','p_higher_orig','p_lower_orig')]
  g=ggplot(res, aes(y = .data[["central_ct"]], x = .data[["nb_ct"]],fill = .data[["coef"]])) +
    theme_minimal() + geom_point(shape=21,color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         limits = c(quantile(res$coef,max.cap) * -1,
                                    quantile(res$coef, max.cap)),
                         oob = scales::squish) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    theme(panel.grid.major.x = element_blank(),
          axis.text = element_text(color = "black")) +
    labs(x = '',y = '', fill = "coefficient")+
    theme(axis.title.x  = element_blank(),
          axis.title.y = element_blank())
  if(return.mat==F){
    g
  }else{
    return(list(coloc_mat=res,plot=g))
  }
}


