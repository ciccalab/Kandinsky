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
#' @param seurat a Seurat object containing Kandinsky data slot
#' @param label name of the metadata column variable containing cell type annotation to be used for join count test
#' @param return.mat boolean, whether returning only colocalization plot (FALSE) or also a dataframe containing colocalization test results
#' @param wc.type input type to use to extract summary weight constants used to calculate join count statistics. Must be one of the following: "listw", "mat"
#' @returns a colocalization plot (when return.mat is FALSE) or a list containing colocalization plot and results in a tabular format
#' @export
jc_coloc = function(seurat=NULL,label=NULL,return.mat=F,
                    wc.type=c('mat','listw')){
  if(!is.factor(seurat@meta.data[[label]])){
    seurat@meta.data[[label]] = as.factor(seurat@meta.data[[label]])
  }
  if(any(methods::slotNames(KanData(seurat))=='listw')){
    jc_mat = multi_jc(seurat@meta.data[[label]],
                      KanData(seurat,'listw'),wc.type=wc.type)
  }else{
    jc_mat = multi_jc(seurat@meta.data[[label]],
                      KanData(seurat,'nb'),wc.type=wc.type)
  }
  #jc_mat = as.data.frame(jc_mat)
  jc_mat$V1 = sapply(rownames(jc_mat),function(x){strsplit(x,split="\\:")[[1]][[1]]})
  jc_mat = jc_mat[jc_mat$V1 != 'Jtot',]
  jc_mat$V2 = sapply(rownames(jc_mat),function(x){strsplit(x,split="\\:")[[1]][[2]]})

  if(is.factor(seurat@meta.data[[label]])){
    levels= levels(seurat@meta.data[[label]])
    jc_mat[["V2"]] = factor(jc_mat[["V2"]],levels=levels)
    jc_mat[["V1"]] = factor(jc_mat[["V1"]],levels=levels)
  }
  jc_mat[["oe_log2ratio"]] = log2(jc_mat$Joincount/jc_mat$Expected)
    #g=ggplot(jc_mat,aes(y=.data[["V2"]],x=.data[["V1"]],fill=.data[["z-value"]]))+
    #  theme_minimal()+geom_tile(color='gray20')+
    #  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,
    #                       limits=c(quantile(jc_mat[['z-value']],0.9)*-1,
    #                                quantile(jc_mat[['z-value']],0.9)),
    #                       oob = scales::squish)+
    #  theme(axis.text.x=element_text(angle=45,hjust=1))+
    #  theme(panel.grid.major.x = element_blank(),axis.text = element_text(color='black'))+
    #  labs(x='Cell Type',y='Cell Type',fill='Coefficient')
  g=ggplot(jc_mat, aes(y = .data[["V2"]], x = .data[["V1"]],fill = .data[["z-value"]],size=abs(.data[["oe_log2ratio"]]))) +
    theme_minimal() + geom_point(shape=21,color = "black") +
    scale_size_continuous(range=c(0,9),
                          breaks=c(seq(round(min(abs(jc_mat$oe_log2ratio))),round(max(abs(jc_mat$oe_log2ratio))),1)),
                          labels = 2^(c(seq(round(min(abs(jc_mat$oe_log2ratio))),round(max(abs(jc_mat$oe_log2ratio))),1))))+
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         limits = c(quantile(jc_mat[["z-value"]],0.9) * -1,
                                    quantile(jc_mat[["z-value"]], 0.9)),
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




