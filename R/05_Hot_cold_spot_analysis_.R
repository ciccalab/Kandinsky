###Calculate Getis-Ord Gi statistic for hotspot analysis of a variable of choice.
####Classify spots/cells in 'Hot' and 'Cold' spatial clusters
#' @title Hotspot analysis using Getis-Ord Gi* statistics
#' @name hotspot_analysis
#' @description
#' Apply Getis-Ord Gi statistic to any numeric variable stored in Seurat object and calculate spatial local hot- and cold-spots for that variable
#' @details
#' The function will use the `localG_perm()` function from the `spdep` package to calculate, for a numeric variable of interest (e.g. expression level of a gene), the local Getis-Ord Gi* statistic for each spot/cell analysed.
#' The `hotspot` function from the same package wil be then use to define local hot-spots or cold-spots on the basis of the spot/cell Gi statistics. These results should be interpreted as follows:
#' - Spots / cells classified as hot-spots show significantly high values of the variable of interest, as well as their neighbouring spots/cells
#' - Spots / cells classified as cold-spots show significantly low values of the variable of interest, as well as their neighbouring spots/cells
#' For more information about Getis-Ord Gi* statistics implemented in `spdep`, see https://r-spatial.github.io/spdep/reference/localmoran.html
#'
#' @param seurat a Seurat object containing Kandinsky data slot(`KanData()`)
#' @param feature character string specifying the name of the feature to use to compute Getis-Ord Gi* statistics
#' @param layer character string indicating which Seurat layer will be considered to search for the variable of interest
#' @param sim number of Monte Carlo simulations to be run for estimating Local Gi* coefficients significance
#' @param padj.thresh numeric value indicating the significance threshold to be applied to the adjusted pvalues resulting from Monte Carlo simulations
#' @param lag integer value indicating the extent of cell/spot neighbours to be considered to calculate Getis-Ord Gi* statistics. `lag = 1` indicates that only 1st order neighbours will be considered, while `lag=2` indicates that all neighbours of each 1st order neighbour will be also considered for each spot/cell, and so on.
#' When the average number or 1st order neighbours is limited as in the case of Visium spots (that is, by default, only spots/bins composing the ring immediately surrounding each spot/bin), increasing the lag might help in reducing the variability of the final results. Default value is set to 1.
#' @param seed numeric, random seed for reproducibility. Default is set to 347548
#' @param prune boolean, whether ignoring any spot/cell annotated as hot or cold spot with no neighbouring hot or cold spots, respectively (i.e., singlets)
#' @returns updated Seurat object with two new columns stored in the meta.data table: `feature_name_GI` and `feature_name_clust`, reporting the local Gi* statistics and the association with a hot- or cold-spot for each spot/cell, respectively
#' @export
hotspot_analysis = function(seurat=NULL,feature=NULL,layer='data',sim=999,padj.thresh=0.05,lag=1,seed=347548,prune=F){
  seurat= update_kandinsky(seurat)
  if(length(KanData(seurat,'nb')) == 0){
    stop('Kandinsky data must contain a non-empty "nb" data slot')
  }
  subset = Seurat::FetchData(seurat,vars=feature,layer=layer,clean=F)
  if(lag > 1){
    seurat = nb_expand(seurat,maxorder=lag,cumul=T)
  }
  listw = KanData(seurat,'nb')
  set.seed(seed)
  gi = spdep::localG_perm(subset[,feature], listw, nsim=sim, zero.policy=T, spChk=NULL,
                          alternative = "two.sided", iseed=NULL,
                          no_repeat_in_row=FALSE)
  padj = attr(gi,'internals') %>% as.data.frame() %>% dplyr::select(.data[["Pr(z != E(Gi)) Sim"]])
  padj = p.adjust(padj[,1],'BH')
  gi = cbind(as.data.frame(gi),(padj))
  colnames(gi) = c('Gi','padj')
  gi = gi %>% dplyr::mutate(cluster = ifelse(.data[["Gi"]] > 0 & .data[["padj"]] < padj.thresh,'Hot',
                                      ifelse(.data[["Gi"]] <0 & .data[["padj"]] < padj.thresh,'Cold','NS')))
  
  # gi[is.na(gi$cluster),]$cluster='NS'
  if (any(is.na(gi$cluster))) {
      gi[is.na(gi$cluster), ]$cluster = "NS"
    }
  
  gi$cluster = factor(gi$cluster,levels=c('Cold','Hot','NS'))
  #hotspot=as.character(spdep::hotspot(gi,Prname='Pr(z != E(Gi)) Sim',cutoff=padj.thresh))
  #hotspot[is.na(hotspot)]='NS'
  #hotspot[hotspot == 'High']='Hot'
  #hotspot[hotspot == 'Low']='Cold'
  #hotspot = factor(hotspot,levels=c('Cold','Hot','NS'))
  colname1 = paste0(feature,"_Gi")
  colname2 = paste0(feature,'_clust')
  seurat@meta.data[,c(colname1,colname2)] = gi[,c('Gi','cluster')]
  if(prune ==T){
    hot = rownames(KanData(seurat,'sf')[rownames(seurat@meta.data[seurat@meta.data[[colname2]]=='Hot',]),])
    if(length(hot) >0){
      hot_nb = suppressWarnings(subset(listw$neighbours,attr(listw$neighbours,'region.id') %in% hot))
      if(sum(spdep::card(hot_nb$neighbours) == 0)>0){
        hot_rm = attributes(hot_nb)$region.id[spdep::card(hot_nb$neighbours) == 0]
        seurat@meta.data[hot_rm,colname2] = 'NS'
      }
    }
    cold = rownames(KanData(seurat,'sf')[rownames(seurat@meta.data[seurat@meta.data[[colname2]]=='Cold',]),])
    if(length(cold) >0){
      cold_nb = suppressWarnings(subset(listw$neighbours,attr(listw$neighbours,'region.id') %in% cold))
      if(sum(spdep::card(cold_nb$neighbours) == 0)>0){
        cold_rm = attributes(cold_nb)$region.id[spdep::card(cold_nb$neighbours) == 0]
        seurat@meta.data[cold_rm,colname2] = 'NS'
      }
    }
  }
  return(seurat)
}

#' @title Count and describe Getis-Ord Gi* hotspot sites
#' @name hotspot_count
#' @description
#' Summarise number and min/max/median size of hot- and cold-spots for a variable of choice.
#' @details
#' This function will look for the presence of already existing hotspot annotation in the meta.data of Seurat object for the variable(s) of interest.
#' If no annotations exist, the `hotspot_analysis()` function will be automatically called.
#'
#' @param data a Seurat object containing Kandinsky data (`KanData()`)
#' @param feature character string specifying the name of the feature(s) used to compute Getis-Ord Gi* statistics
#' @param ... arguments passed to hotspot_analysis function
#' @returns data.frame summarising hotspot and coldspot size and number for each variable of interest
#' @export
hotspot_count = function(data=NULL,feature=NULL,...){
  all_feat = lapply(feature,function(feature){
    col_name = paste0(feature,"_clust")
    if(length(intersect(col_name,colnames(data@meta.data))) == 0){
      data = hotspot_analysis(data,feature=feature,...)
    }
    hot_id = rownames(data@meta.data[data@meta.data[[col_name]]=='Hot',])
    cold_id = rownames(data@meta.data[data@meta.data[[col_name]]=='Cold',])
    if(length(hot_id)>0){
      hot = spdep::mat2listw(as(KanData(data,'nb'),'CsparseMatrix')[hot_id,hot_id],style='B',zero.policy = T)$neighbours
    }
    if(length(cold_id)>0){
      cold = spdep::mat2listw(as(KanData(data,'nb'),'CsparseMatrix')[cold_id,cold_id],style='B',zero.policy = T)$neighbours
    }
    suppressWarnings({
      if(length(hot_id)>0){
        n_hot = spdep::n.comp.nb(hot)$nc
        quantiles_hot = quantile(as.numeric(table(spdep::n.comp.nb(hot)$comp.id)))
        names(quantiles_hot) = c('hot_minsize','hot_q25size','hot_medsize','hot_q75size','hot_maxsize')
        quantiles_hot = t(as.data.frame(quantiles_hot))
      }else{
        n_hot = 0
        quantiles_hot = rep(0,5)
        names(quantiles_hot) = c('hot_minsize','hot_q25size','hot_medsize','hot_q75size','hot_maxsize')
        quantiles_hot = t(as.data.frame(quantiles_hot))
      }
      if(length(cold_id)>0){
        n_cold = spdep::n.comp.nb(cold)$nc
        quantiles_cold = quantile(as.numeric(table(spdep::n.comp.nb(cold)$comp.id)))
        names(quantiles_cold) = c('cold_minsize','cold_q25size','cold_medsize','cold_q75size','cold_maxsize')
        quantiles_cold = t(as.data.frame(quantiles_cold))
      }else{
        n_cold = 0
        quantiles_cold = rep(0,5)
        names(quantiles_cold) = c('cold_minsize','cold_q25size','cold_medsize','cold_q75size','cold_maxsize')
        quantiles_cold = t(as.data.frame(quantiles_cold))
      }
    })
    all_quantiles = cbind(quantiles_hot,quantiles_cold)
    rownames(all_quantiles) = NULL
    return(cbind(data.frame(feature=feature,n_hotspot=n_hot,n_coldspot=n_cold),all_quantiles))
  })
  all_feat= purrr::reduce(all_feat,rbind)
  return(all_feat)
}


#' @title define overlap between hot/cold spots of two variables
#' @name hotspot_overlap
#' @description
#' Given two numeric variables, this function will define hot and cold spots for each variable,
#' and it will annotate cells/spots according to the presence/absence of overlap between the two hot/cold spot annotations.
#' If a hot/cold spot annotation for one or both variables already exist,hotspot analysis will not be rerun for that variable(s)
#'
#' @param data a Seurat object containing Kandinsky data (`KanData()`)
#' @param feat.1 character string specifying the name of the first feature to consider
#' @param feat.2 character string specifying the name of the second feature to consider
#' @param ... arguments passed to hotspot_analysis function
#' @returns updated Seurat object with a new annotation column reporting the type of overlap between hot/cold spot annotations:
#' - 'Hot-Hot': a cell/spot is annotated as 'hot' for both the considered variables
#' - 'Cold-Cold': a cell/spot is annotated as 'cold' for both the considered variables
#' - 'Cold-Hot': a cell/spot is annotated as 'cold' for one variable and 'hot' for the other
#' - 'Hot': a cell/spot is annotated as 'hot' only for one of the two variable (the other is annotated as non-significant 'NS')
#' - 'Cold': a cell/spot is annotated as 'cold' only for one of the two variable (the other is annotated as non-significant 'NS')
#' - 'NS': a cell/spot is not annotated as either 'hot' or 'cold' for any of the two variables
#' @export
hotspot_overlap = function(seurat=NULL,feat.1=NULL,feat.2=NULL,...){
  col_name = paste0(c(feat.1,feat.2),"_clust")
  missing = setdiff(col_name,colnames(seurat@meta.data))
  if(length(missing)>0){
    for(m in missing){
      seurat = hotspot_analysis(seurat,feature=gsub('_clust','',m),...)
    }
  }
  hotspot_overlap = seurat@meta.data[[col_name[1]]] == seurat@meta.data[[col_name[2]]]
  seurat@meta.data[[paste0(feat.1,'_',feat.2,'_clust')]] = ifelse(hotspot_overlap == T & seurat@meta.data[[col_name[1]]] == 'Hot','Hot-Hot',
                               ifelse(hotspot_overlap == T & seurat@meta.data[[col_name[1]]] == 'Cold','Cold-Cold',
                                      ifelse(hotspot_overlap==F & seurat@meta.data[[col_name[1]]] != 'NS' & seurat@meta.data[[col_name[2]]] != 'NS','Cold-Hot',
                                             ifelse(hotspot_overlap==F & (seurat@meta.data[[col_name[1]]] == 'Hot' | seurat@meta.data[[col_name[2]]] == 'Hot'),'Hot',
                                                    ifelse(hotspot_overlap==F & (seurat@meta.data[[col_name[1]]] == 'Cold' | seurat@meta.data[[col_name[2]]] == 'Cold'),'Cold','NS')
                                                    )
                                             )
                                      )
                               )
  seurat@meta.data[[paste0(feat.1,'_',feat.2,'_clust')]] = factor(seurat@meta.data[[paste0(feat.1,'_',feat.2,'_clust')]],levels = c('Cold-Cold','Cold','Cold-Hot','Hot','Hot-Hot','NS'))
  return(seurat)

}

