#' @importFrom ggplot2 ggplot theme_void geom_point geom_sf geom_text geom_tile geom_density annotate theme scale_fill_manual scale_fill_gradientn scale_fill_gradient2 scale_y_reverse guides guide_axis labs aes element_line element_rect element_text element_blank theme_classic theme_minimal stat_density2d_filled scale_color_manual after_stat guide_colorbar guide_legend vars facet_wrap ggtitle geom_bar scale_color_gradientn coord_fixed scale_size_continuous
NULL

#color palettes for categorical data


#' C24 discrete color palette
#'
#' 24-color palette to apply to categorical data
#' @format ## `c24`
#' @export
c24 <- c(
  'dodgerblue2', '#E31A1C', # red
  'green4',
  '#6A3D9A', # purple
  '#FF7F00', # orange
  'gold1',
  'skyblue2', '#FB9A99', # lt pink
  'palegreen2',
  '#CAB2D6', # lt purple
  '#FDBF6F', # lt orange
  'gray70', 'khaki2',
  'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4',
  'darkturquoise', 'green1', 'yellow4', 'yellow3',
  'darkorange4', 'brown'
)

#' C12 discrete color palette
#'
#' 12-color palette to apply to categorical data
#' @format ## `c12`
#' @export
c12 <- c('#88CCEE','#CC6677','#DDCC77','#117733',
         '#332288','#AA4499','#44AA99','#999933',
         '#882255','#661100','#6699CC','#888888')


#' @title Plot Kandinsky images
#' @name SlidePlot
#' @description
#' Visualize raster images stored in Kandinsky data
#'
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param max.px numeric value indicating the number of image pixels to plot. Higher values will give a better resolution and a  heavier image size
#' @param img_type type of raster image stored in Kandinsky data. Must be one of the following: `img`, `mask`
#' @param as.ggplot, boolean, whether returning a plot created with `tidyterra::geom_spatraster_rgb()` function (TRUE) or using `terra::plotRGB()` function
#' @returns Kandinsky image plot, as a base R plot (as.ggplot=F) or as a ggplot object (as.ggplot=T)
#' @export
SlidePlot = function(seurat = NULL,max.px=10^6,img_type=c('img','mask'),as.ggplot=T){
  if(is.null(KanData(seurat))){
    stop('Seurat object does not contain any Kandinsky data. You need to call kandinsky_init() before running this function')
  }
  if(as.ggplot==T){
    ggplot()+theme_void()+tidyterra::geom_spatraster_rgb(data=KanData(seurat,img_type),maxcell = max.px)
  }else{
    terra::plotRGB(KanData(seurat,img_type),maxcell=max.px)
  }
}


#' @title visualize Kandinsky and Seurat spatial data
#' @name KanPlot
#' @description
#' Plot function for Seurat/Kandinsky data, as an alternative to Seurat SpatialDim/SpatialFeature plots
#'
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param feature character string specifying the discrete or continuous variable to plot
#' @param fovs vector of fov identifiers selected for plotting. If NULL (default), all fovs will be considered for plotting
#' @param assay the name of Seurat assay to use to search for the feature to plot
#' @param layer the name of Seurat layer to use to search for the feature to plot
#' @param palette_cont character string or vector specifying the color palette to use for continuous variables. It can be a vector of color names, or the name of a `viridis` palette
#' @param palette_discrete character string or vector specifying the color palette to use for discrete variables.
#' @param alpha numeric value between 0 and 1 specifying the level of transparency to apply to polygon/spots/bins. 0 means completely transparent, 1 mean means completely opaque. If not specified, default value is set to 1.
#' @param show_border boolean, whether plotting (TRUE) or not (FALSE) polygon borders.
#' @param border.size numeric, polygon border thickness passed to `lwd` parameter of geom_sf() function when `show_border` is set to TRUE. If not specified, default value is set to 0.01
#' @param border.col color of polygon borders when `show_border` is set to TRUE.
#' @param max.cap numeric value between 0 and 1 specifying the higher percentile cutoff value for continuous variables. Default value is set to 0.99
#' @param bg_img boolean, whether plotting (TRUE) or not (FALSE) any image associated with Kandinsky data, like Visium H&E image,in the plot background.
#' @param bg_type image type to be plotted in the background when `bg_img = TRUE`. Must be one of the following: `img`, `mask`.
#' @param max.px numeric value indicating the number of pixels to plot from the background image, when `bg_img = TRUE`. higher values will give a better resolution and a  heavier image size
#' @param theme color style of the plot background and text. Must one of the following: `dark`, `light`.
#' @param pt.shape numeric, shape to apply to point geometries. Must be one of the following: 16 (no border), 21 (with border, default)
#' @param pt.stroke numeric, size of point borders. Only applied when pt.shape is set to 21 (default). Default is 0.1
#' @param pt.size numeric, point size to apply to centroids (when no polygons are available)
#' @param lg.pt.size numeric, size of points reported in plot legend (when no polygons are available)
#' @return plot object generated with ggplot
#' @importFrom viridis cividis viridis inferno plasma rocket magma mako turbo
#' @export
KanPlot = function(seurat = NULL,feature = NULL,fovs=NULL,
                      assay=DefaultAssay(seurat),layer='data',
                      palette_cont='inferno',palette_discrete=c24,alpha=NULL,show_border=F,border.size=NULL,border.col='black',max.cap=0.99,
                      bg_img=F,bg_type = c('img','mask'),max.px = 5e+05,theme=c('dark','light'),pt.shape = c(16,21),pt.stroke=0.1,
                   pt.size=0.1,lg.pt.size=2){
  if(bg_img == T){
    if(is.null(bg_type) | length(bg_type) > 1){
      bg_type = 'img'
   }
  }
  if(is.null(theme) | length(theme) > 1){
    theme = 'dark'
  }
  if(is.null(KanData(seurat))){
    stop('Seurat object does not contain any Kandinsky data You need to call kandinsky_init() before running this function')
  }
  if(is.null(feature)){
    stop('Please specify which feature you want to plot')
  }
  seurat = update_kandinsky(seurat)
  KanData(seurat,'sf') = populate_sf(seurat,vars=feature,assay=assay,layer=layer)
  if(stringr::str_detect(feature,'-')){
    colnames(KanData(seurat,'sf'))[colnames(KanData(seurat,'sf')) == gsub('-','\\.',feature)] = feature
  }
  if(is.numeric(sf::st_drop_geometry(KanData(seurat,'sf'))[[feature]]) ==F){
    palette = palette_discrete
    if(is.factor(sf::st_drop_geometry(KanData(seurat,'sf'))[[feature]]) == T){
      levels = sort(levels(sf::st_drop_geometry(KanData(seurat,'sf'))[[feature]]))
    }else{
      levels = sort(unique(sf::st_drop_geometry(KanData(seurat,'sf'))[[feature]]))
    }
    cols = palette[seq_len(length(levels))]
    names(cols) = levels
  }else{
    palette = palette_cont
    #vals = quantile(KanData(seurat,'sf')[[feature]],c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
    minval = quantile(KanData(seurat,'sf')[[feature]],0,na.rm=T)
    maxval = quantile(KanData(seurat,'sf')[[feature]],max.cap,na.rm=T)
    if(length(palette)==1){
      if(palette %in% c('magma','inferno','plasma','cividis','mako','viridis','turbo','rocket')){
        palette = eval(parse(text=paste0('viridis::',palette,'(n=9)')))
      }
    }
  }
  if(is.null(alpha)){
    alph = 1
  }else{
    alph = alpha
  }
  if(show_border ==T & inherits(KanData(seurat,'sf')$geometry,'sfc_POLYGON')){
    color = border.col
    lwd= border.size %||% 0.01
  }else{
    color =NA
    lwd = 0
  }
  if(!is.null(fovs)){
    seurat = suppressWarnings(seurat[,seurat$fov %in% fovs])
    seurat = update_kandinsky(seurat)
  }
  if(length(pt.shape) > 1){
    pt.shape = NULL
  }
  pt.shape = pt.shape %||% 21
  g=ggplot()
  if(bg_img == T){
    if(bg_type == 'mask'){
      rast = KanData(seurat,'mask')
    }else{
      rast = KanData(seurat,'img')
    }
    g = g + tidyterra::geom_spatraster_rgb(data=rast,alpha=ifelse(bg_type =='mask',1,0.8),maxcell=max.px)
  }
  if(inherits(KanData(seurat,'sf')$geometry,'sfc_POINT')){
    if(pt.shape==21){
      g=g+geom_sf(data=KanData(seurat,'sf'),mapping=aes(fill=.data[[feature]]),color='black',shape=pt.shape,stroke=pt.stroke,alpha=alph,size=pt.size)+
        theme_void()
      }else{
        g=g+geom_sf(data=KanData(seurat,'sf'),mapping=aes(color=.data[[feature]]),shape=pt.shape,alpha=alph,size=pt.size)+
          theme_void()
      }
    }else{
      g = g +geom_sf(data=KanData(seurat,'sf'),mapping=aes(fill=.data[[feature]]),shape=21,color=color,lwd=lwd,alpha=alph)+
      theme_void()
    }
  if(theme == 'dark'){
    g=g+theme(plot.background=element_rect(fill='black',colour='black'),
              legend.text = element_text(colour='white'),
              legend.title = element_text(colour='white'),
              legend.position='top')
  }else if(theme == 'light'){
    g=g+theme(plot.background=element_blank(),#element_rect(fill='snow',colour='snow'),
              legend.text = element_text(colour='black'),
              legend.title = element_text(colour='black'),
              legend.position='top')
  }
  if(!inherits(KanData(seurat,'sf')$geometry,'sfc_POINT') | pt.shape==21){
    if(is.numeric(sf::st_drop_geometry(KanData(seurat,'sf'))[,feature]) == F){
      g = g + scale_fill_manual(values = cols) + labs(fill=feature)
    }else{
      g = g + scale_fill_gradientn(colours=palette,na.value=palette[length(palette)],limits=c(minval,maxval))+ labs(fill = feature)
    }
  }
  if(inherits(KanData(seurat,'sf')$geometry,'sfc_POINT') & pt.shape==16){
    if(is.numeric(sf::st_drop_geometry(KanData(seurat,'sf'))[,feature]) == F){
      g = g + scale_color_manual(values = cols) + labs(color=feature)
    }else{
      g = g + scale_color_gradientn(colours=palette,na.value=palette[length(palette)],limits=c(minval,maxval))+ labs(color = feature)
    }
    }
  if(is.numeric(sf::st_drop_geometry(KanData(seurat,'sf'))[[feature]]) ==F){
    if(inherits(KanData(seurat,'sf')$geometry,'sfc_POINT') & pt.shape==16){
      g = g+guides(color=guide_legend(override.aes = list(size = lg.pt.size)))
      }else{
        g = g + guides(fill=guide_legend(override.aes = list(size = lg.pt.size)))
      }
    }
  g
}


#' @title Plot minimum distance between cell types as a spatial heatmap
#' @name fireflyPlot
#' @description
#' Calculate the shortest distance between cells belonging to two different groups, and visualize distance measures of the first cell type
#'
#' @details
#' For two groups A and B of cells of interest, this function will calculate the minimum distance from each cell of type A to cells of type B.
#' Find the closest cell of type B for each cell of type A, then plot the corresponding euclidean distance as spatial heatmap of cells (type A/B)
#' with the color scale based on the euclidean distance of cells A to the closest cell B, and with cells of type B colored in white to resemble a 'firefly' effect
#'
#' @param data a Seurat object containing Kandinsky data (`KanData()`)
#' @param label character string specifying the variable name to be used to defne cell annotation groups
#' @param ct1 cell type A label
#' @param ct2 cell type B label
#' @param ct1.lwd linewidth to use in plot for cell type A boundaries
#' @param ct2.lwd linewidth to use in plot for cell type B boundaries
#' @param fovs list of fov IDs to include in the analysis. If NULL, all fovs will be considered for the analysis
#' @param ct1.colthresh quantile threshold to apply to the distance colorscheme for cells belonging to cell type A in the spatial plot. For instance, if `ct1.coltresh = 0.99`, any distance greater than the 99 percentile will not result in a lighter or darker cell color
#' @returns a list object containing the data frame with the distance measurement plus a spatial plot reporting cells A colored according to their distance from cells B
#' @importFrom magrittr %>%
#' @export
fireflyPlot = function(data=NULL,label=NULL,ct1=NULL,ct1.lwd=0,ct2=NULL,ct2.lwd=0,fovs=NULL,ct1.colthresh = 0.99){
  data = populate_sf(data,vars=label,return.seurat = F)
  if(length(fovs)>0){
    data = data %>% dplyr::filter(.data$fov %in% fovs)
  }
  data[[label]] = stringr::str_replace(sf::st_drop_geometry(data)[[label]],' ','.')
  ct1.label = stringr::str_replace(ct1," ",".")
  ct2.label = stringr::str_replace(ct2," ",".")
  ct1 = data %>% dplyr::filter(.data[[label]] %in% ct1.label)
  ct2 = data %>% dplyr::filter(.data[[label]] %in% ct2.label)
  nearests = sf::st_nearest_feature(ct1$geometry,ct2$geometry)
  nearest_dist = sf::st_distance(x=ct1$geometry,y=ct2[nearests,]$geometry,by_element = T)
  ct1$dist = nearest_dist
  ct2$dist = 0
  if(is(ct1,'sf') == F){
    ct1 = sf::st_as_sf(ct1)
  }
  if(is(ct2,'sf') == F){
    ct2 = sf::st_as_sf(ct2)
  }
    g = ggplot(ct2)+
      geom_sf(aes(fill =.data[[label]]),lwd=ct2.lwd)+
      scale_fill_manual(values='ivory')+
      theme(legend.key = element_rect(fill = "white", colour = "black",linewidth=1))+
      labs(fill= '')+
      theme(legend.position = 'top')+
      guides(fill = guide_legend(title.position = "top"))+
      ggnewscale::new_scale_fill()+
      geom_sf(data= ct1,aes(fill=.data$dist),lwd=ct1.lwd)+
      scale_fill_gradientn(limits=c(0,quantile(ct1$dist,ct1.colthresh)),breaks=c(0,quantile(ct1$dist,ct1.colthresh)),labels=c('0','Max'),
                           na.value = 'midnightblue',
                           colors=c('#EEAA42','#FF6060','#3429C4','#1b129a','midnightblue')
      )+
      labs(fill=paste0(ct1.label,'-->',ct2.label,' Dist'))+
      theme(
        panel.background = element_rect(fill = "black",
                                        colour = "black",
                                        size  = 0.5,
                                        linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5,
                                        linetype = 'solid',
                                        colour = "black"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "black"),
        plot.background = element_rect(fill = 'black',
                                       colour = 'black'),
        legend.background = element_rect(fill='black'),
        legend.text = element_text(colour='white'),
        legend.title=element_text(colour='white'),
        legend.position = 'top'
      )+
      guides(fill = guide_colorbar(title.position = "top"))
    ct2.label = gsub("\\+","",ct2.label)
    colnames(ct1)[colnames(ct1) == 'dist'] = paste0(ct1.label,"_dist")
    return(list(dist = ct1,plot=g))
}


##Plot 2D kernel density field for each class type, like cell types or niches
#SpatialCellDensityPlot = function(seurat,cts=NULL,class=NULL,cols=c24,xcoord='x_global_px',ycoord='y_global_px'){
#  meta = st_as_sf(seurat@meta.data)
#  if(!is.null(cts)){
#    meta = meta[meta[[class]] %in% cts,]
#  }
#  ggplot(meta)+theme_minimal()+
#    geom_sf(aes(fill=.data[[class]]),colour=NA,show.legend = F)+
#    scale_fill_manual(values=cols)+
#    ggnewscale::new_scale_fill()+
#    stat_density2d_filled(data = meta,
#                          geom="density2d_filled",
#                          mapping=aes(x=.data[[xcoord]],
#                                      y=.data[[ycoord]],
#                                      color=.data[[class]],
#                                      alpha=after_stat(level)),
#                          contour=T,size=0.4)+
#    scale_color_manual(values=cols)+
#    theme(panel.background = element_rect(fill='black',colour='black'),
#          axis.line = element_line(colour='black'),
#          panel.grid.major=element_line(colour="black")
#    )+
#    facet_wrap(vars(class))+
#    guides(fill = 'none',alpha='none',colour='none')
#}


#' @title Proportion barplot
#' @name PropPlot
#' @description
#' Create a stacked barplot reporting the proportion covered by var.2 classes within each var.1 class
#' @param seurat a Seurat object containing Kandinsky data
#' @param var.1 character string, name of the first variable
#' @param var.2 character string, name of the second variable
#' @param cols color palette to use for plotting
#' @returns stacked barplot created with ggplot2
#'
#' @export
PropPlot = function(seurat=NULL,var.1=NULL,var.2=NULL,cols=NULL){
  data = FetchData(seurat,vars=c(var.1,var.2))
  cols = cols %||% c24
  data = data %>% dplyr::group_by(.data[[var.1]]) %>% dplyr::count(.data[[var.2]]) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=.data[[var.1]],y=.data[["n"]],fill=.data[[var.2]]))+
    theme_classic()+
    geom_bar(stat='identity',position='fill')+
    labs(y='Proportion')+
    scale_fill_manual(values=cols)+
    theme(axis.text.x=element_text(angle=45,hjust=1))
  data
}
#' @title Chi-square test residual plot
#' @name ResidualPlot
#' @description
#' Create a heatmap reporting standardized residuals from chi-square test between two variables of interest
#' @param seurat a Seurat object containing Kandinsky data
#' @param var.1 character string, name of the first variable
#' @param var.2 character string, name of the second variable
#' @param limits lower and upper residual caps for heatmap colour gradient
#' @param pval boolean, whether returning heatmap with cells annotated for significance or not (based on normal distribution fitting)
#' @param perm numeric, number of var.1 and var.2 reshuffling to estimate enrichment pvalue between var.1 and var.2 via permutation test. Default is set to 1000
#' @param padj.thresh numeric, adjusted pvalue threshold to define significant var.1/var.2 chisquare residuals. Default value is set to 0.01
#' @param alternative character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less".
#' @param return.df boolean, whether returning only residual plot (FALSE) or also a dataframe containing permutation-based chi-squared test results
#' @returns residual heatmap created with ggplot2
#'
#' @export
ResidualPlot = function(seurat=NULL,var.1=NULL,var.2=NULL,limits=c(-20,20),pval = T,perm=1000,padj.thresh=0.01,alternative = c("two.sided", "less", "greater"),return.df=F){
  data = FetchData(seurat,vars=c(var.1,var.2))
  chisq = chisq.test(table(data[[var.1]],data[[var.2]]))$residuals %>%
    as.data.frame()
  if(pval==T){
    if(length(alternative)>1){
      alternative = NULL
    }
    if(is.null(alternative)){
      message('"alternative" parameter is NULL. Setting "greater" as default argument')
      alternative = alternative %||% 'greater'
    }
      set.seed(347548)
      perms = lapply(seq_len(perm),function(x){
        return(chisq.test(table(sample(data[[var.1]],size=length(data[[var.1]])),
                                sample(data[[var.2]],size=length(data[[var.2]]))))$residuals %>%
                 as.data.frame() %>%
                 dplyr::mutate(real_freq = chisq$Freq) %>%
                 dplyr::mutate(higher = ifelse(.data[["Freq"]] > .data[["real_freq"]],1,0),
                               lower = ifelse(.data[["Freq"]] < .data[["real_freq"]],1,0),
                               two = ifelse(abs(.data[["Freq"]]) > abs(.data[["real_freq"]]),1,0)))
      })
      perms = purrr::reduce(perms,rbind) %>%
        dplyr::group_by(.data[["Var1"]],.data[["Var2"]]) %>%
        dplyr::mutate(pval_higher=sum(.data[["higher"]])/perm,
                      pval_lower = sum(.data[["lower"]])/perm,
                      pval_two = sum(.data[["two"]])/perm) %>%
        dplyr::ungroup() %>% dplyr::select(.data[["Var1"]],.data[["Var2"]],.data[["pval_higher"]],.data[["pval_lower"]]) %>%
        unique()
      perms$padj_higher = stats::p.adjust(perms$pval_higher,method='BH')
      perms$padj_lower = stats::p.adjust(perms$pval_lower,method='BH')
      if(alternative=='two.sided'){
        #perms$pval_two =2 * pmin(perms$pval_higher, perms$pval_lower)
        perms$padj_two = stats::p.adjust(perms$pval_two,method='BH')
      }
  #chisq = reshape2::melt(chisq.test(table(data[[var.1]],data[[var.2]]))$residuals) %>%
  #  dplyr::mutate(pval_g = stats::pnorm(value,mean=0,sd=1,lower.tail = F),
  #         pval_l = stats::pnorm(value,mean=0,sd=1,lower.tail = T)) %>%
  #  dplyr::mutate(pval_two = 2*stats::pnorm(abs(value),mean=0,sd=1,lower.tail = F))#2*pmin(pval_g,pval_l))
  #colnames(chisq)[c(1,2)] = c(var.1,var.2)
  #chisq = chisq %>% dplyr::mutate(padj_g = p.adjust(pval_g,method='BH'),
  #                           padj_l = p.adjust(pval_l,method='BH'),
  #                           padj_two = p.adjust(pval_two,method='BH'))
  }
  if(pval==T){
  if(alternative == 'greater'){
    perms$padj = perms$padj_higher
  }
  if(alternative == 'less'){
    perms$padj = perms$padj_lower
  }
  if(alternative=='two.sided'){
    perms$padj = perms$padj_two
  }
  chisq = merge(chisq,perms[,c('Var1','Var2','padj')],by=c('Var1','Var2'))
  chisq$label = ifelse(chisq$padj < padj.thresh,"*","")
  colnames(chisq)[c(1,2,3)] = c(var.1,var.2,'value')
  }
  g = ggplot(chisq,aes(x=.data[[var.1]],y=.data[[var.2]],fill=.data[["value"]]))+
    theme_classic()+geom_tile(color='black')+theme(axis.text.x=element_text(angle=45,hjust=1),axis.text=element_text(color='black'),
                                                   axis.line = element_blank(),axis.ticks = element_blank())+
    labs(fill='Std. residuals')+labs(x=var.1,y=var.2)+coord_fixed(ratio = 1)
    if(pval ==T){
      if(alternative=='greater'){
        g= g+scale_fill_gradient2(low='white',high='red',midpoint=0,limits=limits, oob = scales::squish)
      }else if(alternative=='less'){
      g= g+scale_fill_gradient2(low='blue',high='white',midpoint=0,limits=limits, oob = scales::squish)
      }else{
      g= g+scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,limits=limits, oob = scales::squish)
      }
    }else{
      g= g+scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,limits=limits, oob = scales::squish)
    }
  if(pval==T){
  g= g+geom_text(stat='identity',aes(label=.data[["label"]]),size=4)
  }
  if(return.df==F){
  g
  }else{
    return(list(res = chisq,plot=g))
  }
}

#' @title Combined PropPlot and ResidualPlot
#' @name EnrichPlots
#' @description
#' Merged visualization of proportion barplot and Chi-square residual heatmap between two variables of interest
#' @param seurat a Seurat object containing Kandinsky data
#' @param var.1 character string, name of the first variable
#' @param var.2 character string, name of the second variable
#' @param limits lower and upper residual caps for heatmap colour gradient
#' @param cols color palette to use for plotting
#' @param pval boolean, whether returning heatmap with cells annotated for significance or not (based on permutation test)
#' @param perm numeric, number of var.1 and var.2 reshuffling to estimate enrichment pvalue between var.1 and var.2 via permutation test. Default is set to 1000
#' @param padj.thresh numeric, adjusted pvalue threshold to define significant var.1/var.2 chisquare residuals
#' @param alternative character string specifying the alternative hypothesis, must be one of "two.sided", "greater" (default) or "less".
#' @export
EnrichPlots = function(seurat=NULL,var.1=NULL,var.2=NULL,limits=c(-20,20),cols=NULL,pval=T,perm=1000,padj.thresh=0.05,alternative = c("two.sided", "less", "greater")){
  g1 = PropPlot(seurat=seurat,var.1=var.1,var.2=var.2,cols=cols)
  g2 = ResidualPlot(seurat=seurat,var.1=var.1,var.2=var.2,limits=limits,pval=pval,perm=perm,padj.thresh=padj.thresh,alternative=alternative)
  patchwork::wrap_plots(list(g1,g2),guides = 'collect',ncol=2)+
    patchwork::plot_layout(axis_titles='collect')
}

#' @title density plot of neighbourhood size across cells/spots
#' @name nbSizePlot
#' @description
#' Visualize the distribution of possible neighbourhood sizes observed across cells/spots within the dataset
#' @param seurat a Seurat object containing Kandinsky data
#' @returns density plot of cell/spot neighbourhood sizes
#' @export
nbSizePlot = function(seurat){
  tot_nb = data.frame(size = Matrix::colSums(as(KanData(seurat,'nb'),'CsparseMatrix')))
  labels = unique(c(min(tot_nb$size),median(tot_nb$size),max(tot_nb$size)))
  ggplot(tot_nb,aes(y=.data[["size"]]))+theme_classic()+geom_density()+scale_y_reverse(breaks=(labels),labels=(labels))+
    annotate('segment',x=0,xend= max(density(tot_nb$size)$y), y= median(tot_nb$size),yend=median(tot_nb$size),linetype='dashed',color='red')+theme(axis.text.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x = element_blank())+
    guides(y=guide_axis(cap = "both"))+labs(y='c-NB size',x='')
}




