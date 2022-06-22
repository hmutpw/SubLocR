


#############
# subcellular localization visualization
#############

#' Subcellular visualization
#'
#' @param gene_name Gene symbols used for visualization. At least two genes
#' @param stage Stages used for visualization
#' @param use_data use normalized data?
#' @param plot_type type of plot?
#' @param color_by which column used for visualization?
#' @param show_gene_name Whether to show gene names?
#' @param label_size The gene names size in figure.
#' @param point_size The point size.
#' @param point_shape The point shape.
#' @param continuous_colors The color of points.
#' @param text_size_circle the text color of outlier.
#' @param text_color_circle the color of outlier.
#' @param ... other parameters.
#'
#' @return a ggplot2 object
#' @importFrom purrr map2
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export
#'
#' @examples
#' subLocViz(gene_name=c("RUNX1","MYC"))
#'
subLocViz <- function(gene_name,
                      stage = c("D0","D3","D5"),
                      use_data = c("scaled","ratio","raw"),
                      plot_type = c("point","density","hexagonal","buble"),
                      color_by = NULL,
                      show_gene_name = FALSE,
                      label_size = 2,
                      point_size = 1,
                      point_shape = 16,
                      continuous_colors = continuous_colorbar(),
                      text_size_circle = 3,
                      text_color_circle = "purple",
                      ...){
  #---data preparation
  data(ESC_RSL_data)
  anchors <- c('CE',"ME","SNE","CBNE","PE")
  stage <- match.arg(arg = stage, several.ok = TRUE)
  use_data <- match.arg(arg = use_data)
  plot_type <- match.arg(arg = plot_type)
  data_list <- ESC_RSL_data[stage]
  data_list <- lapply(X = data_list, FUN = function(x, gene){
    x[x$gene_name %in% gene, ,drop=FALSE]
  }, gene = gene_name)
  if(use_data == "ratio"){
    data_list <- lapply(X = data_list, FUN = function(x, anchor){
      ratio <- t(apply(x[,anchors, drop=FALSE],1,function(x){x/sum(x)}))
      x[,anchors] <- ratio
      x
    },anchor = anchors)
  }else if(use_data == "scaled"){
    data_list <- lapply(X = data_list, FUN = function(x, anchor){
      scaled <- t(apply(x[,anchors, drop=FALSE],1,function(x){(x-min(x))/(max(x)-min(x))}))
      x[,anchors] <- scaled
      x
    },anchor = anchors)
  }
  #---data visualization
  rv_list <- lapply(X = data_list,
                    FUN = getRadObj,
                    anchors = anchors,
                    norm_data = FALSE,
                    optim_anchor = FALSE,
                    text_color = text_color_circle,
                    text_size = text_size_circle,
                    ...)
  Radviz_list <- lapply(X = rv_list,
                        FUN = do_RadvizEnrichPlot,
                        color_by = color_by,
                        outline_circle = TRUE,
                        plot_type = plot_type,
                        point_size = point_size,
                        point_shape = point_shape,
                        ...)
  if(show_gene_name){
    gg_lab <- geom_text(mapping = aes(label=gene_name),size = label_size)
    Radviz_list <- lapply(X = Radviz_list,
                          FUN = function(x,gg_lab){
                            x+gg_lab
                          },gg_lab=gg_lab)
  }
  #---add title
  titles <- names(Radviz_list)
  Radviz_list <- purrr::map2(.x = Radviz_list, .y = titles, .f = function(x,y){
    x+ggtitle(label = y)+theme(plot.title = element_text(hjust = 0.5,size = 12))
  })
  ggpubr::ggarrange(plotlist = Radviz_list, nrow = 1)
}
