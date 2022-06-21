
#' Generate Radviz object for visualization
#'
#' @param input_frm A data.frame used to generate the Radviz object. The columns were items around the circle.
#' @param anchors Dimensions used to generate the map. Default: all the columns from the \code{input_frm}.
#' @param norm_data Whether to normalized data into 0 to 1? Default: TRUE.
#' @param trim_ratio The percentage of values trimmed for every anchor. Used when \code{norm_data} is TRUE. Default: 0.005.
#' @param optim_anchor Whether to optimize the anchor order? Default: TRUE.
#' @param recenter User defined the start anchor. Default: NULL.
#' @param text_color Text color of anchors, inherited from \code{\link{do.radviz}}. Default: black.
#' @param text_size Text size of anchors, inherited from \code{\link{do.radviz}}. Default:8.
#' @param ... Other parameters passed to \code{\link{do.radviz}}.
#'
#' @return getRadviz object
#'
#' @export
#' @importFrom Radviz do.L make.S cosine do.optimRadviz get.optim recenter do.radviz

getRadObj <- function(input_frm,
                      anchors = colnames(input_frm),
                      norm_data = TRUE,
                      trim_ratio = 0.005,
                      optim_anchor = TRUE,
                      recenter = NULL,
                      text_color = "black",
                      text_size = 4, ...){
  # check anchors
  anchors <- intersect(anchors,colnames(input_frm))
  if(!all(anchors %in% colnames(input_frm))){
    not_in <- setdiff(anchors, colnames(input_frm))
    stop("Anchors [",paste(not_in, collapse = ","),"] not found from columns of input_frm!")
  }
  exp_mat <- as.matrix(input_frm[,anchors])
  if(!is.numeric(exp_mat)){
    stop("The 'anchors' columns in 'input_frm' must be numeric!")
  }
  # data normalization
  trans <- function(vec, cutoff=trim_ratio){
    Radviz::do.L(vec, fun = function(x) quantile(x,c(cutoff,1-cutoff)))
  }
  if(norm_data){
    exp_mat <- t(apply(exp_mat,1,function(x){(x-min(x))/(max(x)-min(x))}))
    colnames(exp_mat) <- anchors
    other_cols <-setdiff(colnames(input_frm), anchors)
    if(length(other_cols)>0){
      input_frm <- data.frame(exp_mat, input_frm[, other_cols, drop=FALSE], check.names=F)
    }else{
      input_frm <- exp_mat
    }
  }
  #---generate springs
  mat.S <- Radviz::make.S(anchors)
  if(optim_anchor){
    mat.sim <- Radviz::cosine(exp_mat)
    optim.mat <- Radviz::do.optimRadviz(mat.S,mat.sim,iter=10,n=100)
    mat.S <- Radviz::make.S(Radviz::get.optim(optim.mat))
    if(!is.null(recenter) && length(recenter)==1 && (recenter %in% anchors)){
      mat.S <- Radviz::recenter(mat.S, newc = recenter)
    }
  }
  rv <- Radviz::do.radviz(x = input_frm,
                          springs = mat.S,
                          trans = trans,
                          label.color = text_color,
                          label.size = text_size,...)
  rv
}


#############
# add circle for Radviz
#############
addRadvizCircle <- function(rv, color="black"){
  circle <- ggplot2::annotate("path", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi,length.out=100)), color = color)
  rv$proj <- rv$proj+circle
  rv
}

#############
# do Radviz plot
#############
#' Radviz visualization
#' @param rv A \code{\link{getRadviz}} object generate from \code{\link{do.radviz}}.
#' @param color_by Groups used to label the color.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot.
#' @param text_size Text size.
#' @param point_size Label size of point in \code{\link{geom_point()}}.
#' @param point_shape Label shape of point in \code{\link{geom_point()}}.
#' @param ... Other input parameters.
#'
#' @return a ggplot object
#'
#' @import ggplot2
#'
#' @export
do_RadvizEnrichPlot <- function(rv,
                                color_by = NULL,
                                outline_circle=TRUE,
                                plot_type=c("point","density","hexagonal","buble"),
                                point_size=0.5,
                                point_shape=16,
                                alpha=1,
                                continuous_colors = continuous_colorbar(),
                                ...){
  plot_type <- match.arg(plot_type)
  rv_data_frm <- as.data.frame(rv$proj$data)
  if(outline_circle) rv <- addRadvizCircle(rv)
  if(!is.null(color_by)){
    color_by <- intersect(color_by, colnames(rv_data_frm))
    if(length(color_by)!=1) stop("The 'color_by' shold only have one element!")
  }
  p <- rv$proj
  if(plot_type=="point"){
    #---point plot
    if(is.character(point_shape) && point_shape %in% colnames(rv_data_frm)){
      p <- p + geom_point(aes_string(color=color_by, shape=point_shape),
                          size=point_size, alpha=alpha, na.rm = T)
    }else if(is.numeric(point_shape)){
      p <- p + geom_point(aes_string(color=color_by), shape=point_shape,
                          size=point_size, alpha=alpha, na.rm = T)
    }
  }else if(plot_type=="density"){
    #---density
    p <- Radviz::smoothRadviz(rv)+geom_point(shape='.', alpha=alpha)
  }else if(plot_type=="hexagonal"){
    #---hexagonal
    p <- Radviz::hexplot(rv)
  }else if(plot_type=="buble"){
    #---buble
    p <- Radviz::bubbleRadviz(rv, group=color_by)
  }
  if(!is.null(color_by) && mode(rv_data_frm[[color_by]]) == "numeric"){
    p <- p+ggplot2::scale_color_gradientn(colours = rev(continuous_colors))
  }
  p
}

##############
#---calculate entropy from probability
#############
.H <- function(freqs, unit = c("log", "log2", "log10")){
  unit <- match.arg(unit)
  #freqs = freqs/sum(freqs)
  H <- -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
  if (unit == "log2") H <- H/log(2)
  if (unit == "log10") H <- H/log(10)
  H
}


##############
# color bars
##############

continuous_colorbar <- function(type = c('RdYlBu', 'RdYlGn',"PRGn","PiYG", "Oranges", "Purples", "Greens", "Greys"),
                                n = 100){
  type <- match.arg(arg = type)
  col_num <- RColorBrewer::brewer.pal.info[type, 'maxcolors']
  bre <- RColorBrewer::brewer.pal(n = col_num, name = type)
  pal <- colorRampPalette(colors = bre)(n)
  pal
}



