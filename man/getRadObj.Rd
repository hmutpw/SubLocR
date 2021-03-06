% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Radviz.R
\name{getRadObj}
\alias{getRadObj}
\title{Generate Radviz object for visualization}
\usage{
getRadObj(
  input_frm,
  anchors = colnames(input_frm),
  norm_data = TRUE,
  trim_ratio = 0.005,
  optim_anchor = TRUE,
  recenter = NULL,
  text_color = "black",
  text_size = 4,
  ...
)
}
\arguments{
\item{input_frm}{A data.frame used to generate the Radviz object. The columns were items around the circle.}

\item{anchors}{Dimensions used to generate the map. Default: all the columns from the \code{input_frm}.}

\item{norm_data}{Whether to normalized data into 0 to 1? Default: TRUE.}

\item{trim_ratio}{The percentage of values trimmed for every anchor. Used when \code{norm_data} is TRUE. Default: 0.005.}

\item{optim_anchor}{Whether to optimize the anchor order? Default: TRUE.}

\item{recenter}{User defined the start anchor. Default: NULL.}

\item{text_color}{Text color of anchors, inherited from \code{\link{do.radviz}}. Default: black.}

\item{text_size}{Text size of anchors, inherited from \code{\link{do.radviz}}. Default:8.}

\item{...}{Other parameters passed to \code{\link{do.radviz}}.}
}
\value{
getRadviz object
}
\description{
Generate Radviz object for visualization
}
