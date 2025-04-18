% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IObenchmark.R
\name{SCE2AD}
\alias{SCE2AD}
\title{Convert SingleCellExperiemnt to AnnData using sceasy/scDIOR/zellkonverter.}
\usage{
SCE2AD(
  sce.obj,
  method = c("sceasy", "scDIOR", "zellkonverter"),
  out.folder = NULL,
  out.filename = NULL,
  slot = "counts",
  conda.path = NULL
)
}
\arguments{
\item{sce.obj}{A SingleCellExperiment object.}

\item{method}{Method used to perform conversion, choose from "sceasy",    "scDIOR", "zellkonverter". Default: "sceasy".}

\item{out.folder}{Output folder. Default: NULL (current working directory).}

\item{out.filename}{Output file name, e.g., test.h5ad. Default: NULL (\code{seu.obj} name, method, h5ad.).}

\item{slot}{Slot for adata.X, used when \code{method} is "sceasy", "zellkonverter". Default: counts.}

\item{conda.path}{Conda environment path, used when \code{method} is "zellkonverter" or "sceasy". Default: NULL.}
}
\value{
Run log.
}
\description{
Convert SingleCellExperiemnt to AnnData using sceasy/scDIOR/zellkonverter.
}
\examples{
# sceasy
SCE2AD(sce.obj = pbmc3k.sce, method = "sceasy", out.folder = "benchmark", slot = "rawcounts", conda.path = "/path/to/conda")
# scDIOR
pbmc3k.sce.scdior <- pbmc3k.sce
library(SingleCellExperiment)
# scDIOR does not support varm in rowData
rowData(pbmc3k.sce.scdior)$varm <- NULL
SCE2AD(sce.obj = pbmc3k.sce.scdior, method = "scDIOR", out.folder = "benchmark")
# zellkonverter
SCE2AD(sce.obj = pbmc3k.sce, method = "zellkonverter", out.folder = "benchmark", slot = "rawcounts", conda.path = "/path/to/conda")
}
