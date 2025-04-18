% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zenodo.R
\name{ParseZenodo}
\alias{ParseZenodo}
\title{Download Data with Zenodo DOI and Return SeuratObject.}
\usage{
ParseZenodo(
  doi = NULL,
  file.ext = c("rdata", "rds", "h5ad"),
  doi.df = NULL,
  out.folder = NULL,
  timeout = 1000,
  quiet = FALSE,
  parallel = TRUE,
  use.cores = NULL,
  return.seu = FALSE,
  merge = TRUE
)
}
\arguments{
\item{doi}{A vector of Zenodo DOIs to download. Default: NULL.}

\item{file.ext}{The valid file extension for download. When NULL, use all files. Default: c("rdata", "rds", "h5ad").}

\item{doi.df}{DOI dataframe for download. This is useful when something wrong happens in downloading
(e.g. MD5 verification failure, \code{DownloadZenodo} will return a dataframe contains failure terms.). Default: NULL.
It is required to provide either \code{doi} or \code{doi.df}.}

\item{out.folder}{The output folder. Default: NULL (current working directory).}

\item{timeout}{Maximum request time. Default: 1000.}

\item{quiet}{Logical value, whether to show downloading progress. Default: FALSE (show).}

\item{parallel}{Logical value, whether to download parallelly. Default: TRUE. When "libcurl" is available for \code{download.file},
the parallel is done by default (\code{parallel} can be FALSE).}

\item{use.cores}{The number of cores used. Default: NULL (the minimum value of
\code{nrow(doi.df)} (specified by \code{doi.df} or extract with \code{doi}) and \code{parallel::detectCores()}).}

\item{return.seu}{Logical value, whether to load downloaded datasets to Seurat. Valid when rds in \code{file.ext} and all
datasets download successfully. Default: FALSE.}

\item{merge}{Logical value, whether to merge Seurat list when there are multiple rds files,
used when \code{return.seu} is TRUE. Default: FALSE.}
}
\value{
Dataframe contains failed datasets, SeuratObject (\code{return.seu} is TRUE, rds in \code{file.ext}) or
NULL (\code{return.seu} is FALSE or rds not in \code{file.ext}).
}
\description{
Download Data with Zenodo DOI and Return SeuratObject.
}
\examples{
\dontrun{
# need users to provide the output folder
multi.dois.parse <- ParseZenodo(
  doi = c(
    "1111", "10.5281/zenodo.7243603",
    "10.5281/zenodo.7244441"
  ),
  file.ext = c("rdata", "rds"),
  out.folder = "/path/to/outfoder"
)
}
}
