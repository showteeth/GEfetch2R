#' Download Matrix from GEO and Load to Seurat/DESeq2.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Disable when \code{down.supp} is TRUE. Default: NULL (disable).
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or empty,
#' download supplementary files automatically). Default: FALSE.
#' @param supp.idx The index of supplementary files to download. This should be consistent with \code{platform}. Used when \code{supp.type} is 10x/count. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param data.type The data type of the dataset, choose from "sc" (single-cell) and "bulk" (bulk). Default: "sc".
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param file.regex The regex to extract correct count matrix files. Used when \code{supp.type} is count.
#' Default: NULL (when multiple file extensions are available in the downloaded tar file, use the first file extension).
#' @param extra.cols Extra columns to remove, e.g., "Chr", "Start", "End", "Strand", "Length" (featureCounts). Used when \code{supp.type} is count.
#' Default: "chr", "start", "end", "strand", "length", "width", "chromosome", "seqnames", "seqname", "chrom", "chromosome_name", "seqid", "stop".
#' @param transpose Logical value, whether to transpose the matrix. Used when the number of rows is less than the number of columns. Used when \code{supp.type} is count. Default: TRUE.
#' @param out.folder Output folder to save 10x files. Used when \code{supp.type} is 10x/10xSingle. Default: NULL (current working directory).
#' @param accept.fmt Vector of accepted 10x output format, MEX (barcode/feature/gene/matrix), h5 (h5/h5.gz). Used when \code{supp.type} is 10x/10xSingle. Default: c("MEX", "h5").
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Used when \code{supp.type} is 10x/10xSingle. Default: TRUE.
#' @param load2R Logical value, whether to load the count matrix to R. Default: TRUE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple 10x files (\code{supp.type} is 10x).
#' Used when \code{supp.type} is 10x/10xSingle. Default: FALSE.
#' @param meta.data Dataframe contains sample information for DESeqDataSet. Used when \code{data.type} is bulk. Default: NULL.
#' @param fmu Column of \code{meta.data} contains group information. Used when \code{data.type} is bulk. Default: NULL.
#' @param ... Parameters for \code{\link{getGEO}}. Used when \code{down.supp} is FALSE.
#'
#' @return If \code{load2R} is FALSE, return count matrix. If \code{data.type} is "sc", return SeuratObject (if \code{merge} is TRUE), SeuratObject list (if \code{merge} is FALSE), NULL (no SeuratObject detected).
#' If \code{data.type} is "bulk", return DESeqDataSet.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO getGEOSuppFiles
#' @importFrom xml2 read_html xml_text xml_find_all
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom tools file_ext
#' @importFrom utils untar
#' @importFrom R.utils gzip
#' @importFrom data.table fread
#' @importFrom readxl read_excel
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom methods new
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @export
#'
#' @examples
#' \dontrun{
#' # the supp files are count matrix
#' GSE94820.seu <- ParseGEO(acce = "GSE94820", down.supp = TRUE, supp.idx = 1, supp.type = "count")
#' # the supp files are cellranger output files: barcodes, genes/features and matrix
#' # need users to provide the output folder
#' GSE200257.seu <- ParseGEO(
#'   acce = "GSE200257", down.supp = TRUE, supp.idx = 1, supp.type = "10x",
#'   out.folder = "/path/to/output/folder"
#' )
#' # need users to provide the output folder
#' GSE236082.seu <- ParseGEO(
#'   acce = "GSE236082", down.supp = TRUE, supp.type = "10xSingle",
#'   out.folder = "/path/to/output/folder"
#' )
#' }
ParseGEO <- function(acce, platform = NULL, down.supp = FALSE, supp.idx = 1, timeout = 3600, data.type = c("sc", "bulk"),
                     supp.type = c("count", "10x", "10xSingle"), file.regex = NULL,
                     extra.cols = c(
                       "chr", "start", "end", "strand", "length",
                       "width", "chromosome", "seqnames", "seqname",
                       "chrom", "chromosome_name", "seqid", "stop"
                     ),
                     transpose = TRUE, out.folder = NULL, accept.fmt = c("MEX", "h5"), gene2feature = TRUE,
                     load2R = TRUE, merge = TRUE, meta.data = NULL, fmu = NULL, ...) {
  # check parameters
  data.type <- match.arg(arg = data.type)
  supp.type <- match.arg(arg = supp.type)
  # check platform
  if (down.supp) {
    message("Download supplementary files to generate matrix!")
    pf.obj <- NULL
  } else {
    message("Extract expression data from eSets!")
    if (is.null(platform)) {
      stop("Platform is required to extract expression data!")
    }
    # get GEO object
    pf.obj <- GEOobj(acce = acce, platform = platform, ...)
  }
  # change supp type to count when bulk
  if (data.type == "bulk") {
    supp.type <- "count"
  }
  # extract counts matrix
  pf.count <- ExtractGEOExp(
    pf.obj = pf.obj, acce = acce, supp.idx = supp.idx, down.supp = down.supp,
    timeout = timeout, supp.type = supp.type, file.regex = file.regex,
    extra.cols = extra.cols, transpose = transpose,
    out.folder = out.folder, accept.fmt = accept.fmt, gene2feature = gene2feature
  )
  # load to R
  if (load2R) {
    if (data.type == "bulk") {
      de.obj <- Loading2DESeq2(mat = pf.count, meta = meta.data, fmu = fmu)
      return(de.obj)
    } else if (data.type == "sc") {
      # load seurat
      # if (is.null(pf.count) && supp.type == "10x") {
      if (is.null(pf.count) && (supp.type == "10x" || supp.type == "10xSingle")) {
        message("Loading data to Seurat!")
        out.folder <- file.path(out.folder, acce) # optimize output folder
        seu.list <- list()
        # check MEX format
        if ("MEX" %in% accept.fmt) {
          all.bc.files <- list.files(path = out.folder, pattern = "barcodes.tsv.gz$", recursive = TRUE)
          if (length(all.bc.files) > 0) {
            all.samples.folder <- dir(out.folder, full.names = TRUE)
            # check file
            valid.samples.folder <- Check10XFiles(folders = all.samples.folder, gene2feature = gene2feature)
            if (length(valid.samples.folder) > 0) {
              # load to seurat
              seu.list.mex <- sapply(valid.samples.folder, function(x) {
                x.mat <- Seurat::Read10X(data.dir = x)
                if (class(x.mat) == "list") {
                  if ("Gene Expression" %in% names(x.mat)) {
                    message("10X data contains more than one type: ", paste(names(x.mat), collapse = ", "), ". Load Gene Expression to Seurat!")
                    x.mat <- x.mat[["Gene Expression"]]
                  } else {
                    message("10X data contains more than one type: ", paste(names(x.mat), collapse = ", "), ", and no 'Gene Expression', load the first element!")
                    x.mat <- x.mat[[1]]
                  }
                }
                seu.obj <- Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
                seu.obj
              })
            } else {
              message("No valid sample folder (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz/genes.tsv.gz) detected under ", out.folder, ". Please check!")
              seu.list.mex <- list()
            }
          } else {
            message("No barcodes.tsv.gz files detected under: ", out.folder, ". Please check!")
            seu.list.mex <- list()
          }
          seu.list <- c(seu.list, seu.list.mex)
        }
        # check h5 format
        if ("h5" %in% accept.fmt) {
          all.h5.files <- list.files(path = out.folder, pattern = "h5$", recursive = TRUE, full.names = TRUE)
          if (length(all.h5.files) > 0) {
            # load to seurat
            seu.list.h5 <- sapply(all.h5.files, function(x) {
              x.mat <- Seurat::Read10X_h5(filename = x)
              sample.name <- gsub(pattern = ".h5$", replacement = "", x = basename(x))
              seu.obj <- Seurat::CreateSeuratObject(counts = x.mat, project = sample.name)
              seu.obj
            })
          } else {
            message("No h5 files detected under: ", out.folder, ". Please check!")
            seu.list.h5 <- list()
          }
          seu.list <- c(seu.list, seu.list.h5)
        }
        if (length(seu.list) > 0) {
          if (isTRUE(merge)) {
            seu.obj <- mergeExperiments(seu.list)
          } else {
            seu.obj <- seu.list
          }
        } else {
          message("No SeuratObject detected, please check!")
          seu.obj <- NULL
        }
      } else if (!is.null(pf.count) && supp.type == "count") {
        seu.obj <- Seurat::CreateSeuratObject(counts = pf.count, project = acce)
      }
      return(seu.obj)
    }
  } else {
    return(pf.count)
  }
}

#' Download Processed Objects from GEO.
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param file.ext The valid file extension for download (ignore case and gz suffix). When NULL, use all files. Default: c("rdata", "rds", "h5ad").
#' @param out.folder The output folder. Default: NULL (\code{acce} folder under current working directory).
#' @param return.seu Logical value, whether to load downloaded datasets to Seurat. Valid when rds in \code{file.ext} and all
#' datasets download successfully. Default: FALSE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple rds files,
#' used when \code{return.seu} is TRUE. Default: FALSE.
#'
#' @return SeuratObject (\code{return.seu} is TRUE, rds in \code{file.ext}) or
#' NULL (\code{return.seu} is FALSE or rds not in \code{file.ext}).
#' @importFrom xml2 read_html xml_text xml_find_all
#' @importFrom tools file_ext
#' @importFrom utils untar
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom rlang parse_expr
#' @importFrom methods new
#' @export
#'
#' @examples
#' \dontrun{
#' # need users to provide the output folder
#' # suitable for rdata, rdata.gz, rds, rds.gz
#' # return SeuratObject when return.seu = TRUE
#' GSE285723.seu <- ParseGEOProcessed(
#'   acce = "GSE285723", supp.idx = 1,
#'   file.ext = c("rdata", "rds"), return.seu = TRUE,
#'   out.folder = "/path/to/outfoder"
#' )
#' geo.rdata.log <- ParseGEOProcessed(
#'   acce = "GSE311825", supp.idx = 4,
#'   file.ext = c("rdata", "rds"),
#'   out.folder = "/path/to/outfoder"
#' )
#' # suitable for h5ad, h5ad.gz
#' geo.h5ad.log <- ParseGEOProcessed(
#'   acce = "GSE311813", supp.idx = 1,
#'   file.ext = c("h5ad"),
#'   out.folder = "/path/to/outfoder"
#' )
#' # suitable for loom, loom.gz
#' geo.loom.log <- ParseGEOProcessed(
#'   acce = "GSE286325", supp.idx = 1,
#'   file.ext = c("loom"),
#'   out.folder = "/path/to/outfoder"
#' )
#' }
ParseGEOProcessed <- function(acce, timeout = 3600, supp.idx = 1, file.ext = c("rdata", "rds", "h5ad", "loom"),
                              out.folder = NULL, return.seu = FALSE, merge = TRUE) {
  # file.ext: ignore case, tar.gz, gz
  if (is.null(file.ext)) {
    warning("There is no file extension provided, use all valid (rdata, rds, h5ad and loom).")
    file.ext <- c("rdata", "rds", "h5ad", "loom")
  }
  file.ext <- intersect(file.ext, c("rdata", "rds", "h5ad", "loom"))
  if (length(file.ext) == 0) {
    stop("Please provide valid file extension: rdata, rds, h5ad and loom.")
  }
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # set output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  out.folder <- file.path(out.folder, acce) # optimize output folder
  # create folder
  if (!dir.exists(out.folder)) {
    dir.create(path = out.folder, recursive = TRUE)
  }
  # download supplementary file
  supp.down.log <- tryCatch(
    expr = {
      getGEOSuppFilesInner(GEO = acce, baseDir = tmp.folder, index = supp.idx)
    },
    error = function(e) {
      print(e)
      stop("Please check the supp.idx or change the timeout with a larger value.")
    }
  )
  supp.file.path <- rownames(supp.down.log)
  # file unzip
  fext <- tools::file_ext(supp.file.path)
  if (fext == "gz") {
    # gunzip file
    Gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    fext <- tools::file_ext(supp.file.path)
  }
  if (fext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # unzip, move to given format
    all.gz.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "gz$")
    if (length(all.gz.files) > 0) {
      message("Detect files in gz format, gunzip!")
      # unzip
      unzip.log <- sapply(
        all.gz.files,
        function(x) {
          Gunzip(x, overwrite = TRUE)
        }
      )
    }
    valid.pat <- paste0(file.ext, "$", collapse = "|")
    all.used.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat, ignore.case = TRUE)
    if (length(all.used.files) > 0) {
      rename.log <- sapply(all.used.files, function(x) {
        # move file
        copy.tag <- file.copy(from = x, to = file.path(out.folder, basename(x)))
        # remove the original file
        remove.tag <- file.remove(x)
        copy.tag
      })
      # load RDS file to Seurat
      if (isTRUE(return.seu)) {
        seu.obj <- LoadRDS2Seurat(out.folder = out.folder, merge = merge)
        return(seu.obj)
      } else {
        return(NULL)
      }
    } else {
      stop("No file in given format, please check file.ext!")
    }
  } else if (tolower(fext) %in% tolower(file.ext)) {
    # move file
    copy.tag <- file.copy(from = supp.file.path, to = file.path(out.folder, basename(supp.file.path)))
    # remove the original file
    remove.tag <- file.remove(supp.file.path)
    # load RDS file to Seurat
    if (isTRUE(return.seu)) {
      seu.obj <- LoadRDS2Seurat(out.folder = out.folder, merge = merge)
      return(seu.obj)
    } else {
      return(NULL)
    }
  } else {
    stop("No file in given format, please check file.ext!")
  }
}

#' Extract Sample Metadata from GEO.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Default: NULL (all platforms).
#' @param down.supp Logical value, whether to download supplementary files to extract sample metadata. If TRUE, always
#' download supplementary files. If FALSE, extract GEO sample metadata. Default: FALSE.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return Dataframe contains all metadata of provided GEO accession number.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom tools file_ext
#' @importFrom readxl read_excel
#' @importFrom data.table fread
#' @importFrom xml2 read_html xml_text xml_find_all
#' @export
#'
#' @examples
#' \donttest{
#' # users may need to set the size of the connection buffer
#' # Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
#' # extract GEO sample metadata of specified platform
#' GSE200257.meta <- ExtractGEOMeta(acce = "GSE200257", platform = "GPL24676")
#' # extract sample metadata from supplementary files
#' GSE292261.meta.supp <- ExtractGEOMetaSupp(acce = "GSE292261", supp.idx = 4)
#' }
ExtractGEOMeta <- function(acce, platform = NULL, down.supp = FALSE, supp.idx = 1, timeout = 3600, ...) {
  if (down.supp) {
    message("Download supplementary files to extract sample metadata!")
    geo.meta.df <- ExtractGEOMetaSupp(acce = acce, timeout = timeout, supp.idx = supp.idx)
  } else {
    # get GEO object
    if (is.null(platform)) {
      geo.obj <- GEOquery::getGEO(GEO = acce, ...)
      pfs <- sapply(geo.obj, function(x) {
        Biobase::annotation(x)
      })
      names(geo.obj) <- pfs
    } else {
      geo.obj <- list()
      pf.obj <- GEOobj(acce = acce, platform = platform, ...)
      geo.obj[[platform]] <- pf.obj
    }
    # extract metadata
    geo.meta.list <- lapply(names(geo.obj), function(x) {
      # extract general information
      pf.info <- ExtractGEOInfo(pf.obj = geo.obj[[x]], sample.wise = FALSE)
      pf.info$Platform <- x
      # select meta data
      pf.meta <- ExtractGEOSubMeta(pf.obj = geo.obj[[x]])
      pf.meta$Platform <- x
      # merge all dataframe
      pf.all <- merge(pf.meta, pf.info, by = "Platform", all.x = TRUE) %>% as.data.frame()
      pf.all
    })
    # all metadta
    if (length(geo.meta.list) > 1) {
      geo.meta.df <- data.table::rbindlist(geo.meta.list, fill = TRUE) %>% as.data.frame()
    } else {
      geo.meta.df <- geo.meta.list[[1]]
    }
  }
  return(geo.meta.df)
}

#' Extract Sample Metadata from Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#'
#' @return A dataframe.
#'
ExtractGEOMetaSupp <- function(acce, timeout = 3600, supp.idx = 1) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supplementary file
  supp.down.log <- tryCatch(
    expr = {
      getGEOSuppFilesInner(GEO = acce, baseDir = tmp.folder, index = supp.idx)
    },
    error = function(e) {
      print(e)
      stop("Please check the supp.idx or change the timeout with a larger value.")
    }
  )
  supp.file.path <- rownames(supp.down.log)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    Gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  # read metatada
  if (file.ext %in% c("xlsx", "xls")) {
    # read excel file
    sample.meta <- tryCatch(
      {
        readxl::read_excel(path = supp.file.path, col_names = TRUE) %>% as.data.frame()
      },
      error = function(e) {
        message(e)
        # empty dataframe
        data.frame()
      }
    )
  } else if (file.ext %in% c("csv", "tsv", "txt", "tab")) {
    # read text file
    sample.meta <- tryCatch(
      {
        data.table::fread(file = supp.file.path, header = T) %>% as.data.frame()
      },
      error = function(e) {
        message(e)
        # empty dataframe
        data.frame()
      }
    )
  }
  return(sample.meta)
}

# connect to GEO, extract GEO object, extract platform object
GEOobj <- function(acce, platform, ...) {
  # obtain GEO object
  geo.obj <- GEOquery::getGEO(GEO = acce, ...)

  # extract platform
  pfs <- sapply(geo.obj, function(x) {
    Biobase::annotation(x)
  })
  if (!platform %in% pfs) {
    stop(paste("The platform you provides is not valid!", paste(pfs, collapse = ", ")))
  }
  # extract platform data
  pf.idx <- which(pfs == platform[1])
  pf.obj <- geo.obj[[pf.idx]]
  return(pf.obj)
}

# merge cols into one
SimpleCol <- function(df, col) {
  cols <- grep(pattern = col, x = colnames(df), value = T)
  df.cols <- df[cols]
  value <- paste(c(t(df.cols)), collapse = ". ")
  return(value)
}

#' Extract GEO Study Information.
#'
#' @param pf.obj GEO object of platform.
#' @param sample.wise Logical value, whether to extract sample-wise information. Default: FALSE.
#'
#' @return A dataframe.
#'
ExtractGEOInfo <- function(pf.obj, sample.wise = FALSE) {
  # platform information
  pf.info <- Biobase::experimentData(pf.obj)
  # additional information
  pf.info.add.df <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  if (sample.wise) {
    pf.info.final <- cbind(data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds), pf.info.add.df
    )) %>%
      as.data.frame()
  } else {
    used.cols <- c("organism", "molecule", "strategy", "extract_protocol", "data_processing")
    pf.info.add.used <- pf.info.add.df[, grep(pattern = paste(used.cols, collapse = "|"), colnames(pf.info.add.df), value = T)]
    pf.info.add.sim <- apply(pf.info.add.used, 2, function(x) {
      paste(unique(x), collapse = ", ")
    }) %>%
      t() %>%
      as.data.frame()
    # final information
    pf.info.final <- data.frame(
      Title = pf.info@title,
      Type = Biobase::notes(pf.info)$type,
      Organism = SimpleCol(df = pf.info.add.sim, col = "organism"),
      Abstract = pf.info@abstract,
      Design = pf.info@other$overall_design,
      SampleCount = length(Biobase::sampleNames(pf.obj)),
      Molecule = SimpleCol(df = pf.info.add.sim, col = "molecule"),
      ExtractProtocol = SimpleCol(df = pf.info.add.sim, col = "extract_protocol"),
      LibraryStrategy = SimpleCol(df = pf.info.add.sim, col = "strategy"),
      DataProcessing = SimpleCol(df = pf.info.add.sim, col = "data_processing"),
      SupplementaryFile = gsub(pattern = "\n", replacement = ", ", pf.info@other$supplementary_file),
      Contact = paste(pf.info@name, pf.info@contact, sep = "; "),
      PMID = gsub(pattern = "\n", replacement = ", ", pf.info@pubMedIds)
    )
  }
  return(pf.info.final)
}

#' Extract Sample Metadata.
#'
#' @param pf.obj GEO object of platform.
#'
#' @return A dataframe.
#'
ExtractGEOSubMeta <- function(pf.obj) {
  # extract sample detail information
  pf.info <- as.data.frame(Biobase::pData(Biobase::phenoData(pf.obj)))
  # select used basic cols
  valid.cols <- intersect(colnames(pf.info), c(c("title", "geo_accession", "source_name_ch1", "description")))
  pf.info.used <- pf.info[valid.cols]
  # process characteristics
  pf.info.charac <- pf.info[grep(pattern = "^characteristics", x = colnames(pf.info))]
  ## modify colnames
  pf.info.charac.colnames <-
    apply(pf.info.charac, 2, function(x) {
      x.head <- gsub(pattern = "(.*?): (.*)", replacement = "\\1", x = x)
      paste(unique(x.head), collapse = ",")
    })
  colnames(pf.info.charac) <- pf.info.charac.colnames
  ## modify values
  pf.info.charac <- apply(pf.info.charac, 2, function(x) {
    gsub(pattern = "(.*?): (.*)", replacement = "\\2", x = x)
  })
  ## final meta
  pf.meta <- cbind(pf.info.used, pf.info.charac) %>% as.data.frame()

  return(pf.meta)
}


#' Extract Raw Count Matrix from Supplementary Files.
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param file.regex The regex to extract correct count matrix files. Default: NULL.
#' @param extra.cols Extra columns to remove, e.g., "Chr", "Start", "End", "Strand", "Length" (featureCounts).
#' @param transpose Logical value, whether to transpose the matrix. Used when the number of rows is less than the number of columns.
#'
#' @return A dataframe.
#'
ExtractGEOExpSupp <- function(acce, timeout = 3600, supp.idx = 1, file.regex = NULL,
                              extra.cols = c(
                                "chr", "start", "end", "strand", "length",
                                "width", "chromosome", "seqnames", "seqname",
                                "chrom", "chromosome_name", "seqid", "stop"
                              ),
                              transpose = TRUE) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supplementary file
  supp.down.log <- tryCatch(
    expr = {
      getGEOSuppFilesInner(GEO = acce, baseDir = tmp.folder, index = supp.idx)
    },
    error = function(e) {
      print(e)
      stop("Please check the supp.idx or change the timeout with a larger value.")
    }
  )
  supp.file.path <- rownames(supp.down.log)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    Gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    # detect file extension
    all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE)
    all.file.exts <- sapply(all.files, getExtWithoutCompression)
    if (length(unique(all.file.exts)) > 1) {
      message(
        "Detect more than one file extension: ", paste(unique(all.file.exts), collapse = ", "),
        ".", "\n", "Please make sure you have set parameter file.regex to extract the correct file."
      )
      if (is.null(file.regex)) {
        message("The file.regex is NULL, use the first file extension: ", unique(all.file.exts)[1])
        file.regex <- unique(all.file.exts)[1]
        # subset files
        all.files <- all.file.exts[all.file.exts == file.regex] %>% names()
      } else {
        # subset files
        all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = file.regex)
      }
    } else {
      all.files <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = file.regex)
    }
    # unzip
    unzip.log <- sapply(
      grep(pattern = "gz$", x = all.files, value = T),
      function(x) {
        Gunzip(x, overwrite = TRUE)
      }
    )
    all.files <- gsub(pattern = ".gz", replacement = "", x = all.files)
    # read files
    count.list <- lapply(
      all.files,
      function(x) {
        sample.count <- ReadFile(file.path = x, extra.cols = extra.cols, transpose = transpose) %>% tibble::rownames_to_column(var = "GeneName")
        sample.name <- gsub(pattern = "(GSM[0-9]*).*", replacement = "\\1", x = basename(x))
        if (ncol(sample.count) > 2) {
          message(
            "The count matrix has multiple columns, which may be output by STAR!", "\n",
            "Adding colnames with col1, col2, col3..., you can use dplyr::select(dplyr::ends_with('col2')) to extract the correct count matrix (use 'col2' as example)."
          )
          count.col.names <- paste(sample.name, paste0("col", 1:(ncol(sample.count) - 1)), sep = "_")
        } else {
          count.col.names <- sample.name
        }
        colnames(sample.count) <- c("GeneName", count.col.names)
        sample.count
      }
    )
    # create count matrix
    count.mat <- Reduce(f = function(x, y) {
      merge.mat <- merge(x, y, by = "GeneName", all = T)
    }, x = count.list)
    rownames(count.mat) <- count.mat$GeneName
    count.mat$GeneName <- NULL
  } else if (file.ext %in% c("xlsx", "xls", "csv", "tsv", "txt", "tab")) {
    count.mat <- ReadFile(file.path = supp.file.path, extra.cols = extra.cols, transpose = transpose)
  } else {
    stop("Unsupported file extension: ", file.ext)
  }
  return(count.mat)
}

#' Fortmat Supplementary Files to 10x.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param accept.fmt Vector of accepted 10x output format, MEX (barcode/feature/gene/matrix), h5 (h5/h5.gz). Default: c("MEX", "h5").
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}.
#' Default: TURE.
#'
#' @return NULL
#'
ExtractGEOExpSupp10x <- function(acce, supp.idx = 1, timeout = 3600,
                                 out.folder = NULL, accept.fmt = c("MEX", "h5"), gene2feature = TRUE) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supp file
  supp.down.log <- tryCatch(
    expr = {
      getGEOSuppFilesInner(GEO = acce, baseDir = tmp.folder, index = supp.idx)
    },
    error = function(e) {
      print(e)
      stop("Please check the supp.idx or change the timeout with a larger value.")
    }
  )
  supp.file.path <- rownames(supp.down.log)
  # file unzip
  file.ext <- tools::file_ext(supp.file.path)
  if (file.ext == "gz") {
    # gunzip file
    Gunzip(supp.file.path, overwrite = TRUE)
    supp.file.path <- gsub(pattern = "\\.gz", replacement = "", x = supp.file.path)
    # update file extension
    file.ext <- tools::file_ext(supp.file.path)
  }
  if (file.ext == "tar") {
    # untar
    utils::untar(supp.file.path, exdir = file.path(tmp.folder, acce, "sample"))
    process.log <- Process10xFiles(acce = acce, folder = file.path(tmp.folder, acce, "sample"), accept.fmt = accept.fmt, out.folder = out.folder, gene2feature = gene2feature)
  } else {
    stop("Does not support non-tar(.gz) file for 10x mode!")
  }
}

#' Fortmat Supplementary Files to 10x (separate files).
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param accept.fmt Vector of accepted 10x output format, MEX (barcode/feature/gene/matrix), h5 (h5/h5.gz). Default: c("MEX", "h5").
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}.
#' Default: TURE.
#'
#' @return NULL
#'
ExtractGEOExpSupp10xSingle <- function(acce, timeout = 3600, out.folder = NULL, accept.fmt = c("MEX", "h5"), gene2feature = TRUE) {
  # create tmp folder
  tmp.folder <- tempdir()
  # get current timeout
  if (!is.null(timeout)) {
    message("Change Timeout to: ", timeout)
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
  }
  # download supp file
  supp.down.log <- tryCatch(
    expr = {
      GEOquery::getGEOSuppFiles(GEO = acce, baseDir = tmp.folder)
    },
    error = function(e) {
      print(e)
      stop("You can change the timeout with a larger value.")
    }
  )
  process.log <- Process10xFiles(acce = acce, folder = file.path(tmp.folder, acce), accept.fmt = accept.fmt, out.folder = out.folder, gene2feature = gene2feature)
  return(process.log)
}

#' Extract Raw Count Matrix from Supplementary Files or Fortmat Supplementary Files to 10x.
#'
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param file.regex The regex to extract correct count matrix files. Default: NULL.
#' @param extra.cols Extra columns to remove, e.g., "Chr", "Start", "End", "Strand", "Length" (featureCounts).
#' Default: "chr", "start", "end", "strand", "length", "width", "chromosome", "seqnames", "seqname", "chrom", "chromosome_name", "seqid", "stop".
#' @param transpose Logical value, whether to transpose the matrix. Used when the number of rows is less than the number of columns. Default: TRUE.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param accept.fmt Vector of accepted 10x output format, MEX (barcode/feature/gene/matrix), h5 (h5/h5.gz). Default: c("MEX", "h5").
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#' Default: TURE.
#'
#' @return Count matrix (\code{supp.type} is count) or NULL (\code{supp.type} is 10x).
#'
ExtractGEOExpSuppAll <- function(acce, supp.idx = 1, timeout = 3600,
                                 supp.type = c("count", "10x", "10xSingle"), file.regex = NULL,
                                 extra.cols = c(
                                   "chr", "start", "end", "strand", "length",
                                   "width", "chromosome", "seqnames", "seqname",
                                   "chrom", "chromosome_name", "seqid", "stop"
                                 ),
                                 transpose = TRUE, out.folder = NULL, accept.fmt = c("MEX", "h5"), gene2feature = TRUE) {
  if (supp.type == "count") {
    count.mat <- ExtractGEOExpSupp(
      acce = acce, supp.idx = supp.idx, timeout = timeout, file.regex = file.regex,
      extra.cols = extra.cols, transpose = transpose
    )
    return(count.mat)
  } else if (supp.type == "10x") {
    ExtractGEOExpSupp10x(acce = acce, supp.idx = supp.idx, timeout = timeout, out.folder = out.folder, accept.fmt = accept.fmt, gene2feature = gene2feature)
    return(NULL)
  } else if (supp.type == "10xSingle") {
    ExtractGEOExpSupp10xSingle(acce = acce, timeout = timeout, out.folder = out.folder, accept.fmt = accept.fmt, gene2feature = gene2feature)
    return(NULL)
  }
}

#' Extract Raw Count Matrix or Fortmat Supplementary Files to 10x.
#'
#' @param pf.obj GEO object of platform.
#' @param acce GEO accession number.
#' @param supp.idx The index of supplementary files to download. Default: 1.
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or emoty,
#' download supplementary files automatically). Default: FALSE.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param file.regex The regex to extract correct count matrix files. Default: NULL.
#' @param extra.cols Extra columns to remove, e.g., "Chr", "Start", "End", "Strand", "Length" (featureCounts).
#' Default: "chr", "start", "end", "strand", "length", "width", "chromosome", "seqnames", "seqname", "chrom", "chromosome_name", "seqid", "stop".
#' @param transpose Logical value, whether to transpose the matrix. Used when the number of rows is less than the number of columns. Default: TRUE.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param accept.fmt Vector of accepted 10x output format, MEX (barcode/feature/gene/matrix), h5 (h5/h5.gz). Default: c("MEX", "h5").
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#'
#' @return Count matrix (\code{supp.type} is count) or NULL (\code{supp.type} is 10x/10xSingle).
#'
ExtractGEOExp <- function(pf.obj, acce, supp.idx = 1, down.supp = FALSE, timeout = 3600,
                          supp.type = c("count", "10x", "10xSingle"), file.regex = NULL,
                          extra.cols = c(
                            "chr", "start", "end", "strand", "length",
                            "width", "chromosome", "seqnames", "seqname",
                            "chrom", "chromosome_name", "seqid", "stop"
                          ),
                          transpose = TRUE, out.folder = NULL, accept.fmt = c("MEX", "h5"), gene2feature = TRUE) {
  # check parameters
  supp.type <- match.arg(arg = supp.type)
  # download supplementary files
  if (down.supp) {
    exp.data <- ExtractGEOExpSuppAll(
      acce = acce, supp.idx = supp.idx, timeout = timeout,
      supp.type = supp.type, file.regex = file.regex,
      extra.cols = extra.cols, transpose = transpose,
      out.folder = out.folder, accept.fmt = accept.fmt, gene2feature = gene2feature
    )
  } else {
    expr.mat <- Biobase::exprs(pf.obj)
    if (nrow(expr.mat) == 0) {
      message("Matrix not available! Downloading supplementary files.")
      exp.data <- ExtractGEOExpSuppAll(
        acce = acce, supp.idx = supp.idx, timeout = timeout,
        supp.type = supp.type, file.regex = file.regex,
        extra.cols = extra.cols, transpose = transpose,
        out.folder = out.folder, accept.fmt = accept.fmt, gene2feature = gene2feature
      )
    } else {
      if (all(expr.mat %% 1 == 0)) {
        exp.data <- expr.mat
      } else {
        message("Matrix contains non-integer values! Downloading supplementary files.")
        exp.data <- ExtractGEOExpSuppAll(
          acce = acce, supp.idx = supp.idx, timeout = timeout,
          supp.type = supp.type, file.regex = file.regex,
          extra.cols = extra.cols, transpose = transpose,
          out.folder = out.folder, accept.fmt = accept.fmt, gene2feature = gene2feature
        )
      }
    }
  }
  return(exp.data)
}
