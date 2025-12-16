#' Download Matrix from GEO and Load to Seurat/DESeq2.
#'
#' @param acce GEO accession number.
#' @param platform Platform information/field. Disable when \code{down.supp} is TRUE. Default: NULL (disable).
#' @param down.supp Logical value, whether to download supplementary files to create count matrix. If TRUE, always
#' download supplementary files. If FALSE, use \code{ExpressionSet} (If contains non-integer or empty,
#' download supplementary files automatically). Default: FALSE.
#' @param supp.idx The index of supplementary files to download. This should be consistent with \code{platform}. Default: 1.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param data.type The data type of the dataset, choose from "sc" (single-cell) and "bulk" (bulk). Default: "sc".
#' @param supp.type The type of downloaded supplementary files, choose from count (count matrix file or single count matrix file),
#' 10x (cellranger output files in tar/gz supplementary files, contains barcodes, genes/features and matrix, e.g. GSE200257)
#' and 10xSingle (cellranger output files in supplementary files directly, e.g. GSE236082). Default: count.
#' @param file.regex The regex to extract correct count matrix files. Default: NULL (when multiple file extensions are available in the downloaded tar file, use the first file extension).
#' @param extra.cols Extra columns to remove, e.g., "Chr", "Start", "End", "Strand", "Length" (featureCounts).
#' Default: "chr", "start", "end", "strand", "length", "width", "chromosome", "seqnames", "seqname", "chrom", "chromosome_name", "seqid", "stop".
#' @param transpose Logical value, whether to transpose the matrix. Used when the number of rows is less than the number of columns. Default: TRUE.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}. Default: TRUE.
#' @param load2R Logical value, whether to load the count matrix to R. Default: TRUE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple 10x files (\code{supp.type} is 10x). Default: FALSE.
#' @param meta.data Dataframe contains sample information for DESeqDataSet, use when \code{data.type} is bulk. Default: NULL.
#' @param fmu Column of \code{meta.data} contains group information. Default: NULL.
#' @param ... Parameters for \code{\link{getGEO}}.
#'
#' @return If \code{load2R} is FALSE, return count matrix. If \code{data.type} is "sc", return Seurat object (if \code{merge} is TRUE) or Seurat object list (if \code{merge} is FALSE).
#' If \code{data.type} is "bulk", return DESeqDataSet.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO getGEOSuppFiles
#' @importFrom xml2 read_html xml_text xml_find_all
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom tools file_ext
#' @importFrom utils untar
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
        all.samples.folder <- dir(out.folder, full.names = TRUE)
        # check file
        valid.samples.folder <- Check10XFiles(folders = all.samples.folder, gene2feature = gene2feature)
        if (length(valid.samples.folder) == 0) {
          stop("No valid sample folder detected under ", out.folder, ". Please check!")
        }
        # load to seurat
        seu.list <- sapply(valid.samples.folder, function(x) {
          x.mat <- Seurat::Read10X(data.dir = x)
          seu.obj <- Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
          seu.obj
        })
        if (isTRUE(merge)) {
          seu.obj <- mergeExperiments(seu.list)
        } else {
          seu.obj <- seu.list
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
    if("MEX" %in% accept.fmt){
      # recognize valid files: barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz and features.tsv.gz
      valid.pat.gz <- "barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$"
      all.files.gz <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat.gz)
    }else{
      all.files.gz = c()
    }
    if("h5" %in% accept.fmt){
      h5.gz.file = list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "h5.gz$")
      if(length(h5.gz.file) > 0){
        # gunzip file
        unzip.log <- sapply(
          h5.gz.file,
          function(x) {
            Gunzip(x, overwrite = TRUE)
          }
        )
      }
      all.files.h5 <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = "h5$")
    }else{
      all.files.h5 = c()
    }
    # prepare out folder
    if (is.null(out.folder)) {
      out.folder <- getwd()
    }
    # rename and move files
    if(length(all.files.gz) > 0){
      message("Detect ", length(all.files.gz), " files in MEX(barcode/feature/gene/matrix) format.")
      # change file name
      if (gene2feature) {
        change.name.log <- sapply(all.files.gz, function(x) {
          if (grepl(pattern = "genes.tsv.gz$", x = x)) {
            new.name <- gsub(pattern = "genes.tsv.gz$", replacement = "features.tsv.gz", x = x)
            file.rename(from = x, to = new.name)
          }
        })
        all.files.gz <- list.files(file.path(tmp.folder, acce, "sample"), full.names = TRUE, pattern = valid.pat.gz)
      }
      # get folder
      all.sample.folder <- sapply(all.files.gz, function(x) {
        # get basename and dirname
        file.name <- basename(x)
        dir.name <- dirname(x)
        # remove file type tag
        file.name <- gsub(pattern = valid.pat.gz, replacement = "", x = file.name)
        # remove possible _ and .
        file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
        file.folder <- file.path(out.folder, file.name)
      })
      # create folder and move file
      move.file.log <- sapply(all.files.gz, function(x) {
        # get folder name
        folder.name <- all.sample.folder[x]
        # create folder
        if (!dir.exists(folder.name)) {
          dir.create(path = folder.name, recursive = TRUE)
        }
        new.file.name <- gsub(pattern = paste0(".*(barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$)"), replacement = "\\1", x = x)
        # move file
        copy.tag <- file.copy(from = x, to = file.path(folder.name, new.file.name))
        # remove the original file
        remove.tag <- file.remove(x)
        copy.tag
      })
    }
    if(length(all.files.h5) > 0){
      message("Detect ", length(all.files.h5), " files in h5 format.")
      # get folder
      all.sample.folder <- sapply(all.files.h5, function(x) {
        # get basename and dirname
        file.name <- basename(x)
        dir.name <- dirname(x)
        # remove file type tag
        file.name <- gsub(pattern = "h5$", replacement = "", x = file.name)
        # remove possible _ and .
        file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
        # remove fix prefix
        file.name <- gsub(pattern = "filtered_feature_bc_matrix$|filtered_gene_bc_matrix$|raw_feature_bc_matrix$|raw_gene_bc_matrix$", replacement = "", x = file.name)
        # remove possible _ and .
        file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
        file.folder <- file.path(out.folder, file.name)
      })
      # create folder and move file
      move.file.log <- sapply(all.files.h5, function(x) {
        # get folder name
        folder.name <- all.sample.folder[x]
        # create folder
        if (!dir.exists(folder.name)) {
          dir.create(path = folder.name, recursive = TRUE)
        }
        # move file
        copy.tag <- file.copy(from = x, to = file.path(folder.name, basename(x)))
        # remove the original file
        remove.tag <- file.remove(x)
        copy.tag
      })
    }
    if(length(all.files.gz) == 0 && length(all.files.h5) == 0){
      stop("No valid 10x format (barcode/feature/gene/matrix, h5) files detected! Please check.")
    }
    message("Process 10x fiels done! All files are in ", out.folder)
  } else {
    stop("Does not support non-tar file for 10x mode!")
  }
}

#' Fortmat Supplementary Files to 10x (separate files).
#'
#' @param acce GEO accession number.
#' @param timeout Timeout for \code{\link{download.file}}. Default: 3600.
#' @param out.folder Output folder to save 10x files. Default: NULL (current working directory).
#' @param gene2feature Logical value, whether to rename \code{genes.tsv.gz} to \code{features.tsv.gz}.
#' Default: TURE.
#'
#' @return NULL
#'
ExtractGEOExpSupp10xSingle <- function(acce, timeout = 3600, out.folder = NULL, gene2feature = TRUE) {
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
  # get valid files
  valid.pat <- "barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$"
  all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat)
  # change file name
  if (gene2feature) {
    change.name.log <- sapply(all.files, function(x) {
      if (grepl(pattern = "genes.tsv.gz$", x = x)) {
        new.name <- gsub(pattern = "genes.tsv.gz$", replacement = "features.tsv.gz", x = x)
        file.rename(from = x, to = new.name)
      }
    })
    all.files <- list.files(file.path(tmp.folder, acce), full.names = TRUE, pattern = valid.pat)
  }
  # prepare out folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # get folder
  all.sample.folder <- sapply(all.files, function(x) {
    # get basename and dirname
    file.name <- basename(x)
    dir.name <- dirname(x)
    # remove file type tag
    file.name <- gsub(pattern = valid.pat, replacement = "", x = file.name)
    # remove possible _ and .
    file.name <- gsub(pattern = "[_.]$", replacement = "", x = file.name)
    file.folder <- file.path(out.folder, file.name)
  })
  # create folder and move file
  move.file.log <- sapply(all.files, function(x) {
    # get folder name
    folder.name <- all.sample.folder[x]
    # create folder
    if (!dir.exists(folder.name)) {
      dir.create(path = folder.name, recursive = TRUE)
    }
    new.file.name <- gsub(pattern = paste0(".*(barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$)"), replacement = "\\1", x = x)
    # move file
    copy.tag <- file.copy(from = x, to = file.path(folder.name, new.file.name))
    # remove the original file
    remove.tag <- file.remove(x)
    copy.tag
  })
  message("Process 10x fiels done! All files are in ", out.folder)
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
    ExtractGEOExpSupp10xSingle(acce = acce, timeout = timeout, out.folder = out.folder, gene2feature = gene2feature)
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
