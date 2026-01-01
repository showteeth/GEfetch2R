#' Load Objects, Count Matrix, and Metadata from RData.
#'
#' @param rdata Path to RData file.
#' @param accept.fmt Vector, the format of objects for loading. Default: c("Seurat", "seurat", "SingleCellExperiment",
#' "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"). "Seurat" for Seurat v3, v4; "seurat" for Seurat v2.
#' @param show.object Logical value, whether to show the class of available objects. Default: TRUE.
#' @param return.obj Logical value, whether to load the available objects in \code{accept.fmt} to global environment. Default: TRUE.
#' @param slot Vector, the type of count matrix to pull. Default: c("counts", "data", "scale.data").
#' 'counts': raw, un-normalized counts, 'data': normalized data, scale.data: z-scored/variance-stabilized data.
#'
#' @return List contains count matrix and metadta or NULL.
#' @importFrom magrittr %>%
#' @importFrom SeuratObject Assays Version GetAssayData
#' @importFrom SingleCellExperiment altExpNames objectVersion altExp colData
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom Biobase exprs pData
#' @importFrom DESeq2 counts
#' @export
#'
#' @examples
#' \dontrun{
#' # download RData from GEO
#' ParseGEOProcessed(acce = "GSE244572", timeout = 360000, supp.idx = 1, out.folder = "/path/to/RData", file.ext = c("rdata", "rds", "h5ad", "loom"))
#' # with six assays: "RNA", "ADT", "nADT", "SCT", "integrated", "IADT"
#' multi.assyas.list <- LoadRData(
#'   rdata = "/path/to/GSE244572_RPE_CITESeq.RData",
#'   accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
#'   slot = "counts"
#' )
#' # only one object in GSE244572_RPE_CITESeq.RData
#' names(multi.assyas.list)
#' # metadata
#' multi.assyas.list$obj$meta.data %>% head()
#' # assay names
#' names(multi.assyas.list$obj$count.mat)
#' # raw count matrix (slot)
#' multi.assyas.list$obj$count.mat$RNA$counts[1:5, 1:5]
#' }
LoadRData <- function(rdata, accept.fmt = c(
                        "Seurat", "seurat", "SingleCellExperiment",
                        "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"
                      ),
                      show.object = TRUE, return.obj = TRUE, slot = c("counts", "data", "scale.data")) {
  # create new environment
  load.env <- new.env()
  load(file = rdata, envir = load.env)
  obj.names <- ls(envir = load.env)
  # get all object class
  obj.class <- sapply(mget(ls(envir = load.env), envir = load.env), function(x) {
    x.class <- paste(class(x), collapse = ", ")
    x.class
  }) %>% as.data.frame()
  colnames(obj.class) <- "Class"
  message("The object classes stored in RData: ", paste(unique(obj.class$Class), collapse = ", "), ".")
  if (show.object) {
    print(obj.class)
  }
  # filter class
  valid.obj <- obj.class[obj.class$Class %in% accept.fmt, , drop = FALSE]
  if (nrow(valid.obj) > 0) {
    message("Detect ", nrow(valid.obj), " object(s) in given format(s): ", paste(accept.fmt, collapse = ", "), ".")
    obj.count.list <- list()
    for (nm in rownames(valid.obj)) {
      nm.class <- valid.obj[nm, "Class"]
      if (return.obj) {
        message("Load object: ", nm, " (", nm.class, ")", " to global environment!")
        assign(nm, get(nm, envir = load.env), envir = globalenv())
      }
      message("Extract count matrix and metadata (if available) from: ", nm, " (", nm.class, ").")
      obj.count.list[[nm]] <- ExtractObject(obj = get(nm, envir = load.env), obj.class = nm.class, slot = slot)
    }
    return(obj.count.list)
  } else {
    message("No valid object in given format(s): ", paste(accept.fmt, collapse = ", "), ". Now we will guess the type!")
    message("The slot parameter does not work here!")
    obj.count.list <- list()
    obj.meta.list <- list()
    for (nm in rownames(obj.class)) {
      nm.class <- obj.class[nm, "Class"]
      if (return.obj) {
        message("Load object: ", nm, " (", nm.class, ")", " to global environment!")
        assign(nm, get(nm, envir = load.env), envir = globalenv())
      }
      nm.obj <- get(nm, envir = load.env)
      if (nm.class %in% c("matrix, array", "data.frame")) {
        nm.obj <- as.data.frame(nm.obj)
        nm.obj.colclass <- sapply(nm.obj, class)
        nm.obj.colclass.stat <- table(nm.obj.colclass)
        if (ncol(nm.obj) >= 3) {
          if (length(nm.obj.colclass.stat) == 1 && ("numeric" == names(nm.obj.colclass.stat) || "integer" == names(nm.obj.colclass.stat))) {
            message(nm, " has ", ncol(nm.obj), " columns and each column is numerical! Most likely a count matrix!")
            obj.count.list[[nm]] <- nm.obj
          } else if (all(nm.obj.colclass[2:length(nm.obj.colclass)] == "numeric") || all(nm.obj.colclass[2:length(nm.obj.colclass)] == "integer")) {
            message(nm, " has ", ncol(nm.obj), " columns, all of which are numerical except for the first column! Maybe a count matrix!")
            obj.count.list[[nm]] <- nm.obj
          } else {
            if (any(grepl(pattern = "sample|name|cell|id|library|well|barcode|index|type|condition|treat|group", x = colnames(nm.obj), ignore.case = TRUE))) {
              possible.meta.key <- grep(pattern = "sample|name|cell|id|library|well|barcode|index|type|condition|treat|group", x = colnames(nm.obj), ignore.case = TRUE, value = TRUE)
              message("Detect possible sample metadata keys: ", paste(possible.meta.key, collapse = ", "), " in ", nm, ". Maybe metadata/annotation!")
              obj.meta.list[[nm]] <- nm.obj
            } else {
              message("Can not determine if ", nm, " is a count matrix or metadata/annotation. Load to the global environment, please manually check!")
              print(str(nm.obj))
              if (!return.obj) {
                assign(nm, get(nm, envir = load.env), envir = globalenv())
              }
            }
          }
        } else if (any(grepl(pattern = "sample|name|cell|id|library|well|barcode|index|type|condition|treat|group", x = colnames(nm.obj), ignore.case = TRUE))) {
          possible.meta.key <- grep(pattern = "sample|name|cell|id|library|well|barcode|index|type|condition|treat|group", x = colnames(nm.obj), ignore.case = TRUE, value = TRUE)
          message("Detect possible sample metadata keys: ", paste(possible.meta.key, collapse = ", "), ". Maybe metadata/annotation!")
          obj.meta.list[[nm]] <- nm.obj
        } else {
          message("Can not determine if ", nm, " is metadata/annotation. Load to the global environment, please manually check!")
          print(str(nm.obj))
          if (!return.obj) {
            assign(nm, get(nm, envir = load.env), envir = globalenv())
          }
        }
      } else if (nm.class %in% c("dgCMatrix", "dgRMatrix", "dgTMatrix")) {
        message(nm, " is a sparse matrix. Most likely a count matrix!")
        obj.count.list[[nm]] <- nm.obj
      } else {
        message(nm, " is ", nm.class, ".")
        print(head(nm.obj))
      }
    }
    return(list(count = obj.count.list, meta = obj.meta.list))
  }
}

# Extract count matrix and metadata from object
ExtractObject <- function(obj, obj.class, slot) {
  if (obj.class == "Seurat") {
    # v3, v4 (not test on Seurat v5)
    # get all assays
    all.assays <- SeuratObject::Assays(obj)
    message("Detect Seurat version: ", SeuratObject::Version(obj), ", with assay(s): ", paste(all.assays, collapse = ", "), ".")
    count.mat.list <- lapply(all.assays, function(x) {
      x.list <- list()
      if ("counts" %in% slot) {
        count.mat <- SeuratObject::GetAssayData(obj, assay = x, slot = "counts")
        x.list[["counts"]] <- count.mat
      }
      if ("data" %in% slot) {
        data.mat <- SeuratObject::GetAssayData(obj, assay = x, slot = "data")
        x.list[["data"]] <- data.mat
      }
      if ("scale.data" %in% slot) {
        scale.mat <- SeuratObject::GetAssayData(obj, assay = x, slot = "scale.data")
        x.list[["scale.data"]] <- scale.mat
      }
      return(x.list)
    })
    names(count.mat.list) <- all.assays
    meta.data <- obj@meta.data
  } else if (obj.class == "seurat") {
    # v2
    message("Detect Seurat version: ", obj@version, ".")
    count.mat.list <- list()
    if ("counts" %in% slot) {
      count.mat.list[["counts"]] <- obj@raw.data
    }
    if ("data" %in% slot) {
      count.mat.list[["data"]] <- obj@data
    }
    if ("scale.data" %in% slot) {
      count.mat.list[["scale.data"]] <- obj@scale.data
    }
    meta.data <- obj@meta.data
  } else if (obj.class %in% c("SingleCellExperiment", "cell_data_set")) {
    if (obj.class == "cell_data_set") {
      message("The cell_data_set class is derived from the SingleCellExperiment class!")
    }
    # mainExpName appears in 1.13.4: https://github.com/drisso/SingleCellExperiment/commit/2989ffd4e28bfd8710749afcb8d219cf76e56020
    if (packageVersion(pkg = "SingleCellExperiment") >= "1.13.4") {
      main.exp.name <- SingleCellExperiment::mainExpName(obj)
    } else {
      main.exp.name <- NULL
    }
    if (is.null(main.exp.name)) {
      message("Main experiment name is NULL, use main instead.")
      main.exp.name <- "main"
    }
    all.exps <- c(main.exp.name, SingleCellExperiment::altExpNames(obj))
    message("Detect SingleCellExperiment version: ", SingleCellExperiment::objectVersion(obj), ", with experiment(s): ", paste(all.exps, collapse = ", "), ".")
    all.exps.list <- list(obj)
    if (length(all.exps) > 1) {
      for (altn in SingleCellExperiment::altExpNames(obj)) {
        altn.obj <- SingleCellExperiment::altExp(obj, altn)
        all.exps.list <- c(all.exps.list, altn.obj)
      }
    }
    names(all.exps.list) <- all.exps
    count.mat.list <- lapply(all.exps, function(x) {
      x.obj <- all.exps.list[[x]]
      x.assays <- SummarizedExperiment::assayNames(x.obj)
      message("Detect slot(s): ", paste(x.assays, collapse = ", "), " under experiment: ", x, ".")
      x.list <- list()
      if ("counts" %in% slot) {
        if ("counts" %in% x.assays) {
          x.list[["counts"]] <- SummarizedExperiment::assay(x.obj, "counts")
        } else {
          message("No counts detected under experiment: ", x, ".")
          x.list[["counts"]] <- matrix()
        }
      }
      if ("data" %in% slot) {
        if ("logcounts" %in% x.assays) {
          x.list[["data"]] <- SummarizedExperiment::assay(x.obj, "logcounts")
        } else {
          message("No logcounts detected under experiment: ", x, ".")
          x.list[["data"]] <- matrix()
        }
      }
      if ("scale.data" %in% slot) {
        if ("scaledata" %in% x.assays) {
          x.list[["scale.data"]] <- SummarizedExperiment::assay(x.obj, "scaledata")
        } else if ("scale.data" %in% x.assays) {
          x.list[["scale.data"]] <- SummarizedExperiment::assay(x.obj, "scale.data")
        } else {
          message("No scaledata/scale.data detected under experiment: ", x, ".")
          x.list[["scale.data"]] <- matrix()
        }
      }
      return(x.list)
    })
    names(count.mat.list) <- all.exps
    meta.data <- as.data.frame(SingleCellExperiment::colData(obj))
  } else if (obj.class == "CellDataSet") {
    message("Detect CellDataSet version: ", paste(obj@.__classVersion__$CellDataSet, collapse = "."), ".")
    count.mat <- Biobase::exprs(obj)
    meta.data <- as.data.frame(Biobase::pData(obj))
    count.mat.list <- list()
    if (all(count.mat %% 1 == 0)) {
      if ("counts" %in% slot) {
        count.mat.list[["counts"]] <- count.mat
      }
      if ("data" %in% slot) {
        message("Detect data in slot and integer matrix, there is no normalized data, return empty matrix!")
        count.mat.list[["data"]] <- matrix()
      }
      if ("scale.data" %in% slot) {
        message("Detect scale.data in slot and integer matrix, there is no scaled data, return empty matrix!")
        count.mat.list[["scale.data"]] <- matrix()
      }
    } else if (min(count.mat) < 0) {
      if ("counts" %in% slot) {
        message("Detect counts in slot and matrix containing negative values, there is no raw data, return empty matrix!")
        count.mat.list[["counts"]] <- matrix()
      }
      if ("data" %in% slot) {
        message("Detect data in slot and matrix containing negative values, there is no normalized data, return empty matrix!")
        count.mat.list[["data"]] <- matrix()
      }
      if ("scale.data" %in% slot) {
        count.mat.list[["scale.data"]] <- count.mat
      }
    } else {
      if ("counts" %in% slot) {
        message("Detect counts in slot and matrix containing non-negative and non-integer values, there is no raw data, return empty matrix!")
        count.mat.list[["counts"]] <- matrix()
      }
      if ("data" %in% slot) {
        count.mat.list[["data"]] <- count.mat
      }
      if ("scale.data" %in% slot) {
        message("Detect scale.data in slot and matrix containing non-negative and non-integer values, there is no scaled data, return empty matrix!")
        count.mat.list[["scale.data"]] <- matrix()
      }
    }
  } else if (obj.class == "DESeqDataSet") {
    message("Detect DESeqDataSet version: ", obj@metadata$version, ".")
    count.mat.list <- list()
    if ("counts" %in% slot) {
      count.mat.list[["counts"]] <- DESeq2::counts(obj, normalized = FALSE)
    }
    if ("data" %in% slot) {
      count.mat.list[["data"]] <- tryCatch(
        expr = {
          DESeq2::counts(obj, normalized = TRUE)
        },
        error = function(e) {
          print(e)
          message("Run estimateSizeFactors and get the normalized data!")
          obj <- DESeq2::estimateSizeFactors(obj)
          DESeq2::counts(obj, normalized = TRUE)
        }
      )
    }
    if ("scale.data" %in% slot) {
      message("DESeqDataSet does not contain scaled data!")
      count.mat.list[["scale.data"]] <- matrix()
    }
    meta.data <- as.data.frame(obj@colData)
  } else if (obj.class == "DGEList") {
    message("Detect DGEList object!")
    count.mat.list <- list()
    if ("counts" %in% slot) {
      count.mat.list[["counts"]] <- obj$counts
    }
    if ("data" %in% slot) {
      # require edgeR
      if (!require("edgeR", quietly = TRUE, character.only = TRUE)) {
        stop(
          "Please install edgeR package! You can try BiocManager::install('edgeR)!"
        )
      } else {
        message("Calculate CPM with edgeR::calcNormFactors(method = 'TMM') and edgeR::cpm()!")
      }
      obj <- edgeR::calcNormFactors(obj, method = "TMM")
      count.mat.list[["data"]] <- edgeR::cpm(obj)
    }
    if ("scale.data" %in% slot) {
      message("DGEList does not contain scaled data!")
      count.mat.list[["scale.data"]] <- matrix()
    }
    meta.data <- as.data.frame(obj$samples)
  }
  return(list(count.mat = count.mat.list, meta.data = meta.data))
}
