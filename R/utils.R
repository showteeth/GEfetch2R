# merge Seurat object, modify from: https://github.com/dosorio/rPanglaoDB/blob/master/R/mergeExperiments.R
mergeExperiments <- function(experimentList) {
  el.df <- lapply(experimentList, FUN = function(x) {
    dim(x)
  })
  el.df <- as.data.frame(t(as.data.frame(el.df)))
  rownames(el.df) <- 1:nrow(el.df)
  el.df.vec <- apply(el.df, 1, function(row) all(row != 0))
  empty.seu.index <- names(el.df.vec)[el.df.vec == FALSE]
  if (length(empty.seu.index) > 0) {
    message("Detect empty SeuratObject: ", paste0(empty.seu.index, collapse = ", "), ". Skip these!")
  }
  experimentList <- experimentList[el.df.vec]
  for (i in seq_along(experimentList)[-1]) {
    experimentList[[1]] <- suppressWarnings(merge(experimentList[[1]], experimentList[[i]]))
    experimentList[[i]] <- methods::new("Seurat")
  }
  experimentList <- experimentList[[1]]
  return(experimentList)
}

# check columns existence
CheckColumns <- function(df, columns) {
  if (!all(columns %in% colnames(df))) {
    miss.cols <- setdiff(columns, colnames(df))
    stop(paste0(paste(miss.cols, collapse = ", "), " does not exist, Please Check!"))
  }
}

# used in UCSCCellBrowser and cellxgene, merge multiple attributes
PasteAttr <- function(df, attr) {
  for (at in attr) {
    df[[at]] <- sapply(df[[at]], function(x) {
      paste0(x, collapse = ", ")
    })
  }
  return(df)
}

# used in UCSCCellBrowser, recursively extract samples
ExtractSample <- function(df, base.url, json.folder, quiet) {
  # prepare json
  if (base.url != json.folder) {
    df.json <- file.path(base.url, df$name, "dataset.json")
    names(df.json) <- df$name
    df.json.folder <- file.path(json.folder, df$name)
    names(df.json.folder) <- df$name
    df.desc <- file.path(base.url, df$name, "desc.json")
    names(df.desc) <- df$name
    dird <- sapply(df.json.folder, function(x) {
      dir.create(x, showWarnings = FALSE, recursive = TRUE)
    })
    down.status <- lapply(df$name, function(x) {
      utils::download.file(url = df.json[x], destfile = file.path(df.json.folder[x], "dataset.json"), quiet = quiet, mode = "wb", method = "wget", extra = "--no-check-certificate")
      utils::download.file(url = df.desc[x], destfile = file.path(df.json.folder[x], "desc.json"), quiet = quiet, mode = "wb", method = "wget", extra = "--no-check-certificate")
    })
  }
  # process
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu.json.folder <- file.path(json.folder, cf$name)
    cul <- lapply(file.path(cu.json.folder, "dataset.json"), function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSample(df = cu.df, base.url = base.url, json.folder = json.folder, quiet = quiet)), fill = TRUE))
  }
}

# used in UCSCCellBrowser, recursively extract samples online
ExtractSampleOnline <- function(df) {
  base.url <- "https://cells.ucsc.edu/"
  if (!"isCollection" %in% colnames(df)) {
    return(df)
  } else {
    cf <- df[!is.na(df$isCollection), ]
    sf <- df[is.na(df$isCollection), ]
    cu <- file.path(base.url, cf$name, "dataset.json")
    cul <- lapply(cu, function(x) {
      x.json <- jsonlite::fromJSON(txt = x)
      x.df <- jsonlite::flatten(x.json$datasets)
      colnames(x.df) <- gsub(pattern = ".*\\.", replacement = "", x = colnames(x.df))
      x.df
    })
    cu.df <- data.table::rbindlist(cul, fill = TRUE)
    # df = data.table::rbindlist(list(sf, cu.df), fill = TRUE)
    # return(list(sf, ExtractSample(cu.df)))
    return(data.table::rbindlist(list(sf, ExtractSampleOnline(cu.df)), fill = TRUE))
  }
}

# used in UCSCCellBrowser, inherit attributes from parents
InheritParient <- function(df, attr) {
  for (at in attr) {
    df[[at]] <- ifelse(df[[at]] == "", ifelse(df[[paste0("parent_", at)]] == "", "", paste0(df[[paste0("parent_", at)]], "|parent")), df[[at]])
  }
  return(df)
}

# used in UCSCCellBrowser, extract sample attribute
ExtractDesc <- function(lst, attr) {
  at.list <- list()
  for (atn in names(attr)) {
    at <- attr[atn]
    if (at %in% names(lst)) {
      at.value <- paste0(lst[[at]], collapse = ", ")
    } else {
      at.value <- ""
    }
    at.list[[atn]] <- at.value
  }
  return(as.data.frame(at.list))
}

# used in UCSCCellBrowser, filter attributes and return index
CheckParas <- function(df, column, para.value, fuzzy.match = TRUE) {
  # convert to lower case to avoid case-sensitive
  all.values <- gsub("\\|parent$", replacement = "", x = tolower(df[[column]]))
  if (is.null(para.value)) {
    message("Use all ", column, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    para.value <- tolower(para.value)
    # deal with fuzzy match
    if (fuzzy.match) {
      value.list <- sapply(para.value, function(x) grep(pattern = x, x = all.values, fixed = TRUE))
      value.list.len <- sapply(value.list, function(x) length(x))
      invalid.value <- names(value.list.len[value.list.len == 0])
      value <- unique(unlist(value.list))
    } else {
      if (column %in% c("body_parts", "diseases", "organisms", "projects")) {
        # value contains dot
        value.list <- list()
        for (pv in para.value) {
          pv.vec <- c()
          for (avi in 1:length(all.values)) {
            av <- all.values[avi]
            av.vec <- strsplit(x = av, split = ", ")[[1]]
            if (pv %in% av.vec) {
              pv.vec <- c(pv.vec, avi)
            }
          }
          value.list[[pv]] <- pv.vec
        }
        # get invalid value
        invalid.value <- setdiff(para.value, names(value.list))
        value <- unique(unlist(value.list))
      } else {
        # value doesn't contain dot
        value.list <- lapply(para.value, function(x) {
          which(all.values %in% x)
        })
        names(value.list) <- para.value
        value.list.len <- sapply(value.list, function(x) length(x))
        invalid.value <- names(value.list.len[value.list.len == 0])
        value <- unique(unlist(value.list))
      }
    }
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", column)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(para.value, collapse = ", "), " in ", column, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# used in UCSCCellBrowser, create seurat object (add coord to metadata)
Load2Seurat <- function(exp.file, barcode.url = NULL, feature.url = NULL,
                        meta.file, coord.file = NULL, name = NULL, obs.value.filter = NULL,
                        obs.keys = NULL, include.genes = NULL) {
  # source: https://cellbrowser.readthedocs.io/en/master/load.html
  # read matrix
  if (is.null(barcode.url)) {
    mat <- data.table::fread(exp.file, check.names = FALSE)
    # get genes
    genes <- gsub(".+[|]", "", mat[, 1][[1]])
    # modify mat genes
    mat <- data.frame(mat[, -1], row.names = genes, check.names = FALSE)
  } else {
    # with default parameters
    mat <- Read10XOnline(matrix.url = exp.file, barcode.url = barcode.url, feature.url = feature.url)
  }
  # read metadata
  meta <- data.frame(data.table::fread(meta.file, check.names = FALSE), row.names = 1, check.names = FALSE)
  # filter dataset metadata
  ## filter cell's metadata values
  if (!is.null(obs.value.filter)) {
    # TODO: test key of obs.value.filter in colnames(meta)
    tryCatch(
      {
        raw.meta <- meta
        meta <- meta %>% dplyr::filter(eval(rlang::parse_expr(obs.value.filter)))
        if (nrow(meta) == 0) {
          message("Please check the value of obs.value.filter, e.g. the 'oligodendrocyte' in cell_type == 'oligodendrocyte'!")
          meta <- raw.meta
        }
      },
      error = function(cond) {
        message("Please check obs.value.filter: ", cond)
      }
    )
  }
  ## filter cell's metadata colnames
  if (!is.null(obs.keys)) {
    valid.obs.keys <- intersect(colnames(meta), obs.keys)
    invalid.obs.keys <- setdiff(obs.keys, colnames(meta))
    if (length(invalid.obs.keys) > 0) {
      message("Detected invalid obs.keys: ", paste(invalid.obs.keys, collapse = ", "))
    }
    if (length(valid.obs.keys) > 0) {
      meta <- meta[valid.obs.keys]
    } else {
      message("No valid obs.keys, return all obs.keys!")
      meta <- meta
    }
  }
  ## filter genes
  if (!is.null(include.genes)) {
    mat <- mat[include.genes, rownames(meta)]
  } else {
    mat <- mat[, rownames(meta)]
  }
  if (is.null(coord.file)) {
    seu.obj <- Seurat::CreateSeuratObject(counts = mat, project = name, meta.data = meta)
  } else {
    # prepare coord file
    coord.list <- lapply(1:length(coord.file), function(x) {
      coord.name <- gsub(pattern = ".coords.tsv.gz", replacement = "", x = basename(coord.file[x]))
      coord.df <- data.frame(data.table::fread(coord.file[x], check.names = FALSE), row.names = 1, check.names = FALSE)
      colnames(coord.df) <- paste(coord.name, 1:ncol(coord.df), sep = "_")
      coord.df$Barcode <- rownames(coord.df)
      return(coord.df)
    })
    if (length(coord.file) == 1) {
      all.coord.df <- coord.list[[1]]
    } else {
      all.coord.df <- coord.list %>% purrr::reduce(dplyr::full_join, by = "Barcode")
    }
    # merge metadata
    meta <- merge(meta, all.coord.df, by.x = 0, by.y = "Barcode", all.x = TRUE) %>%
      tibble::column_to_rownames(var = "Row.names")
    seu.obj <- Seurat::CreateSeuratObject(counts = mat, project = name, meta.data = meta)
  }
  return(seu.obj)
}

Loading2DESeq2 <- function(mat, meta, fmu) {
  # loadding into DESeq2
  if (is.null(meta)) {
    meta <- data.frame(condition = colnames(mat))
    rownames(meta) <- colnames(mat)
    meta$condition <- as.factor(meta$condition)
    fmu <- "condition"
  } else {
    if (all(rownames(meta) != colnames(mat))) {
      stop("The columns of the count matrix and the rows of the meta.data are not in the same order!")
    }
  }
  if (is.null(fmu)) {
    message("The condition column (fmu) is empty, use the first column!")
    fmu <- colnames(meta)[1]
  }
  fmu.used <- stats::formula(paste("~", fmu))
  de.obj <- DESeq2::DESeqDataSetFromMatrix(
    countData = mat,
    colData = meta, design = fmu.used
  )
  return(de.obj)
}

# used in UCSCCellBrowser
# source: https://github.com/satijalab/seurat/blob/master/R/utilities.R#L1949
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

# used in UCSCCellBrowser, load cellranger output to matrix
Read10XOnline <- function(matrix.url, barcode.url, feature.url, gene.column = 2,
                          cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
  # load matrix
  data <- Matrix::readMM(file = gzcon(url(matrix.url)))
  # load barcode
  cell.barcodes <- as.data.frame(data.table::fread(barcode.url, header = FALSE))
  cn <- ifelse(ncol(x = cell.barcodes) > 1, cell.column, 1)
  cell.names <- cell.barcodes[, cn]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # matrix colnames
  colnames(x = data) <- cell.names
  # load feature
  feature.names <- as.data.frame(data.table::fread(feature.url, header = FALSE))
  # modify gene column
  gene.column <- min(ncol(feature.names), gene.column)
  if (any(is.na(x = feature.names[, gene.column]))) {
    warning(
      "Some features names are NA. Replacing NA names with ID from the opposite column requested",
      call. = FALSE,
      immediate. = TRUE
    )
    na.features <- which(x = is.na(x = feature.names[, gene.column]))
    replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
    feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
  }

  # modify matrix rownames
  if (unique.features) {
    rownames(x = data) <- make.unique(names = feature.names[, gene.column])
  } else {
    rownames(x = data) <- feature.names[, gene.column]
  }
  # In cell ranger 3.0, a third column specifying the type of data was added
  # and we will return each type of data as a separate matrix
  if (ncol(x = feature.names) > 2) {
    data_types <- factor(x = feature.names$V3)
    lvls <- levels(x = data_types)
    if (length(x = lvls) > 1) {
      message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
    }
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) { # Return Gene Expression first
      lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
    }
    data <- lapply(
      X = lvls,
      FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      }
    )
    names(x = data) <- lvls
  } else {
    data <- list(data)
  }
  # convert to dgCMatrix
  final.data <- Seurat::as.sparse(data[[1]])
  return(final.data)
}

# used in UCSCCellBrowser, extract desc.json
ParseCBdesc <- function(desc.files) {
  desc.list <- lapply(names(desc.files), function(x) {
    sd.json <- jsonlite::fromJSON(txt = desc.files[x])
    used.attr <- c(
      "title" = "title", "paper" = "paper_url", "abstract" = "abstract", "unit" = "unitDesc",
      "coords" = "coordFiles", "methods" = "methods", "geo" = "geo_series"
    )
    sd.df <- ExtractDesc(lst = sd.json, attr = used.attr)
    sd.df$name <- x
    sd.df
  })
  desc.df <- do.call(rbind, desc.list)
  return(desc.df)
}

# used in UCSCCellBrowser, extract desc.json
ParseCBdataset <- function(datasets.files) {
  datasets.list <- lapply(names(datasets.files), function(x) {
    x.json <- jsonlite::fromJSON(txt = datasets.files[x])
    x.file.df <- jsonlite::flatten(as.data.frame(x.json$fileVersions))
    mat.name <- ifelse("outMatrix.fname" %in% colnames(x.file.df), basename(x.file.df[, "outMatrix.fname"]), "")
    barcode.name <- ifelse("barcodes.fname" %in% colnames(x.file.df), basename(x.file.df[, "barcodes.fname"]), "")
    feature.name <- ifelse("features.fname" %in% colnames(x.file.df), basename(x.file.df[, "features.fname"]), "")
    x.file.df <- data.frame(matrix = mat.name, barcode = barcode.name, feature = feature.name)
    x.file.df$name <- x
    x.file.df
  })
  datasets.df <- do.call(rbind, datasets.list)
  return(datasets.df)
}

# used in UCSCCellBrowser, get dataset information from link
ParseCBlink <- function(link) {
  dataset.vec <- gsub(pattern = "https://cells.ucsc.edu/?ds=", replacement = "", x = link, fixed = TRUE)
  dataset.vec <- gsub(pattern = "&.*$", replacement = "", x = dataset.vec)
  dataset.name <- gsub(pattern = "+", replacement = "/", x = dataset.vec, fixed = TRUE)
  # resolve collection and dataset
  base.url <- "https://cells.ucsc.edu/"
  dataset.info.list <- lapply(dataset.name, function(x) {
    x.url <- file.path(base.url, x, "dataset.json")
    x.json <- jsonlite::fromJSON(txt = x.url)
    if ("datasets" %in% names(x.json)) {
      x.df <- data.frame(name = x, isCollection = TRUE)
      x.samples.df <- ExtractSampleOnline(x.df) %>% as.data.frame()
      x.samples.df[, c("name", "shortLabel")]
    } else {
      data.frame(name = x.json$name, shortLabel = x.json$shortLabel)
    }
  })
  dataset.info.df <- do.call(rbind, dataset.info.list) %>% dplyr::distinct()
  # add matrix information
  datasets.files <- file.path(base.url, dataset.info.df$name, "dataset.json")
  names(datasets.files) <- dataset.info.df$name
  datasets.df <- ParseCBdataset(datasets.files)
  datasets.df$matrixType <- ifelse(datasets.df$barcode == "", "matrix", "10x")
  # add sample information
  desc.files <- file.path(base.url, datasets.df$name, "desc.json")
  names(desc.files) <- datasets.df$name
  desc.df <- ParseCBdesc(desc.files)
  meta <- merge(datasets.df, desc.df, by = "name") %>% as.data.frame()
  return(meta)
}

# used in GEO, check the integrity of 10x files
Check10XFiles <- function(folders, gene2feature) {
  folders.flag <- sapply(folders, function(x) {
    if (gene2feature) {
      if (file.exists(file.path(x, "matrix.mtx.gz")) &&
        file.exists(file.path(x, "barcodes.tsv.gz")) &&
        file.exists(file.path(x, "features.tsv.gz"))) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      if (file.exists(file.path(x, "matrix.mtx.gz")) &&
        file.exists(file.path(x, "barcodes.tsv.gz"))) {
        if (file.exists(file.path(x, "features.tsv.gz")) ||
          file.exists(file.path(x, "genes.tsv.gz"))) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      } else {
        return(FALSE)
      }
    }
  })
  valid.folders <- folders[folders.flag]
  drop.folders <- setdiff(folders, valid.folders)
  if (length(drop.folders) > 0) {
    if (gene2feature) {
      message(paste0(drop.folders, collapse = ", "), " don't contain matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz.")
    } else {
      message(paste0(drop.folders, collapse = ", "), " don't contain matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz/genes.tsv.gz.")
    }
  }
  return(valid.folders)
}

# check pair-end or single-end
CheckBam <- function(bam, samtools.path) {
  samtools.cmd <- paste(samtools.path, "view -h", bam, "2>/dev/null |head -n 100000 |", samtools.path, "view -c -f 1 -")
  # run command
  message(paste("Check pair/single end: ", samtools.cmd))
  samtools.status <- system(samtools.cmd, intern = TRUE)
  samtools.status <- as.numeric(samtools.status)
  if (samtools.status == 0) {
    return(FALSE)
  } else if (samtools.status > 0) {
    return(TRUE)
  }
}

# used in cellxxgene, extract content from url
URLRetrieval <- function(url) {
  url.page <- curl::curl_fetch_memory(url)
  url.content <- jsonlite::fromJSON(rawToChar(url.page$content))
  return(url.content)
}

# used in cellxgene, merge multiple attributes
PasteAttrCXG <- function(df, attr, col) {
  for (at in attr) {
    df[[at]] <- sapply(df[[at]], function(x) {
      paste0(x[[col]], collapse = ", ")
    })
  }
  return(df)
}

# used in cellxgene, filter datasets
cellxgeneAttrFilter <- function(df, attr, attr.value) {
  if (is.null(attr.value)) {
    message("Use all ", attr, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    # lower case
    attr.value <- tolower(attr.value)
    all.values <- df[[attr]]
    # value contains dot
    value.list <- list()
    for (pv in attr.value) {
      pv.vec <- c()
      for (avi in 1:length(all.values)) {
        av <- all.values[avi]
        av.vec <- strsplit(x = av, split = ", ")[[1]] %>% tolower()
        if (pv %in% av.vec) {
          pv.vec <- c(pv.vec, avi)
        }
      }
      value.list[[pv]] <- pv.vec
    }
    # get invalid value
    invalid.value <- setdiff(attr.value, names(value.list))
    value <- unique(unlist(value.list))
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", attr)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(attr.value, collapse = ", "), " in ", attr, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# used in cellxgene, get download urls
PostDatasetURL <- function(url) {
  response <- httr::POST(url)
  httr::stop_for_status(response)
  result <- httr::content(response, as = "text", encoding = "UTF-8")
  result.list <- jsonlite::fromJSON(result)
  presigned.url <- result.list$presigned_url
  names(presigned.url) <- result.list$file_name
  return(presigned.url)
}

# used in cellxgene, prepare download urls with metadata
# reference: https://gist.github.com/ivirshup/f1a1603db69de3888eacb4bdb6a9317a
PrepareCELLxGENEUrls <- function(df, fe) {
  # file extension column
  fe.id <- paste0(fe, "_id")
  CheckColumns(df = df, columns = c("dataset_id", fe.id, "dataset_description"))
  invalid.df <- df[is.na(df[[fe.id]]) | is.na(df$dataset_id) | df$dataset_id == "" | df[[fe.id]] == "", ]
  # message("Detect ", nrow(invalid.df), " invalid metadata (", fe.id, "/dataset_id is empty or NA).")
  message("There is no file in dataset_id: ", paste(invalid.df$dataset_id, collapse = ", "), " with extension rds.")
  valid.df <- df[!(is.na(df[[fe.id]]) | is.na(df$dataset_id) | df$dataset_id == "" | df[[fe.id]] == ""), ]
  valid.urls <- df[[fe.id]]
  valid.names <- make.names(valid.df$dataset_description, unique = TRUE)
  valid.filenames <- paste0(valid.names, ".", fe)
  names(valid.urls) <- valid.filenames
  return(list(df = valid.df, urls = valid.urls))
}

# used in cellxgene, deal with NULL value
CXGnullCol <- function(value) {
  if (is.null(value)) {
    return(NA)
  } else {
    return(value)
  }
}

# used in cellxgene, parse dataset information of single collection
ParseCXGcollection <- function(collection_id) {
  # base url
  cellxgene.collections.url <- "https://api.cellxgene.cziscience.com/curation/v1/collections/"
  # single collection url
  cellxgene.sc.url <- paste0(cellxgene.collections.url, collection_id)
  cellxgene.sc.content <- URLRetrieval(cellxgene.sc.url)
  # check collection_id
  if (length(names(cellxgene.sc.content)) == 4 && all(names(cellxgene.sc.content) == c("detail", "status", "title", "type"))) {
    print(as.data.frame(cellxgene.sc.content))
    stop("Error occurred when accessing collection: ", collection_id, ". If CheckAPI() returns OK, please check collection_id (link) you provided!")
  }
  cellxgene.sc.datasets <- jsonlite::flatten(cellxgene.sc.content$datasets)
  colnames(cellxgene.sc.datasets) <- gsub(pattern = "^title$", replacement = "dataset_description", x = colnames(cellxgene.sc.datasets))
  # prepare metadata
  cellxgene.sc.df <- data.frame(
    title = cellxgene.sc.content$name, description = cellxgene.sc.content$description, doi = CXGnullCol(cellxgene.sc.content$doi),
    contact = CXGnullCol(cellxgene.sc.content$contact_name), contact_email = CXGnullCol(cellxgene.sc.content$contact_email),
    collection_id = cellxgene.sc.content$collection_id, collection_url = cellxgene.sc.content$collection_url,
    consortia = paste(unlist(cellxgene.sc.content$consortia), collapse = ", "), curator_name = CXGnullCol(cellxgene.sc.content$curator_name),
    visibility = cellxgene.sc.content$visibility
  )
  cellxgene.sc.df <- cbind(cellxgene.sc.df, cellxgene.sc.datasets) %>% as.data.frame()
  # modify metadata
  label.col <- c(
    "assay", "cell_type", "organism", "self_reported_ethnicity", "sex", "tissue",
    "disease", "development_stage"
  )
  valid.label.col <- intersect(colnames(cellxgene.sc.df), label.col)
  cellxgene.sc.df <-
    PasteAttrCXG(
      df = cellxgene.sc.df,
      attr = valid.label.col, col = "label"
    )
  list.col <- c("batch_condition", "suspension_type", "donor_id")
  valid.list.col <- intersect(colnames(cellxgene.sc.df), list.col)
  cellxgene.sc.df <-
    PasteAttr(df = cellxgene.sc.df, attr = valid.list.col)
  # add h5ad and rds information
  cellxgene.sc.list <- lapply(1:nrow(cellxgene.sc.df), function(x) {
    x.df <- cellxgene.sc.df[x, ]
    x.df.dataset <- x.df$assets[[1]]
    # x.df$dataset_id <- unique(x.df.dataset$dataset_id)
    if ("RDS" %in% unique(x.df.dataset$filetype)) {
      x.rds.idx <- which(x.df.dataset$filetype == "RDS")
      x.df$rds_id <- x.df.dataset$url[x.rds.idx]
    } else {
      x.df$rds_id <- NA
    }
    if ("H5AD" %in% unique(x.df.dataset$filetype)) {
      x.h5ad.idx <- which(x.df.dataset$filetype == "H5AD")
      x.df$h5ad_id <- x.df.dataset$url[x.h5ad.idx]
    } else {
      x.df$h5ad_id <- NA
    }
    x.df
  })
  cellxgene.sc.final <- data.table::rbindlist(cellxgene.sc.list, fill = TRUE) %>%
    as.data.frame()
  return(cellxgene.sc.final)
}

# used in cellxgene, parse dataset information
ParseCXGdataset <- function(dataset_id) {
  dataset.url <- paste0("https://api.cellxgene.cziscience.com/curation/v1/datasets/", dataset_id, "/versions")
  dataset.content <- URLRetrieval(dataset.url)
  # check dataset_id
  if (length(names(dataset.content)) == 4 && all(names(dataset.content) == c("detail", "status", "title", "type"))) {
    print(as.data.frame(dataset.content))
    stop("Error occurred when accessing collection: ", dataset_id, ". If CheckAPI() returns OK, please check dataset_id (link) you provided!")
  }
  # get first dataset, newer
  collection.id <- dataset.content[1, "collection_id"]
  collection.info <- ParseCXGcollection(collection_id = collection.id)
  used.collection.info <- collection.info[collection.info$dataset_id == dataset_id, ]
  return(used.collection.info)
}

# used in hca, recursively extract projects
# reference: https://bioconductor.org/packages/release/bioc/html/hca.html
RecurURLRetrieval <- function(url) {
  url.content <- URLRetrieval(url)
  next.url <- url.content$pagination$`next`
  if (!is.null(next.url)) {
    # return(c(url.content, RecurURLRetrieval(next.url)))
    return(data.table::rbindlist(list(url.content$hits, RecurURLRetrieval(next.url)), fill = TRUE))
  } else {
    return(url.content$hits)
  }
}

# used in hca, two-level list, final is vector
HCAPasteCol <- function(df, col) {
  if (col %in% colnames(df)) {
    col.value <- paste0(sapply(
      df[[col]],
      function(x) {
        ifelse(is.null(x), "",
          paste0(x, collapse = "|")
        )
      }
    ), collapse = ", ")
  } else {
    col.value <- ""
  }
  return(col.value)
}

# used in hca, dataframe, check column exists
HCAPasteColdf <- function(df, col = NULL) {
  if (col %in% colnames(df)) {
    return(paste0(df[[col]], collapse = ", "))
  } else {
    return("")
  }
}

# used in hca, filter proojects
HCAAttrFilter <- function(df, attr, attr.value) {
  if (is.null(attr.value)) {
    message("Use all ", attr, " as input!")
    # return row index
    value <- 1:nrow(df)
  } else {
    # lower case
    attr.value <- tolower(attr.value)
    all.values <- df[[attr]]
    # value contains dot
    value.list <- list()
    for (pv in attr.value) {
      pv.vec <- c()
      for (avi in 1:length(all.values)) {
        av <- all.values[avi]
        av.vec <- strsplit(x = strsplit(x = av, split = ", ")[[1]], split = "\\|") %>%
          unlist() %>%
          tolower()
        if (pv %in% av.vec) {
          pv.vec <- c(pv.vec, avi)
        }
      }
      value.list[[pv]] <- pv.vec
    }
    # get invalid value
    invalid.value <- setdiff(attr.value, names(value.list))
    value <- unique(unlist(value.list))
    # print invalid value
    if (length(invalid.value) > 0) {
      message(paste0(invalid.value, collapse = ", "), " are not valid for ", attr)
    }
    # deal with empty value
    if (length(value) == 0) {
      warning("There is no valid value under ", paste(attr.value, collapse = ", "), " in ", attr, ". This filter is invalid!")
      value <- 1:nrow(df)
    }
  }
  return(value)
}

# get filter value for "PanglaoDB", "UCSC", "CELLxGENE", "HCA"
CheckFilter <- function(df, filter, all.filter, database, combine) {
  if (combine) {
    # extract dataframe
    filter.df <- df[, all.filter[filter]]
    # filter all empty row
    filter.df <- dplyr::filter(filter.df, !dplyr::if_all(
      dplyr::everything(),
      function(x) {
        x == "" | is.na(x)
      }
    ))
    # for UCSC, remove parent
    if (database == "UCSC") {
      filter.df <- apply(filter.df, 2, FUN = function(x) {
        gsub("\\|parent$", replacement = "", x = tolower(x))
      }) %>% as.data.frame()
    }
    # split dataframe, a dataset may contain multiple value, eg: multiple organ
    ## sep is ,
    for (col in colnames(filter.df)) {
      filter.df <- tidyr::separate_rows(filter.df, tidyr::all_of(col), sep = ", ")
    }
    ## sep is |
    for (col in colnames(filter.df)) {
      filter.df <- tidyr::separate_rows(filter.df, tidyr::all_of(col), sep = "\\|")
    }
    if (database != "PanglaoDB") {
      filter.df <- apply(filter.df, 2, tolower) %>% as.data.frame()
    }
    # summarise
    filter.df.stat <- filter.df %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(Num = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(Num))
    return(filter.df.stat)
  } else {
    filter.list <- lapply(filter, function(x) {
      filter.values <- df[[all.filter[x]]]
      if (database == "UCSC") {
        filter.values <- gsub("\\|parent$", replacement = "", x = tolower(filter.values))
      }
      if (database == "PanglaoDB") {
        vf.df <- strsplit(x = unlist(strsplit(x = filter.values, split = ", ")), split = "\\|") %>%
          unlist() %>%
          table() %>%
          as.data.frame()
      } else {
        vf.df <- strsplit(x = unlist(strsplit(x = filter.values, split = ", ")), split = "\\|") %>%
          unlist() %>%
          tolower() %>%
          table() %>%
          as.data.frame()
      }
      colnames(vf.df) <- c("Value", "Num")
      vf.df <- vf.df[order(vf.df$Num, decreasing = TRUE), ]
      vf.df$Key <- x
      rownames(vf.df) <- NULL
      return(vf.df)
    })
    names(filter.list) <- filter
    return(filter.list)
  }
}

# used in hca, extract data information from contributedAnalyses and matrices
HCAExtactData <- function(df) {
  # unlist
  df.vec <- unlist(df)
  # create dataframe
  df.unlist <- data.frame(meta = names(df.vec), value = df.vec)
  rownames(df.unlist) <- NULL
  # data columns
  data.cols <- c(
    "contentDescription", "format", "isIntermediate", "name", "sha256", "size",
    "fileSource", "uuid", "version", "matrixCellCount", "drs_uri", "url"
  )
  data.col.pattern <- paste0(data.cols, collapse = "|")
  type.pattern <- paste0("(.*)\\.(", data.col.pattern, ")([0-9]*)")
  # add col
  df.unlist$type <- gsub(pattern = type.pattern, replacement = "\\2", x = df.unlist$meta)
  df.unlist$num <- gsub(pattern = type.pattern, replacement = "\\3", x = df.unlist$meta)
  df.unlist$meta <- gsub(pattern = type.pattern, replacement = "\\1", x = df.unlist$meta)
  df.unlist$meta <- paste0(df.unlist$meta, ".", df.unlist$num)
  df.final <- tidyr::spread(data = df.unlist[c("meta", "type", "value")], key = "type", value = "value")
  return(df.final)
}

# used in hca, parse dataset information of project
ParseHCAdataset <- function(project.id) {
  # get all catalogs
  catalog.vec <- ExtractHCACatalogs()
  # get project detail
  projects.cat.list <- lapply(catalog.vec, function(x) {
    project.url <- paste0(
      "https://service.azul.data.humancellatlas.org/index/projects?catalog=", x,
      "&filters=%7B%0A%20%20%22projectId%22%3A%20%7B%0A%20%20%20%20%22is%22%3A%20%5B%0A%20%20%20%20%20%20%22",
      project.id, "%22%0A%20%20%20%20%5D%0A%20%20%7D%0A%7D"
    )
    projects.list <- URLRetrieval(project.url)
    projects.info <- projects.list$hits %>% as.data.frame()
    if (nrow(projects.info) != 0) {
      projects.info$catalog <- x
    }
    return(projects.info)
  })
  projects.cat.df <- data.table::rbindlist(projects.cat.list, fill = TRUE) %>% as.data.frame()
  if (nrow(projects.cat.df) == 0) {
    stop("Error occurred when accessing project: ", project.id, ". If CheckAPI() returns OK, please check project.id (link) you provided!")
  }
  # remove duplicated information in different catalogs
  projects.cat.df <- projects.cat.df %>% dplyr::distinct(.data[["entryId"]], .keep_all = TRUE)
  # process the project information
  projects.cat.final <- ProcessHCAProject(x.df = projects.cat.df)
  return(projects.cat.final)
}

# used in hca, process dataset information of project
ProcessHCAProject <- function(x.df) {
  # entryid and catalog
  entryId <- x.df$entryId
  catalog <- x.df$catalog

  # procotol information
  x.df.protocol <- x.df$protocols[[1]]
  workflow <- HCAPasteCol(x.df.protocol, col = "workflow")
  libraryConstructionApproach <- HCAPasteCol(x.df.protocol, col = "libraryConstructionApproach")
  nucleicAcidSource <- HCAPasteCol(x.df.protocol, col = "nucleicAcidSource")
  instrumentManufacturerModel <- HCAPasteCol(x.df.protocol, col = "instrumentManufacturerModel")
  pairedEnd <- HCAPasteCol(x.df.protocol, col = "pairedEnd")

  # source
  x.df.source <- x.df$sources[[1]]
  sourceId <- HCAPasteColdf(x.df.source, col = "sourceId")
  sourceSpec <- HCAPasteColdf(x.df.source, col = "sourceSpec")

  # project
  x.df.projects <- x.df$projects[[1]]
  projectId <- HCAPasteColdf(x.df.projects, col = "projectId")
  projectTitle <- HCAPasteColdf(x.df.projects, col = "projectTitle")
  projectShortname <- HCAPasteColdf(x.df.projects, col = "projectShortname")
  laboratory <- HCAPasteCol(x.df.projects, col = "laboratory")
  estimatedCellCount <- HCAPasteColdf(x.df.projects, col = "estimatedCellCount")
  projectDescription <- HCAPasteColdf(x.df.projects, col = "projectDescription")
  publications <- HCAPasteColdf(x.df.projects$publications[[1]], col = "publicationTitle")
  accessions <- HCAPasteColdf(x.df.projects$accessions[[1]], col = "accession")
  accessible <- HCAPasteColdf(x.df.projects, col = "accessible")

  # sample
  x.df.samples <- x.df$samples[[1]]
  sampleEntityType <- HCAPasteCol(df = x.df.samples, col = "sampleEntityType")
  organ <- HCAPasteCol(df = x.df.samples, col = "effectiveOrgan")
  sampleID <- HCAPasteCol(df = x.df.samples, col = "id")
  organPart <- HCAPasteCol(df = x.df.samples, col = "organPart")
  disease <- HCAPasteCol(df = x.df.samples, col = "disease")
  preservationMethod <- HCAPasteCol(df = x.df.samples, col = "preservationMethod")

  # cell line
  x.df.cellLines <- x.df$cellLines[[1]]
  cellLineID <- HCAPasteCol(df = x.df.cellLines, col = "id")
  cellLineType <- HCAPasteCol(df = x.df.cellLines, col = "cellLineType")
  cellLinemodelOrgan <- HCAPasteCol(df = x.df.cellLines, col = "modelOrgan")

  # Organism
  x.df.organisms <- x.df$donorOrganisms[[1]]
  donorCount <- ifelse(is.null(x.df.organisms$donorCount), "", x.df.organisms$donorCount)
  developmentStage <- HCAPasteCol(df = x.df.organisms, col = "developmentStage")
  genusSpecies <- HCAPasteCol(df = x.df.organisms, col = "genusSpecies")
  biologicalSex <- HCAPasteCol(df = x.df.organisms, col = "biologicalSex")

  # organoids
  x.df.organoids <- x.df$organoids[[1]]
  organoidsID <- HCAPasteCol(df = x.df.organoids, col = "id")
  organoidsmodelOrgan <- HCAPasteCol(df = x.df.organoids, col = "modelOrgan")
  organoidsmodelOrganPart <- HCAPasteCol(df = x.df.organoids, col = "modelOrganPart")

  # cellSuspensions
  x.df.cellSuspensions <- x.df$cellSuspensions[[1]]
  selectedCellType <- HCAPasteCol(df = x.df.cellSuspensions, col = "selectedCellType")

  # date
  x.df.date <- x.df$dates[[1]]
  lastModifiedDate <- HCAPasteColdf(df = x.df.date, col = "lastModifiedDate")

  # uuid, file name and file formats
  x.df.dataset <- data.frame()
  if (ncol(x.df.projects$matrices) > 0) {
    x.mat.df <- HCAExtactData(x.df.projects$matrices)
    x.mat.df$source <- "matrices"
    x.df.dataset <- data.table::rbindlist(list(x.df.dataset, x.mat.df), fill = TRUE) %>% as.data.frame()
  }
  if (ncol(x.df.projects$contributedAnalyses) > 0) {
    x.ca.df <- HCAExtactData(x.df.projects$contributedAnalyses)
    x.ca.df$source <- "contributedAnalyses"
    x.df.dataset <- data.table::rbindlist(list(x.df.dataset, x.ca.df), fill = TRUE) %>% as.data.frame()
  }
  if (nrow(x.df.dataset) > 0) {
    x.df.dataset <- x.df.dataset[!is.na(x.df.dataset$uuid), ] %>% dplyr::distinct(.data[["uuid"]], .keep_all = TRUE)
    dataMeta <- HCAPasteColdf(df = x.df.dataset, col = "meta")
    dataDescription <- HCAPasteColdf(df = x.df.dataset, col = "contentDescription")
    dataFormat <- HCAPasteColdf(df = x.df.dataset, col = "format")
    dataName <- HCAPasteColdf(df = x.df.dataset, col = "name")
    dataUUID <- HCAPasteColdf(df = x.df.dataset, col = "uuid")
    dataVersion <- HCAPasteColdf(df = x.df.dataset, col = "version")
  } else {
    dataMeta <- NA
    dataDescription <- NA
    dataFormat <- NA
    dataName <- NA
    dataUUID <- NA
    dataVersion <- NA
  }
  # return final dataframe
  project.df <- data.frame(
    projectTitle = projectTitle, projectId = projectId, projectShortname = projectShortname,
    projectDescription = projectDescription, publications = publications, laboratory = laboratory,
    accessions = accessions, accessible = accessible, estimatedCellCount = estimatedCellCount,
    sampleEntityType = sampleEntityType, organ = organ, organPart = organPart, sampleID = sampleID,
    disease = disease, preservationMethod = preservationMethod, donorCount = donorCount,
    developmentStage = developmentStage, genusSpecies = genusSpecies, biologicalSex = biologicalSex,
    selectedCellType = selectedCellType, catalog = catalog, entryId = entryId, sourceId = sourceId, sourceSpec = sourceSpec,
    workflow = workflow, libraryConstructionApproach = libraryConstructionApproach, nucleicAcidSource = nucleicAcidSource,
    instrumentManufacturerModel = instrumentManufacturerModel, pairedEnd = pairedEnd, cellLineID = cellLineID,
    cellLineType = cellLineType, cellLinemodelOrgan = cellLinemodelOrgan, organoidsID = organoidsID,
    organoidsmodelOrgan = organoidsmodelOrgan, organoidsmodelOrganPart = organoidsmodelOrganPart,
    lastModifiedDate = lastModifiedDate, dataMeta = dataMeta, dataDescription = dataDescription, dataFormat = dataFormat,
    dataName = dataName, dataUUID = dataUUID, dataVersion = dataVersion
  )
  return(project.df)
}

# used in CELLxGENE, Zenodo
LoadRDS2Seurat <- function(out.folder, merge, obs.value.filter = NULL, obs.keys = NULL, include.genes = NULL) {
  rds.files <- list.files(path = out.folder, pattern = "rds$", full.names = TRUE, ignore.case = TRUE)
  if (length(rds.files) > 0) {
    message("There is rds in file.ext and return.seu is TRUE, return SeuratOnject!")
    seu.list <- sapply(X = rds.files, FUN = function(x) {
      tryCatch(
        {
          x.rds <- readRDS(x)
          if (class(x.rds) == "Seurat") {
            if (!is.null(obs.value.filter) || !is.null(obs.keys) || !is.null(include.genes)) {
              x.rds.df <- x.rds@meta.data
              # filter dataset metadata
              ## filter cell's metadata values
              if (!is.null(obs.value.filter)) {
                # TODO: test key of obs.value.filter in colnames(x.rds.df)
                tryCatch(
                  {
                    x.rds.df <- x.rds.df %>% dplyr::filter(eval(rlang::parse_expr(obs.value.filter)))
                    if (nrow(x.rds.df) == 0) {
                      message("Please check the value of obs.value.filter, e.g. the 'oligodendrocyte' in cell_type == 'oligodendrocyte'!")
                      x.rds.df <- x.rds@meta.data
                    }
                  },
                  error = function(cond) {
                    message("Please check obs.value.filter: ", cond)
                  }
                )
              }
              ## filter cell's metadata colnames
              if (!is.null(obs.keys)) {
                valid.obs.keys <- intersect(colnames(x.rds.df), obs.keys)
                invalid.obs.keys <- setdiff(obs.keys, colnames(x.rds.df))
                if (length(invalid.obs.keys) > 0) {
                  message("Detected invalid obs.keys: ", paste(invalid.obs.keys, collapse = ", "))
                }
                if (length(valid.obs.keys) > 0) {
                  x.rds.df <- x.rds.df[valid.obs.keys]
                } else {
                  message("No valid obs.keys, return all obs.keys!")
                  x.rds.df <- x.rds.df
                }
              }
              x.rds <- subset(x = x.rds, cells = rownames(x.rds.df), features = include.genes)
            }
            x.rds
          } else {
            message(x, " is not SeuratObject, skip!")
            NULL
          }
        },
        error = function(cond) {
          message("Reading ", x, " error: ", cond)
          NULL
        }
      )
    })
    if (is.null(seu.list)) {
      seu.obj <- NULL
    } else {
      empty.seu <- seu.list[sapply(seu.list, is.null)]
      if (length(empty.seu) > 0) {
        message(paste(names(empty.seu), collapse = ", "), " is NULL, drop!")
      }
      # remove NULL element
      seu.list[sapply(seu.list, is.null)] <- NULL
      if (isTRUE(merge)) {
        seu.obj <- mergeExperiments(seu.list)
      } else {
        seu.obj <- seu.list
      }
    }
    return(seu.obj)
  } else {
    message("There is no rds file under ", out.folder)
    return(NULL)
  }
}

# parse ENA xml (fastq, bam, sra)
ParseENAxml <- function(run, df.type = c("fastq", "bam")) {
  xml.url <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", run)
  # parse xml
  xml.raw <- curl::curl_fetch_memory(xml.url)
  xml.content <- rawToChar(xml.raw$content)
  if (df.type == "fastq") {
    # extract link
    fastq.text.url <- gsub(pattern = ".*<DB>ENA-FASTQ-FILES</DB>\n +<ID><!\\[CDATA\\[(.*?)\\]\\]></ID>.*", replacement = "\\1", x = xml.content)
    # read file info
    fastq.text <- data.table::fread(fastq.text.url, showProgress = F)
    fastq.text <- fastq.text %>%
      tidyr::separate_rows(.data[["fastq_ftp"]], .data[["fastq_md5"]], .data[["fastq_bytes"]],
        sep = ";", convert = TRUE
      ) %>%
      as.data.frame()
    valid.fastq.text <- fastq.text %>% dplyr::filter(!is.na(.data[["fastq_ftp"]]))
    if (nrow(valid.fastq.text) > 0) {
      return(valid.fastq.text)
    } else {
      message("There is no valid fastq files under run: ", run)
      return(NULL)
    }
  } else if (df.type == "bam") {
    # extract link
    submitted.text.url <- gsub(pattern = ".*<DB>ENA-SUBMITTED-FILES</DB>\n +<ID><!\\[CDATA\\[(.*?)\\]\\]></ID>.*", replacement = "\\1", x = xml.content)
    # read file info
    submitted.text <- data.table::fread(submitted.text.url, showProgress = FALSE)
    submitted.text <- submitted.text %>%
      tidyr::separate_rows(.data[["submitted_ftp"]], .data[["submitted_md5"]], .data[["submitted_bytes"]],
        .data[["submitted_format"]],
        sep = ";", convert = TRUE
      ) %>%
      as.data.frame()
    valid.submitted.text <- submitted.text %>%
      dplyr::filter(!is.na(.data[["submitted_ftp"]])) %>%
      dplyr::filter(grepl(pattern = "bam|bai", x = .data[["submitted_format"]], ignore.case = TRUE))
    if (nrow(valid.submitted.text) > 0) {
      return(valid.submitted.text)
    } else {
      message("There is no valid bam|bai files under run: ", run)
      return(NULL)
    }
  }
}

# download files with url (download.file and ascp)
DownloadMethod <- function(rn, url.vec, name.vec = NULL, out.folder = NULL, download.method = c("download.file", "ascp"),
                           quiet = FALSE, timeout = 3600, ascp.path = NULL, max.rate = "300m", rename = FALSE) {
  # download
  if (download.method == "download.file") {
    # set timeout
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
    down.status <- utils::download.file(url = url.vec, destfile = file.path(out.folder, name.vec), quiet = quiet, mode = "wb")
    # failed download
    fail.status <- which(down.status != 0)
    if (length(fail.status) == 0) {
      message("Download successful: ", rn)
      return(NULL)
    } else {
      warning("Run download.file error: ", rn)
      return(rn)
    }
  } else {
    # get ascp path
    if (is.null(ascp.path)) {
      # specify ascp path
      ascp.path <- Sys.which("ascp")
      if (ascp.path == "") {
        stop("Can not find ascp automatically, please specify the path!")
      }
    } else {
      ascp.path <- ascp.path
    }
    # get asperaweb_id_dsa.openssh
    ascp.pubkey <- gsub("bin/ascp$", "etc/asperaweb_id_dsa.openssh", ascp.path)
    ascp.cmd <- paste(ascp.path, "-QT -k1 -l", max.rate, "-P33001", "-i", ascp.pubkey, url.vec, out.folder)
    ascp.cmd <- paste(ascp.cmd, collapse = " && ")
    message(paste("Calling ascp:", ascp.cmd))
    ascp.status <- system(ascp.cmd, intern = TRUE)
    ascp.status.code <- attr(ascp.status, "status")
    if (!is.null(ascp.status.code)) {
      warning("Run ascp error: ", rn)
      do.call(file.remove, list(list.files(out.folder, full.names = TRUE, pattern = "partial$|aspera-ckpt$")))
      return(rn)
    } else {
      message("Download successful: ", rn)
      # rename file
      if (isTRUE(rename)) {
        file.rename(file.path(out.folder, basename(url.vec)), file.path(out.folder, name.vec))
      }
      return(NULL)
    }
  }
}

# gunzip a file
# source: https://github.com/seandavi/GEOquery/blob/026c655561dbcf1b99551a8093642750c48ed038/R/getGEOfile.R
Gunzip <- function(filename, destname = gsub("[.]gz$", "", filename), overwrite = FALSE, remove = TRUE, BFR.SIZE = 1e7) {
  if (filename == destname) {
    stop(sprintf("Argument 'filename' and 'destname' are identical: %s", filename))
  }
  if (!overwrite && file.exists(destname)) {
    stop(sprintf("File already exists: %s", destname))
  }

  inn <- gzfile(filename, "rb")
  on.exit(if (!is.null(inn)) close(inn))

  out <- file(destname, "wb")
  on.exit(close(out), add = TRUE)

  nbytes <- 0
  repeat {
    bfr <- readBin(inn, what = raw(0), size = 1, n = BFR.SIZE)
    n <- length(bfr)
    if (n == 0) {
      break
    }
    nbytes <- nbytes + n
    writeBin(bfr, con = out, size = 1)
  }

  if (remove) {
    close(inn)
    inn <- NULL
    file.remove(filename)
  }

  invisible(nbytes)
}

# used in ParseGEO, read count matrix file of bulk and smart-seq2 RNA-seq.
ReadFile <- function(file.path, extra.cols = c(
                       "chr", "start", "end", "strand", "length",
                       "width", "chromosome", "seqnames", "seqname",
                       "chrom", "chromosome_name", "seqid", "stop"
                     ),
                     transpose = TRUE) {
  file.ext <- tools::file_ext(file.path)
  if (file.ext %in% c("xlsx", "xls")) {
    # read excel file
    count.mat <- tryCatch(
      {
        readxl::read_excel(path = file.path) %>% as.data.frame()
      },
      error = function(e) {
        message(e)
        # empty dataframe
        data.frame()
      }
    )
  } else if (file.ext %in% c("csv", "tsv", "txt", "tab")) {
    # read text file
    count.mat <- tryCatch(
      {
        data.table::fread(file = file.path) %>% as.data.frame()
      },
      error = function(e) {
        message(e)
        # empty dataframe
        data.frame()
      }
    )
  }
  if (nrow(count.mat) > 0) {
    # the first column must be gene
    # check unique
    if (length(count.mat[, 1]) != length(unique(count.mat[, 1]))) {
      message(
        "The row names are not unique, run 'make.unique'.", "\n",
        "You can extract duplicated gene names with grep(pattern = '_dup', x = rownames(count.mat), value = T)."
      )
      # ensure unique
      count.mat[, 1] <- make.unique(count.mat[, 1], sep = "_dup")
    }
    rownames(count.mat) <- count.mat[, 1]
    count.mat[, 1] <- NULL
    # remove extra columns
    count.mat.cols <- colnames(count.mat) %>% tolower()
    extra.cols.idx <- sapply(count.mat.cols, function(x) {
      x %in% extra.cols
    })
    extra.cols.remove <- count.mat.cols[extra.cols.idx]
    if (length(extra.cols.remove) > 0) {
      message("Detect and remove extra columns: ", paste(extra.cols.remove, collapse = ", "))
    }
    count.mat <- count.mat[, !extra.cols.idx, drop = FALSE]
    # check col type
    count.mat.ct <- sapply(count.mat, class)
    # get character columns
    count.mat.ct.c <- count.mat.ct[count.mat.ct == "character"] %>% names()
    if (length(count.mat.ct.c) > 0) {
      message("Detect and remove character columns: ", paste(count.mat.ct.c, collapse = ", "))
      count.mat[count.mat.ct.c] <- NULL
    }
    # check row and col number
    if (nrow(count.mat) < ncol(count.mat)) {
      warning(
        "The number of rows: ", nrow(count.mat), "(", paste(head(rownames(count.mat)), collapse = ", "), ")",
        " is smaller than the number of columns: ", ncol(count.mat), ". May be a transposed matrix, please check and set parameter transpose (Default: TRUE)!"
      )
      if (isTRUE(transpose)) {
        message("The transpose is set to TRUE, transpose the matrix!")
        count.mat <- t(count.mat) %>% as.data.frame()
      }
    }
  }
  return(count.mat)
}

# modified from https://github.com/seandavi/GEOquery/blob/c85b1115cdee87515ed0d72bfde86f3986a0b2b7/R/getGEOSuppFiles.R
getDirListing <- function(url) {
  # Takes a URL and returns a character vector of filenames
  a <- xml2::read_html(url)
  fnames <- grep("^G", xml2::xml_text(xml2::xml_find_all(a, "//a/@href")), value = TRUE)
  return(fnames)
}
# get file name and url
getGEOFileIndexName <- function(GEO, index) {
  geotype <- toupper(substr(GEO, 1, 3))
  stub <- gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  if (geotype == "GSM") {
    url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/",
      stub, GEO
    )
  }
  if (geotype == "GSE") {
    url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",
      stub, GEO
    )
  }
  if (geotype == "GPL") {
    url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/",
      stub, GEO
    )
  }
  fnames <- try(getDirListing(url), silent = TRUE)
  # check status
  if (inherits(fnames, "try-error")) {
    stop("No supplemental files found. Check URL manually if in doubt: ", url)
  }
  # check index
  if (index > length(fnames)) {
    stop("Please provide valid supplementary file index. Total length: ", length(fnames))
  }
  used.fnames <- fnames[index]
  return(c(used.fnames, url))
}
# https://github.com/seandavi/GEOquery/blob/c85b1115cdee87515ed0d72bfde86f3986a0b2b7/R/getGEOSuppFiles.R
getGEOSuppFilesInner <- function(GEO, makeDirectory = TRUE, baseDir = getwd(), index = 1) {
  # get used file
  index.url <- getGEOFileIndexName(GEO = GEO, index = index)
  storedir <- baseDir
  fileinfo <- list()
  # create output folder
  if (makeDirectory) {
    suppressWarnings(dir.create(storedir <- file.path(baseDir, GEO)))
  }
  # download
  destfile <- file.path(storedir, index.url[1])
  result <- tryCatch(
    {
      if (!file.exists(destfile)) {
        res <- download.file(
          paste(file.path(index.url[2], index.url[1]),
            "tool=geoquery",
            sep = "?"
          ),
          destfile = destfile,
          mode = "wb", method = getOption("download.file.method.GEOquery")
        )
      } else {
        message(sprintf(
          "Using locally cached version of supplementary file(s) %s found here:\n%s ",
          GEO, destfile
        ))
        res <- 0
      }
      res == 0
    },
    error = function(e) {
      return(FALSE)
    },
    warning = function(w) {
      return(FALSE)
    }
  )
  if (!result) {
    if (file.exists(destfile)) {
      file.remove(destfile)
    }
    stop(sprintf("Failed to download %s!", destfile))
  }
  return(file.info(destfile))
}

# extract file extension without compression
getExtWithoutCompression <- function(filename) {
  # Define common compression extensions
  compression.exts <- c("gz", "bz2", "xz", "zip", "tgz", "tar.gz") # expanded list
  # Iterate and remove compression extensions
  for (ext in compression.exts) {
    if (grepl(paste0("\\.", ext, "$"), filename)) {
      filename <- sub(paste0("\\.", ext, "$"), "", filename)
    }
  }
  # Now, extract the remaining extension using tools::file_ext
  return(tools::file_ext(filename))
}

# recursively move files to top directory
MoveFileRecursively <- function(folder) {
  move.folder.log <- sapply(folder, function(x) {
    x.files <- list.files(path = x, full.names = TRUE, recursive = TRUE)
    move.file.log <- sapply(x.files, function(y) {
      # get file name
      y.name <- basename(y)
      y.last.folder <- basename(dirname(y))
      # get gse or gsm
      gname <- gsub(pattern = "(GS[EM][0-9]+).*", replacement = "\\1", x = basename(x))
      new.file <- file.path(dirname(x), paste(gname, y.last.folder, y.name, sep = "_"))
      if (y != new.file) {
        # move file
        copy.tag <- file.copy(from = y, to = new.file)
        # remove the original file
        remove.tag <- file.remove(y)
        copy.tag
      }
    })
  })
  return(move.folder.log)
}

# untar and recursively move files to top directory
ProcessTAR <- function(tar.files, remove.cf = TRUE) {
  untar.log <- sapply(
    tar.files,
    function(x) {
      x.folder <- gsub(pattern = ".tar", replacement = "", x = x)
      utils::untar(x, exdir = x.folder)
    }
  )
  # recursively move files
  all.tar.folders <- gsub(pattern = ".tar$", replacement = "", x = tar.files)
  move.folder.log <- MoveFileRecursively(all.tar.folders)
  # remove tar.gz files
  if (remove.cf) {
    rm.log <- sapply(tar.files, file.remove)
  }
  return(move.folder.log)
}

# process compressed files in zip, tar, tar.gz
ProcessCompressedFiles <- function(all.files, remove.cf = TRUE) {
  # deal with zip files
  all.zip.files <- grep(pattern = "zip$", x = all.files, value = TRUE)
  if (length(all.zip.files) > 0) {
    message("Detect files in zip format, extract!")
    unzip.log <- sapply(
      all.zip.files,
      function(x) {
        x.folder <- gsub(pattern = ".zip", replacement = "", x = x)
        unzip(x, exdir = x.folder, overwrite = TRUE)
      }
    )
    # recursively move files
    all.zip.folders <- gsub(pattern = ".zip$", replacement = "", x = all.zip.files)
    move.folder.log <- MoveFileRecursively(all.zip.folders)
    # remove zip files
    if (remove.cf) {
      rm.log <- sapply(all.zip.files, file.remove)
    }
  }
  # deal with tar.gz files
  all.targz.files <- grep(pattern = "tar.gz$", x = all.files, value = TRUE)
  if (length(all.targz.files) > 0) {
    message("Detect files in tar.gz format, extract!")
    # unzip
    unzip.log <- sapply(
      all.targz.files,
      function(x) {
        Gunzip(x, overwrite = TRUE)
      }
    )
    # untar
    all.tar.files <- gsub(pattern = ".gz$", replacement = "", x = all.targz.files)
    process.tar <- ProcessTAR(tar.files = all.tar.files, remove.cf = remove.cf)
  }
  # deal with tar files
  all.tar.files2 <- grep(pattern = "tar$", x = all.files, value = TRUE)
  if (length(all.tar.files2) > 0) {
    message("Detect files in tar format, extract!")
    process.tar2 <- ProcessTAR(tar.files = all.tar.files2, remove.cf = remove.cf)
  }
}

# process 10x files: deal with compressed files, identify h5 and MEX format, rename and move
Process10xFiles <- function(acce, folder, accept.fmt, out.folder, gene2feature) {
  # deal with compressed files
  all.files <- list.files(folder, full.names = TRUE)
  process.compressed.log <- ProcessCompressedFiles(all.files = all.files, remove.cf = TRUE)
  # identify accept files in given format
  if ("MEX" %in% accept.fmt) {
    valid.pat.mex <- "barcodes.*.tsv$|genes.*.tsv$|matrix.*.mtx$|features.*.tsv$"
    all.files.mex <- list.files(folder, full.names = TRUE, pattern = valid.pat.mex)
    gzip.log <- sapply(
      all.files.mex,
      function(x) {
        R.utils::gzip(filename = x, remove = TRUE)
      }
    )
    # sample name is after, rename
    valid.pat.gz <- "barcodes.+.tsv.gz$|genes.+.tsv.gz$|matrix.+.mtx.gz$|features.+.tsv.gz$"
    all.files.gz <- list.files(folder, full.names = TRUE, pattern = valid.pat.gz)
    if (length(all.files.gz) > 0) {
      rename.log <- sapply(all.files.gz, function(x) {
        x.new <- gsub(pattern = "(.*)(barcodes|genes|matrix|features)(.*)(.tsv.gz|.mtx.gz)", replacement = "\\1\\3\\2\\4", x = x)
        file.rename(from = x, to = x.new)
      })
    }
    # recognize valid files: barcodes.tsv.gz, genes.tsv.gz, matrix.mtx.gz and features.tsv.gz
    valid.pat.gz <- "barcodes.tsv.gz$|genes.tsv.gz$|matrix.mtx.gz$|features.tsv.gz$"
    all.files.gz <- list.files(folder, full.names = TRUE, pattern = valid.pat.gz)
  } else {
    all.files.gz <- c()
  }
  if ("h5" %in% accept.fmt) {
    h5.gz.file <- list.files(folder, full.names = TRUE, pattern = "h5.gz$")
    if (length(h5.gz.file) > 0) {
      # gunzip file
      unzip.log <- sapply(
        h5.gz.file,
        function(x) {
          Gunzip(x, overwrite = TRUE)
        }
      )
    }
    all.files.h5 <- list.files(folder, full.names = TRUE, pattern = "h5$")
  } else {
    all.files.h5 <- c()
  }
  # prepare out folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  out.folder <- file.path(out.folder, acce) # optimize output folder
  # rename and move files
  if (length(all.files.gz) > 0) {
    message("Detect ", length(all.files.gz), " files in MEX(barcode/feature/gene/matrix) format.")
    # change file name
    if (gene2feature) {
      change.name.log <- sapply(all.files.gz, function(x) {
        if (grepl(pattern = "genes.tsv.gz$", x = x)) {
          new.name <- gsub(pattern = "genes.tsv.gz$", replacement = "features.tsv.gz", x = x)
          file.rename(from = x, to = new.name)
        }
      })
      all.files.gz <- list.files(folder, full.names = TRUE, pattern = valid.pat.gz)
    }
    # get folder
    all.sample.folder <- sapply(all.files.gz, function(x) {
      # get basename and dirname
      file.name <- basename(x)
      dir.name <- dirname(x)
      # remove file type tag
      file.name <- gsub(pattern = valid.pat.gz, replacement = "", x = file.name)
      # remove possible _ and .
      file.name <- gsub(pattern = "[_-.]$", replacement = "", x = file.name)
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
  if (length(all.files.h5) > 0) {
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
  if (length(all.files.gz) == 0 && length(all.files.h5) == 0) {
    stop("No valid 10x format (barcode/feature/gene/matrix, h5) files detected! Please check.")
  }
  message("Process 10x fiels done! All files are in ", out.folder)
  return(NULL)
}

# used in CheckAPIs, whether the url exists: https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
CheckURL <- function(url_in, t = 2) {
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = T)[1])
  suppressWarnings(try(close.connection(con), silent = T))
  ifelse(is.null(check), TRUE, FALSE)
}
