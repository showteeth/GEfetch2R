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

# used in hca, recursively extract projects (limit size to 100, get all projects when the size is greater than 100)
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
