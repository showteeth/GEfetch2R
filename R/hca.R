#' Show All Available Catalogs in Human Cell Atlas.
#'
#' @return Named vector of catalogs.
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @examples
#' \donttest{
#' # all available catalogs
#' all.hca.catalogs <- ExtractHCACatalogs()
#' }
ExtractHCACatalogs <- function() {
  # base urls
  hca.base.url <- "https://service.azul.data.humancellatlas.org"
  # get all catalogs
  catalog.url <- paste0(hca.base.url, "/index/catalogs")
  catalog.content <- URLRetrieval(catalog.url)
  catalog.vec <- sapply(names(catalog.content$catalogs), function(x) {
    if (!catalog.content$catalogs[[x]]$internal) x
  }) %>% unlist()
  return(catalog.vec)
}

# extract all projects
ExtractHCAProjects <- function(catalog = NULL) {
  # base urls
  hca.base.url <- "https://service.azul.data.humancellatlas.org"
  # get all catalogs
  catalog.vec <- ExtractHCACatalogs()
  if (!is.null(catalog)) {
    catalog.vec <- intersect(catalog, catalog.vec)
    if (length(catalog.vec) == 0) {
      stop("Please check catalog you proided! All available catalogs can be accessed via ExtractHCACatalogs().")
    }
  }
  # get catalog project list
  hca.projects.list <- lapply(catalog.vec, function(x) {
    cat.prj.url <- paste0(hca.base.url, "/index/projects?catalog=", x, "&size=75")
    cat.prj <- RecurURLRetrieval(cat.prj.url) %>% as.data.frame()
    cat.prj$catalog <- x
    return(cat.prj)
  })
  # get catalog project df
  hca.projects.df <- data.table::rbindlist(hca.projects.list, fill = TRUE) %>% as.data.frame()
  # remove duplicated projects in different catalogs
  hca.projects.df <- hca.projects.df %>% dplyr::distinct(.data[["entryId"]], .keep_all = TRUE)
  return(hca.projects.df)
}

#' Show All Available Projects in Human Cell Atlas.
#'
#' @param catalog The catalog of the projects (one or multiple values). Different catalogs may share some projects.
#' All available catalogs can be accessed via \code{ExtractHCACatalogs}. Default: NULL (all catalogs, remove duplicated projects).
#'
#' @return Dataframe contains all available projects.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @importFrom dplyr distinct
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' # all available projects
#' all.hca.projects <- ShowHCAProjects()
#' }
ShowHCAProjects <- function(catalog = NULL) {
  # get all projects information
  hca.projects.df <- ExtractHCAProjects(catalog = catalog)
  # get project detail information
  hca.projects.detail.list <- lapply(1:nrow(hca.projects.df), function(x) {
    x.df <- hca.projects.df[x, ]
    ProcessHCAProject(x.df)
  })
  hca.projects.detail.df <- data.table::rbindlist(hca.projects.detail.list, fill = TRUE) %>% as.data.frame()
  return(hca.projects.detail.df)
}

#' Extract Metadata of Human Cell Atlas Projects with Attributes.
#'
#' @param all.projects.df All detail information of HCA projects, obtained with \code{ShowHCAProjects}.
#' @param organism The organism of the projects, choose from "Homo sapiens", "Mus musculus",
#' "Macaca mulatta", "canis lupus familiaris", one or multiple values. Default: NULL (All).
#' @param sex The sex of the projects, choose from "female", "male", "mixed", "unknown",
#' one or multiple values. Default: NULL (All).
#' @param organ The organ of the projects (e.g. brain), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param organ.part The organ part of the projects (e.g. cortex), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param disease The disease of the projects (e.g. normal), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param sample.type The sex of the projects, choose from "specimens", "organoids", "cellLines",
#' one or multiple values. Default: NULL (All).
#' @param preservation.method The preservation method of the projects (e.g. fresh), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param protocol The protocol of the projects (e.g. 10x 3' v2), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param suspension.type The suspension type of the projects, choose from "single cell", "single nucleus", "bulk cell", "bulk nuclei",
#' one or multiple values. Default: NULL (All).
#' @param cell.type The cell type of the projects (e.g. neuron), obtain available values with \code{StatDBAttribute},
#' one or multiple values. Default: NULL (All).
#' @param cell.num Cell number filter. If NULL, no filter; if one value, lower filter; if two values, low and high filter.
#' Deault: NULL(without filtering).
#' @param sequencing.type The sequencing instrument type of the projects (e.g. illumina hiseq 2500),
#' obtain available values with \code{StatDBAttribute}, one or multiple values. Default: NULL (All).
#' @param remove.nodata Logical value, whether to remove project with no downloadable data. Default: TRUE.
#'
#' @return Dataframe contains filtered projects.
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @export
#' @references https://bioconductor.org/packages/release/bioc/html/hca.html
#'
#' @examples
#' \donttest{
#' # all available projects
#' all.hca.projects <- ShowHCAProjects()
#' # all human projects
#' all.human.projects <- ExtractHCAMeta(all.projects.df = all.hca.projects, organism = "Homo sapiens")
#' # all human and 10x 3' v2
#' all.human.10x.projects <- ExtractHCAMeta(
#'   all.projects.df = all.hca.projects,
#'   organism = "Homo sapiens",
#'   protocol = c("10x 3' v2", "10x 3' v3")
#' )
#' }
ExtractHCAMeta <- function(all.projects.df, organism = NULL, sex = NULL, organ = NULL, organ.part = NULL, disease = NULL,
                           sample.type = NULL, preservation.method = NULL, protocol = NULL,
                           suspension.type = NULL, cell.type = NULL, cell.num = NULL, sequencing.type = NULL, remove.nodata = TRUE) {
  # all projects detail dataframe
  hca.projects.detail.df <- all.projects.df
  # extract row index under different filter
  organism.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "genusSpecies", attr.value = organism)
  sex.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "biologicalSex", attr.value = sex)
  organ.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "organ", attr.value = organ)
  organ.part.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "organPart", attr.value = organ.part)
  disease.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "disease", attr.value = disease)
  sample.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "sampleEntityType", attr.value = sample.type)
  preservation.method.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "preservationMethod", attr.value = preservation.method)
  protocol.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "libraryConstructionApproach", attr.value = protocol)
  suspension.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "nucleicAcidSource", attr.value = suspension.type)
  cell.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "selectedCellType", attr.value = cell.type)
  sequencing.type.idx <- HCAAttrFilter(df = hca.projects.detail.df, attr = "instrumentManufacturerModel", attr.value = sequencing.type)
  if (is.null(cell.num)) {
    cnum.idx <- 1:nrow(hca.projects.detail.df)
  } else if (length(cell.num) == 1) {
    cnum.idx <- which(hca.projects.detail.df$estimatedCellCount > as.numeric(cell.num))
  } else {
    cnum.idx <- which(hca.projects.detail.df$estimatedCellCount > as.numeric(cell.num[1]) &
      hca.projects.detail.df$estimatedCellCount < as.numeric(cell.num[2]))
  }
  # filter on the whole dataset
  valid.idx <- Reduce(intersect, list(
    organism.idx, sex.idx, organ.idx, organ.part.idx, disease.idx, sample.type.idx,
    preservation.method.idx, protocol.idx, suspension.type.idx, cell.type.idx, sequencing.type.idx, cnum.idx
  ))
  used.sample.df <- hca.projects.detail.df[valid.idx, ]
  rownames(used.sample.df) <- NULL
  # remove project with no downloadable data
  if (remove.nodata) {
    remove.sample.df <- used.sample.df[is.na(used.sample.df$dataUUID), ]
    if (nrow(remove.sample.df) > 0) {
      message("Remove ", nrow(remove.sample.df), " project(s) with no downloadable data!")
      used.sample.df <- used.sample.df[!is.na(used.sample.df$dataUUID), ]
    }
  }
  return(used.sample.df)
}

#' Download Human Cell Atlas Datasets.
#'
#' @param meta Metadata used to download, can be from \code{ExtractHCAMeta}.
#' Skip when \code{link} is not NULL. Default: NULL.
#' @param link Vector contains project link(s), e.g. "https://explore.data.humancellatlas.org/projects/902dc043-7091-445c-9442-d72e163b9879".
#' Skip when \code{meta} is not NULL. Default: NULL.
#' @param file.ext The valid file extension for download. When NULL, use "rds", "rdata", "h5", "h5ad", "loom", "tsv".
#' Default: c("rds", "rdata", "h5", "h5ad", "loom", "tsv").
#' @param out.folder The output folder. Default: NULL (current working directory).
#' @param timeout Maximum request time. Default: 3600.
#' @param quiet Logical value, whether to show downloading progress. Default: FALSE (show).
#' @param parallel Logical value, whether to download parallelly. Default: TRUE. When "libcurl" is available for \code{download.file},
#' the parallel is done by default (\code{parallel} can be FALSE).
#' @param use.cores The number of cores used. Default: NULL (the minimum value of
#' extracted \code{length(download.urls)} and \code{parallel::detectCores()}).
#' @param return.seu Logical value, whether to load downloaded datasets to Seurat. Valid when rds in \code{file.ext} and all
#' datasets download successfully. Default: FALSE.
#' @param merge Logical value, whether to merge Seurat list when there are multiple rds files,
#' used when \code{return.seu} is TRUE. Default: FALSE.
#'
#' @return SeuratObject (\code{return.seu} is TRUE, rds in \code{file.ext}) or
#' list contains files' metadata of downloaded successfully (down.meta) and failed (fail.meta).
#' @importFrom magrittr %>%
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom data.table rbindlist
#' @importFrom parallel detectCores mclapply
#' @importFrom utils download.file
#' @importFrom tidyr spread
#' @importFrom dplyr distinct relocate any_of last_col select everything filter
#' @importFrom rlang .data parse_expr
#' @export
#'
#' @examples
#' \dontrun{
#' # all available projects
#' all.hca.projects <- ShowHCAProjects()
#' # all human and 10x 3' v2
#' all.human.10x.projects <- ExtractHCAMeta(
#'   all.projects.df = all.hca.projects,
#'   organism = "Homo sapiens",
#'   protocol = c("10x 3' v2", "10x 3' v3")
#' )
#' # download, need users to provide the output folder
#' ParseHCA(meta = all.human.10x.projects, out.folder = "/path/to/output")
#' # download given projects
#' ParseHCA(
#'   link = c(
#'     "https://explore.data.humancellatlas.org/projects/902dc043-7091-445c-9442-d72e163b9879",
#'     "https://explore.data.humancellatlas.org/projects/cdabcf0b-7602-4abf-9afb-3b410e545703"
#'   ),
#'   out.folder = "/path/to/output"
#' )
#' }
ParseHCA <- function(meta = NULL, link = NULL, file.ext = c("rds", "rdata", "h5", "h5ad", "loom", "tsv"), out.folder = NULL,
                     timeout = 3600, quiet = FALSE, parallel = TRUE, use.cores = NULL, return.seu = FALSE, merge = TRUE) {
  # file.ext: ignore case, tar.gz, gz
  if (is.null(file.ext)) {
    warning("There is no file extension provided, use all valid (rds, rdata, h5, h5ad, loom and tsv).")
    file.ext <- c("rds", "rdata", "h5", "h5ad", "loom", "tsv")
  }
  file.ext <- intersect(file.ext, c("rds", "rdata", "h5", "h5ad", "loom", "tsv"))
  if (length(file.ext) == 0) {
    stop("Please provide valid file extension: rds, rdata, h5, h5ad and loom.")
  }
  # check meta
  if (is.null(meta)) {
    if (is.null(link)) {
      stop("The meta and link are NULL, please provide at least one valid value!")
    } else {
      projects.vec <- sapply(link, basename)
      meta.list <- lapply(projects.vec, ParseHCAdataset)
      meta <- data.table::rbindlist(meta.list, fill = TRUE) %>%
        as.data.frame()
    }
  }
  # check columns
  CheckColumns(df = meta, columns = c("dataUUID", "dataFormat", "dataName", "dataMeta", "dataDescription"))
  # restore project dataset
  projects.datasets.list <- lapply(1:nrow(meta), function(x) {
    x.df <- meta[x, ]
    x.entryId <- x.df$entryId
    x.catalog <- x.df$catalog
    x.dataMeta <- strsplit(x = x.df$dataMeta, split = ", ")[[1]]
    x.dataDescription <- strsplit(x = x.df$dataDescription, split = ", ")[[1]]
    x.dataUUID <- strsplit(x = x.df$dataUUID, split = ", ")[[1]]
    x.dataFormat <- strsplit(x = x.df$dataFormat, split = ", ")[[1]]
    x.dataName <- strsplit(x = x.df$dataName, split = ", ")[[1]]
    x.data.df <- data.frame(meta = x.dataMeta, contentDescription = x.dataDescription, uuid = x.dataUUID, format = x.dataFormat, name = x.dataName)
    x.data.df$entryId <- x.entryId
    x.data.df$catalog <- x.catalog
    # filter file extension
    if (nrow(x.data.df) > 0) {
      x.data.df$entryId <- x.df$entryId
      # filter with file.ext
      file.ext <- c(file.ext, paste0(file.ext, ".tar.gz"), paste0(file.ext, ".gz"))
      x.data.df$lowerformat <- tolower(x.data.df$format)
      x.data.valid.df <- x.data.df[x.data.df$lowerformat %in% file.ext, ]
      if (nrow(x.data.valid.df) == 0) {
        message(
          "There is no file in entryId: ", x.df$entryId, " with extension: ", paste(file.ext, collapse = ", "), ". Available file.ext: ",
          paste(unique(x.data.df$lowerformat), collapse = ", "), "."
        )
      }
    } else {
      message("There is no file to download in entryId: ", x.df$entryId, ".")
      x.data.valid.df <- x.data.df
    }
    return(x.data.valid.df)
  })
  projects.datasets.df <- data.table::rbindlist(projects.datasets.list, fill = TRUE) %>% as.data.frame()
  if (nrow(projects.datasets.df) == 0) {
    stop(
      "There is no file in entryId: ", paste(projects.valid$entryId, collapse = ", "), " with extension: ",
      paste(file.ext, collapse = ", "), ". Please check the file.ext!"
    )
  } else {
    # get url with uuid
    projects.datasets.df$url <- paste0("https://service.azul.data.humancellatlas.org/repository/files/", projects.datasets.df$uuid)
    # add name
    projects.datasets.df$name <- sapply(1:nrow(projects.datasets.df), function(x) {
      x.pdvd <- projects.datasets.df[x, ]
      ifelse(is.null(x.pdvd$name),
        paste0(make.names(x.pdvd$meta), ".", x.pdvd$format),
        ifelse(is.na(x.pdvd$name),
          paste0(make.names(x.pdvd$meta), ".", x.pdvd$format),
          ifelse(x.pdvd$name == "",
            paste0(make.names(x.pdvd$meta), ".", x.pdvd$format),
            x.pdvd$name
          )
        )
      )
    })
    # add metadata
    projects.datasets.df <- merge(meta[c(
      "projectTitle", "projectDescription", "publications",
      "sampleEntityType", "organPart", "disease", "preservationMethod", "biologicalSex",
      "nucleicAcidSource", "entryId", "catalog"
    )], projects.datasets.df, by = c("entryId", "catalog"))
    projects.datasets.df <- projects.datasets.df %>%
      dplyr::relocate(dplyr::any_of(c("entryId", "catalog")), .after = dplyr::last_col()) %>%
      dplyr::select(dplyr::any_of(c("meta", "contentDescription", "name")), dplyr::everything())
    projects.datasets.df$lowerformat <- NULL
    # change colnames
    colnames(projects.datasets.df) <- c(
      "dataMeta", "dataDescription", "dataName", "projectTitle", "projectDescription",
      "publications", "sampleEntityType", "organPart", "disease", "preservationMethod",
      "biologicalSex", "nucleicAcidSource", "dataUUID", "dataFormat", "url",
      "entryId", "catalog"
    )
    # prepare download urls
    download.urls <- projects.datasets.df$url
    names(download.urls) <- projects.datasets.df$dataName
    # prepare output folder
    if (is.null(out.folder)) {
      out.folder <- getwd()
    }
    if (!dir.exists(out.folder)) {
      message(out.folder, " does not exist, create automatically!")
      dir.create(out.folder, recursive = TRUE)
    }
    names(download.urls) <- file.path(out.folder, names(download.urls))
    # download urls
    # set timeout
    env.timeout <- getOption("timeout")
    on.exit(options(timeout = env.timeout)) # restore timeout
    options(timeout = timeout)
    message("Start downloading!")
    if (isTRUE(parallel)) {
      # prepare cores
      cores.used <- min(parallel::detectCores(), length(download.urls), use.cores)
      down.status <- parallel::mclapply(X = 1:length(download.urls), FUN = function(x) {
        utils::download.file(url = download.urls[x], destfile = names(download.urls)[x], quiet = quiet, mode = "wb")
      }, mc.cores = cores.used)
    } else {
      down.status <- utils::download.file(url = download.urls, destfile = names(download.urls), quiet = quiet, mode = "wb")
    }
    # process failed datasets
    down.status <- unlist(down.status)
    fail.status <- which(down.status != 0)
    if (length(fail.status) == 0) {
      message("All datasets downloaded successfully!")
      if (isTRUE(return.seu)) {
        rds.gz.files <- list.files(path = out.folder, pattern = "rds.gz$", full.names = TRUE, ignore.case = TRUE)
        if (length(rds.gz.files) > 0) {
          message("Detect zip files: ", paste(rds.gz.files, collapse = ", "), ". Unzip!")
          # unzip
          unzip.log <- sapply(
            rds.gz.files,
            function(x) {
              Gunzip(x, overwrite = TRUE)
            }
          )
        }
        seu.obj <- LoadRDS2Seurat(out.folder = out.folder, merge = merge)
        return(seu.obj)
      } else {
        res.list <- list(down.meta = projects.datasets.df, fail.meta = NULL)
        return(res.list)
      }
    } else {
      message(length(fail.status), " files downloaded failed, please re-run with fail.meta (meta)")
      fail.entry.id <- projects.datasets.df[fail.status, "entryId"] %>% unique()
      # for re-run
      fail.meta <- meta[meta$entryId %in% fail.entry.id, ]
      success.meta <- projects.datasets.df[!projects.datasets.df$entryId %in% fail.entry.id, ]
      res.list <- list(down.meta = success.meta, fail.meta = fail.meta)
      return(res.list)
    }
  }
}
