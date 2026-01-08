#' Check the Availability of APIs.
#'
#' @param database Vector of databases for API checking. Default: c("SRA/ENA", "GEO", "PanglanDB", "UCSC Cell Browser", "Zenodo", "CELLxGENE", "Human Cell Atlas").
#'
#' @return NULL.
#' @importFrom GEOfastq crawl_gsms
#' @importFrom data.table fread
#' @importFrom GEOquery getGEO getGEOSuppFiles
#' @importFrom rPanglaoDB getSampleList getSampleComposition getSamples
#' @importFrom curl curl_fetch_memory
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \donttest{
#' # check all databases
#' CheckAPI()
#' }
CheckAPI <- function(database = c("SRA/ENA", "GEO", "PanglanDB", "UCSC Cell Browser", "Zenodo", "CELLxGENE", "Human Cell Atlas")) {
  if ("SRA/ENA" %in% database) {
    # check SRR accession
    srr.acce <- tryCatch(
      {
        GEOfastq::crawl_gsms("GSM5628756", max.workers = 1)
      },
      error = function(e) {
        warning("Error occurred when accessing the SRA accession via GEOfastq.")
        NULL
      }
    )
    if (!is.null(srr.acce)) {
      message("The API to access the SRA accession is OK!")
      # check ena link
      ena.search.url <- "https://www.ebi.ac.uk/ena/portal/api/search?query=run_accession="
      sra.url <- paste0(ena.search.url, srr.acce$run[1], "&result=read_run&fields=run_accession,sra_ftp,sra_md5,sra_bytes")
      # read file info
      sra.text <- data.table::fread(sra.url, showProgress = F)
      if (nrow(sra.text) == 0) {
        warning("Error occurred when accessing ENA download links.")
      } else {
        message("The API to access ENA download links is OK!")
      }
    }
  }
  if ("GEO" %in% database) {
    message("start checking APIs to access GEO!")
    # check GEO object
    GEO.obj <- tryCatch(
      {
        Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
        GEOquery::getGEO(GEO = "GSE302912")
      },
      error = function(e) {
        warning("Error occurred when accessing GEO object via GEOquery.")
        NULL
      }
    )
    if (!is.null(GEO.obj)) {
      message("The API to access the GEO object is OK!")
    }
    # check supplementary files
    GEO.supp <- tryCatch(
      {
        GEOquery::getGEOSuppFiles(GEO = "GSE302912", baseDir = tempdir())
      },
      error = function(e) {
        warning("Error occurred when accessing supplementary files via GEOquery.")
        NULL
      }
    )
    if (!is.null(GEO.supp)) {
      message("The API to access supplementary files is OK!")
    }
  }
  if ("PanglanDB" %in% database) {
    message("start checking APIs to access PanglanDB!")
    # check available metadata
    panglandb.meta <- tryCatch(
      {
        rPanglaoDB::getSampleList()
      },
      error = function(e) {
        warning("Error occurred when accessing metadata via rPanglaoDB.")
        NULL
      }
    )
    if (!is.null(panglandb.meta)) {
      message("The API to access all available samples is OK!")
    }
    # check cell type composition
    panglandb.ctc <- tryCatch(
      {
        ExtractPanglaoDBComposition(sra = "SRA429320", srs = "SRS1467249", local.data = FALSE)
      },
      error = function(e) {
        warning("Error occurred when accessing cell type composition via rPanglaoDB.")
        NULL
      }
    )
    if (!is.null(panglandb.ctc)) {
      message("The API to access cell type composition is OK!")
    }
    # check data download
    panglandb.seu <- tryCatch(
      {
        rPanglaoDB::getSamples(sra = "SRA429320", srs = "SRS1467249")
      },
      error = function(e) {
        warning("Error occurred when accessing available data via rPanglaoDB.")
        NULL
      }
    )
    if (!is.null(panglandb.seu)) {
      message("The API to access available data is OK!")
    }
  }
  if ("UCSC Cell Browser" %in% database) {
    message("start checking APIs to access UCSC Cell Browser!")
    # check available projects
    all.dataset.info <- tryCatch(
      {
        URLRetrieval("https://cells.ucsc.edu/dataset.json")
      },
      error = function(e) {
        warning("Error occurred when accessing API: https://cells.ucsc.edu/dataset.json")
        NULL
      }
    )
    if (!is.null(all.dataset.info)) {
      message("The API to access all available projects is OK!")
    }
    # check dataset detail
    dataset.info <- tryCatch(
      {
        URLRetrieval("https://cells.ucsc.edu/p-leidyi-adult/dataset.json")
      },
      error = function(e) {
        warning("Error occurred when accessing API: https://cells.ucsc.edu/p-leidyi-adult/dataset.json")
        NULL
      }
    )
    if (!is.null(dataset.info)) {
      message("The API to access detailed information of a given dataset is OK!")
    }
    # check available data
    dataset.data <- tryCatch(
      {
        URLRetrieval("https://cells.ucsc.edu/p-leidyi-adult/desc.json")
      },
      error = function(e) {
        warning("Error occurred when accessing API: https://cells.ucsc.edu/p-leidyi-adult/desc.json")
        NULL
      }
    )
    if (!is.null(dataset.data)) {
      message("The API to access the available data of a given dataset is OK!")
    }
    # check download link
    data.link <- "https://cells.ucsc.edu/p-leidyi-adult/matrix.mtx.gz"
    if (CheckURL(data.link)) {
      message("The API to access available files is OK!")
    } else {
      warning("Error occurred when accessing API: ", data.link)
    }
  }
  if ("Zenodo" %in% database) {
    message("start checking APIs to access Zenodo!")
    # check record detail
    record.api <- "https://zenodo.org/api/records/7243603"
    record.info <- URLRetrieval(record.api)
    if (length(names(record.info)) == 2 && all(names(record.info) == c("status", "message"))) {
      print(as.data.frame(record.info))
      warning("Error occurred when accessing API: ", record.api)
    } else {
      message("The API to access detailed information of a given doi is OK!")
      # check download link
      record.file.link <- record.info$files$links$self[1]
      if (CheckURL(record.file.link)) {
        message("The API to access available files is OK!")
      } else {
        warning("Error occurred when accessing API: ", record.file.link)
      }
    }
  }
  if ("CELLxGENE" %in% database) {
    message("start checking APIs to access CELLxGENE!")
    # check all collections
    collections.url <- "https://api.cellxgene.cziscience.com/curation/v1/collections/"
    collections.info <- URLRetrieval(collections.url)
    if (length(names(collections.info)) == 4 && all(names(collections.info) == c("detail", "status", "title", "type"))) {
      print(as.data.frame(collections.info))
      warning("Error occurred when accessing API: ", collections.url)
    } else {
      message("The API to access all available collections is OK!")
      # check collection detail
      sg.collections.url <- paste0(collections.url, collections.info$collection_id[1])
      sg.collections.info <- URLRetrieval(sg.collections.url)
      if (length(names(sg.collections.info)) == 4 && all(names(sg.collections.info) == c("detail", "status", "title", "type"))) {
        print(as.data.frame(sg.collections.info))
        warning("Error occurred when accessing API: ", sg.collections.url)
      } else {
        message("The API to access detailed information of a given collection is OK!")
        # check download link
        download.link <- sg.collections.info$datasets$assets[[1]]$url
        if (CheckURL(download.link)) {
          message("The API to access available files is OK!")
        } else {
          warning("Error occurred when accessing API: ", download.link)
        }
      }
    }
    # check dataset
    dataset.url <- "https://api.cellxgene.cziscience.com/curation/v1/datasets/5fde5c9c-5b1e-4df5-982a-3f8e7635161f/versions"
    dataset.content <- URLRetrieval(dataset.url)
    # check dataset_id
    if (length(names(dataset.content)) == 4 && all(names(dataset.content) == c("detail", "status", "title", "type"))) {
      print(as.data.frame(dataset.content))
      warning("Error occurred when accessing API: ", dataset.url)
    } else {
      "The API to access detailed information of a given dataset is OK!"
    }
  }
  if ("Human Cell Atlas" %in% database) {
    message("start checking APIs to access Human Cell Atlas!")
    # check available catalogs
    hca.catalogs <- tryCatch(
      {
        ExtractHCACatalogs()
      },
      error = function(e) {
        warning("Error occurred when accessing API: https://service.azul.data.humancellatlas.org/index/catalogs")
        NULL
      }
    )
    if (!is.null(hca.catalogs)) {
      message("The API to access all available catalogs is OK!")
      # check available projects
      project.url <- paste0("https://service.azul.data.humancellatlas.org/index/projects?catalog=", hca.catalogs[1], "&size=75")
      project.info <- URLRetrieval(project.url)
      if (length(names(project.info)) == 2 && all(names(project.info) == c("Code", "Message"))) {
        print(as.data.frame(project.info))
        warning("Error occurred when accessing API: ", project.url)
      } else {
        message("The API to access all available projects is OK!")
        # check download link
        uuids <- sapply(project.info$hits$projects, function(x) {
          x.info <- x$contributedAnalyses %>% unlist()
          x.info[grepl(pattern = "uuid$", x = names(x.info))]
        }) %>%
          unlist() %>%
          unique()
        if (length(uuids) > 0) {
          uuid.url <- paste0("https://service.azul.data.humancellatlas.org/repository/files/", uuids[1])
          if (CheckURL(uuid.url)) {
            message("The API to access available files is OK!")
          } else {
            warning("Error occurred when accessing API: ", uuid.url)
          }
        } else {
          warning("There are no uuids in the first 75 projects, so it is uncertain whether the API is OK.")
        }
      }
    }
  }
}
