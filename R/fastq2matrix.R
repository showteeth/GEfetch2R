#' Run CellRanger on Downloaded FASTQ Files.
#'
#' @param sample.dir Directory contains all samples.
#' @param ref Path of folder containing 10x-compatible
#' transcriptome reference (\code{--transcriptome}).
#' @param localcores Set max cores the pipeline may request at one time (\code{--localcores}).
#' Only applies to local jobs (\code{--jobmode=local}). Default: 4.
#' @param localmem Set max GB the pipeline may request at one time (\code{--localmem}).
#' Only applies to local jobs (\code{--jobmode=local}). Default: 16.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param cr.path Path to cellranger. Default: NULL (conduct automatic detection).
#' @param cr.paras Parameters for \code{cellranger}.
#' Default: "--chemistry=auto --jobmode=local".
#'
#' @return Vector contains failed samples or NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' RunCellRanger(
#'   sample.dir = "/path/to/fastq", # GSMXXXX
#'   ref = "/path/to/cellranger_tiny_ref/3.0.0",
#'   out.folder = "/path/to/results",
#'   cr.path = "/path/to/cellranger-x.x.x/cellranger"
#' )
#' }
RunCellRanger <- function(sample.dir, ref, localcores = 4, localmem = 16, out.folder = NULL,
                          cr.path = NULL, cr.paras = "--chemistry=auto --jobmode=local") {
  # check valid sample dir
  dir.ext.flag <- dir.exists(sample.dir)
  if (!all(dir.ext.flag)) {
    old.sample.dir <- sample.dir
    sample.dir <- sample.dir[dir.ext.flag]
    if (length(sample.dir) == 0) {
      stop("The provided sample.dir: ", paste0(old.sample.dir, collapse = ", "), " doesn't exist, please check and re-run!")
    } else {
      message("Valid sample.dir: ", paste(sample.dir, collapse = ", "), ". Invalid sample.dir: ", paste(setdiff(old.sample.dir, sample.dir), collapse = ", "))
    }
  }
  # sample.fq.dir <- dir(path = sample.dir, full.names = TRUE)
  all.sample.cr <- sapply(sample.dir, function(x) {
    RunCellRangerSingle(
      fq.dir = x, transcriptome = ref, localcores = localcores,
      localmem = localmem, out.folder = out.folder,
      cr.path = cr.path, cr.paras = cr.paras
    )
  })
  # select fail samples
  fail.flag <- sapply(names(all.sample.cr), function(x) {
    !is.null(all.sample.cr[[x]])
  })
  fail.sample <- all.sample.cr[fail.flag]
  if (length(fail.sample) > 0) {
    message(
      "CellRanger failed on: ",
      paste0(fail.sample, collapse = ","), ". Please check and re-run!"
    )
    return(fail.sample)
  } else {
    message("CellRanger run successfully on all samples!")
    return(NULL)
  }
}

# run CellRanger on single sample
RunCellRangerSingle <- function(fq.dir, transcriptome, localcores = 4, localmem = 16, out.folder = NULL,
                                cr.path = NULL, cr.paras = "--chemistry=auto --jobmode=local") {
  # get cellranger path
  if (is.null(cr.path)) {
    # specify cellranger path
    cr.path <- Sys.which("cellranger")
    if (cr.path == "") {
      stop("Can not find cellranger automatically, please specify the path!")
    }
  } else {
    cr.path <- cr.path
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (!dir.exists(out.folder)) {
    message(out.folder, " does not exist, create automatically!")
    dir.create(out.folder, recursive = TRUE)
  }
  # check the path
  if (!dir.exists(fq.dir)) {
    message(fq.dir, " doesn't exist, please check and re-run!")
    return(basename(fq.dir))
  } else {
    fq.files <- list.files(path = fq.dir, pattern = "fastq.gz$", recursive = T)
    if (length(fq.files) == 0) {
      message("There is no fastq.gz files under ", fq.dir, " , please check and re-run!")
      return(basename(fq.dir))
    } else {
      fq1.file <- grep(pattern = ".*_S[0-9]_L[0-9]{3}_R1_[0-9]{3}.fastq.gz", x = fq.files, value = T)
      fq2.file <- grep(pattern = ".*_S[0-9]_L[0-9]{3}_R2_[0-9]{3}.fastq.gz", x = fq.files, value = T)
      if (length(fq1.file) > 0 && length(fq2.file) > 0) {
        if (length(fq1.file) != length(fq2.file)) {
          message("The number of R1 and R2 fastq.gz files under ", fq.dir, " is differ , please check and re-run!")
          return(basename(fq.dir))
        } else {
          sample.id <- basename(fq.dir)
          # merge multiple run as a single sample
          sample.name <- paste0(dir(path = fq.dir), collapse = ",")
          # sample.name <- basename(fq.dir)
          # check additional paras
          if (grepl(pattern = "--id|--transcriptome|--fastqs|--sample|--localcores|--localmem", cr.paras)) {
            message(
              cr.paras, " overlaps with built-in paras (--id, --transcriptome, --fastqs, --sample,",
              " --localcores, --localmem), please check and re-run!"
            )
          }
          # cellranger command
          cr.cmd <- paste0(
            "cd ", out.folder, " && ", cr.path, " count --id=", sample.id, " --transcriptome=", transcriptome,
            " --fastqs=", fq.dir, " --sample=", sample.name, " --localcores=", localcores,
            " --localmem=", localmem, " ", cr.paras
          )
          # run command
          message(paste("Calling CellRanger:", cr.cmd))
          cr.status <- system(cr.cmd, intern = TRUE)
          cr.status.code <- attr(cr.status, "status")
          if (!is.null(cr.status.code)) {
            cr.status.msg <- paste(cr.status, collapse = " ")
            warning(
              "Run CellRanger error on: ", sample.name, ". Error message :",
              cr.status.msg, " .Please check and re-run!"
            )
            return(sample.name)
          } else {
            message("Finish CellRanger: ", sample.name)
            return(NULL)
          }
        }
      } else {
        message("There is no R1 or R2 fastq.gz files under ", fq.dir, " , please check and re-run!")
        return(basename(fq.dir))
      }
    }
  }
}

#' Run STAR on Downloaded FASTQ Files.
#'
#' @param sample.dir Directory contains all samples.
#' @param ref Path of folder containing STAR version-compatible reference.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param thread The number of threads to use. Default: 4.
#' @param star.path Path to \code{STAR}. Default: NULL (conduct automatic detection).
#' @param star.paras Parameters for \code{STAR}.
#' Default: "--outBAMsortingThreadN 4 --twopassMode None".
#'
#' @return Vector contains failed samples or NULL.
#' @export
#'
#' @examples
#' \dontrun{
#' RunSTAR(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/star/index",
#'   out.folder = "/path/to/star/mapping",
#'   star.path = "/path/to/bin/STAR"
#' )
#' }
RunSTAR <- function(sample.dir, ref, out.folder = NULL, thread = 4, star.path = NULL,
                    star.paras = "--outBAMsortingThreadN 4 --twopassMode None") {
  # check valid sample dir
  dir.ext.flag <- dir.exists(sample.dir)
  if (!all(dir.ext.flag)) {
    old.sample.dir <- sample.dir
    sample.dir <- sample.dir[dir.ext.flag]
    if (length(sample.dir) == 0) {
      stop("The provided sample.dir: ", paste0(old.sample.dir, collapse = ", "), " doesn't exist, please check and re-run!")
    } else {
      message("Valid sample.dir: ", paste(sample.dir, collapse = ", "), ". Invalid sample.dir: ", paste(setdiff(old.sample.dir, sample.dir), collapse = ", "))
    }
  }
  # sample.fq.dir <- dir(path = sample.dir, full.names = TRUE)
  all.sample.star <- sapply(sample.dir, function(x) {
    RunSTARSingle(
      fq.dir = x, ref = ref, out.folder = out.folder, thread = thread,
      star.path = star.path, star.paras = star.paras
    )
  })
  # select fail samples
  fail.flag <- sapply(names(all.sample.star), function(x) {
    !is.null(all.sample.star[[x]])
  })
  fail.sample <- all.sample.star[fail.flag]
  if (length(fail.sample) > 0) {
    message(
      "STAR failed on: ",
      paste0(fail.sample, collapse = ","), ". Please check and re-run!"
    )
    return(fail.sample)
  } else {
    message("STAR run successfully on all samples!")
    return(NULL)
  }
}

# run STAR on single sample
RunSTARSingle <- function(fq.dir, ref, out.folder = NULL, thread = 4, star.path = NULL,
                          star.paras = "--outBAMsortingThreadN 4 --twopassMode None") {
  # get STAR path
  if (is.null(star.path)) {
    # specify STAR path
    star.path <- Sys.which("STAR")
    if (star.path == "") {
      stop("Can not find STAR automatically, please specify the path!")
    }
  } else {
    star.path <- star.path
  }
  # prepare output folder
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  if (!dir.exists(out.folder)) {
    message(out.folder, " does not exist, create automatically!")
    dir.create(out.folder, recursive = TRUE)
  }
  # check the path
  if (!dir.exists(fq.dir)) {
    message(fq.dir, " doesn't exist, please check and re-run!")
    return(basename(fq.dir))
  } else {
    # fq.files <- list.files(path = fq.dir, pattern = "fastq.gz$")
    fq.files <- list.files(path = fq.dir, pattern = "fastq.gz$", recursive = TRUE, full.names = TRUE)
    if (length(fq.files) == 0) {
      message("There is no fastq.gz files under ", fq.dir, " , please check and re-run!")
      return(basename(fq.dir))
    } else {
      # GSMXXXX
      sample.name <- basename(fq.dir)
      out.folder <- file.path(out.folder, sample.name, "")
      # run: SRRXXXX
      sub.samples <- dir(path = fq.dir)
      r1.pattern <- paste(paste0(sub.samples, "_1.fastq.gz"), collapse = "|")
      r2.pattern <- paste(paste0(sub.samples, "_2.fastq.gz"), collapse = "|")
      single.pattern <- paste(paste0(sub.samples, ".fastq.gz"), collapse = "|")
      pair.r1 <- grep(pattern = r1.pattern, x = fq.files, value = T)
      pair.r2 <- grep(pattern = r2.pattern, x = fq.files, value = T)
      single.read <- grep(pattern = single.pattern, x = fq.files, value = T)
      # check additional paras
      if (grepl(pattern = "--runThreadN|--genomeDir|--readFilesIn|--outSAMtype|--readFilesCommand|--quantMode|--outFileNamePrefix", star.paras)) {
        message(
          star.paras, " overlaps with built-in paras (--runThreadN, --genomeDir, --readFilesIn, --outSAMtype,",
          " --readFilesCommand, --quantMode, --outFileNamePrefix), please check and re-run!"
        )
      }
      # prepare command
      if (length(pair.r1) == length(pair.r2)) {
        if (length(pair.r1) > 1) {
          message(sample.name, " has multiple runs: ", paste(sub.samples, collapse = ", "))
          pair.r1 <- paste(pair.r1, collapse = ",")
          pair.r2 <- paste(pair.r2, collapse = ",")
        } else {
          message(sample.name, " has single run: ", paste(sub.samples, collapse = ", "))
        }
        # STAR command
        star.cmd <- paste(
          star.path, "--runThreadN", thread, "--genomeDir", ref,
          "--readFilesIn", pair.r1, pair.r2,
          "--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode GeneCounts",
          "--outFileNamePrefix", out.folder, star.paras, "&& mv",
          paste0(out.folder, "ReadsPerGene.out.tab"), paste0(out.folder, sample.name, ".txt")
        )
      } else if (length(single.read) > 0) {
        if (length(single.read) > 1) {
          message(sample.name, " has multiple runs: ", paste(sub.samples, collapse = ", "))
          single.read <- paste(single.read, collapse = ",")
        } else {
          message(sample.name, " has single run: ", paste(sub.samples, collapse = ", "))
        }
        # STAR command
        star.cmd <- paste(
          star.path, "--runThreadN", thread, "--genomeDir", ref,
          "--readFilesIn", single.read,
          "--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --quantMode GeneCounts",
          "--outFileNamePrefix", out.folder, star.paras, "&& mv",
          paste0(out.folder, "ReadsPerGene.out.tab"), paste0(out.folder, sample.name, ".txt")
        )
      } else {
        message(
          "There is no valid fastq.gz files under ", fq.dir,
          ". For pair-end,the fastq files should be: ", paste0(sample.name, sub.samples, "_[12].fastq.gz"),
          ". For single-end, the fastq file should be: ", paste0(sample.name, sub.samples, ".fastq.gz"),
          ". Please check and re-run!"
        )
        return(sample.name)
      }
      # run command
      message(paste("Calling STAR:", star.cmd))
      star.status <- system(star.cmd, intern = TRUE)
      star.status.code <- attr(star.status, "status")
      if (!is.null(star.status.code)) {
        star.status.msg <- paste(star.status, collapse = " ")
        warning(
          "Run STAR error on: ", sample.name, ". Error message :",
          star.status.msg, " .Please check and re-run!"
        )
        return(sample.name)
      } else {
        message("Finish STAR: ", sample.name)
        return(NULL)
      }
    }
  }
}

#' Pipe FASTQ files to SeuratObject and DESeqDataSet.
#'
#' @param sample.dir Directory contains all samples.
#' @param ref Path of folder containing 10x-compatible transcriptome \code{\link{RunCellRanger}}
#' STAR \code{\link{RunSTAR}} reference.
#' @param method Mapping methods, choose from CellRanger (10x Genomics) and STAR (Smart-seq2 or bulk RNA-seq).
#' Default: CellRanger.
#' @param localcores The max cores/thread used \code{\link{RunCellRanger}}/\code{\link{RunSTAR}}. Default: 4.
#' @param localmem The max memory (GB) used \code{\link{RunCellRanger}}. Default: 16.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param st.path Path to \code{STAR} or \code{cellranger}. Default: NULL (conduct automatic detection).
#' @param st.paras Parameters for \code{STAR} or \code{cellranger}.
#' Default: "--chemistry=auto --jobmode=local".
#' @param merge Logical, whether to merge the SeuratObjects, use when \code{method} is CellRanger. Default: TRUE.
#' @param count.col Column contains used count data (2: unstranded; 3: \code{stranded=yes}; 4: \code{stranded=reverse}),
#' use when \code{method} is STAR. Default: 2.
#' @param meta.data Dataframe contains sample information for DESeqDataSet, use when \code{method} is STAR. Default: NULL.
#' @param fmu Column of \code{meta.data} contains group information. Default: NULL.
#'
#' @return SeuratObject, DESeqDataSet or NULL.
#' @importFrom magrittr %>%
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom methods new
#' @importFrom data.table fread
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @export
#'
#' @examples
#' \dontrun{
#' # run CellRanger (10x Genomics)
#' # the sample.dir corresponding to sra.folder (SplitSRA) or out.folder (DownloadFastq)
#' seu <- Fastq2R(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/10x/ref",
#'   method = "CellRanger",
#'   out.folder = "/path/to/results",
#'   st.path = "/path/to/cellranger"
#' )
#' # run STAR (Smart-seq2 or bulk RNA-seq)
#' # the sample.dir corresponding to sra.folder/GSMXXXX (SplitSRA) or out.folder/GSMXXXX (DownloadFastq)
#' deobj <- Fastq2R(
#'   sample.dir = "/path/to/fastq",
#'   ref = "/path/to/star/ref",
#'   method = "STAR",
#'   out.folder = "/path/to/results",
#'   st.path = "/path/to/STAR",
#'   st.paras = "--outBAMsortingThreadN 4 --twopassMode None"
#' )
#' }
Fastq2R <- function(sample.dir, ref, method = c("CellRanger", "STAR"), localcores = 4, localmem = 16,
                    out.folder = NULL, st.path = NULL, st.paras = "--chemistry=auto --jobmode=local",
                    merge = TRUE, count.col = 2, meta.data = NULL, fmu = NULL) {
  # check parameters
  method <- match.arg(arg = method)
  if (method == "CellRanger") {
    res <- RunCellRanger(
      sample.dir = sample.dir, ref = ref, localcores = localcores,
      localmem = localmem, out.folder = out.folder,
      cr.path = st.path, cr.paras = st.paras
    )
    if (is.null(res)) {
      all.cr.folder <- dir(out.folder, full.names = TRUE)
      all.samples.folder <- file.path(all.cr.folder, "outs", "filtered_feature_bc_matrix")
      # check file
      valid.samples.folder <- Check10XFiles(folders = all.samples.folder, gene2feature = TRUE)
      if (length(valid.samples.folder) == 0) {
        stop("No valid sample folder detected under ", out.folder, ". Please check!")
      }
      # load to seurat
      seu.list <- sapply(valid.samples.folder, function(x) {
        x.mat <- Seurat::Read10X(data.dir = x)
        seu.obj <- Seurat::CreateSeuratObject(counts = x.mat, project = basename(x))
        seu.obj
      })
      # merge SeuratObject
      if (isTRUE(merge)) {
        out.obj <- mergeExperiments(seu.list)
      } else {
        out.obj <- seu.list
      }
      return(out.obj)
    } else {
      message("Some samples failed to run, skipping loading into Seurat!")
      return(NULL)
    }
  } else if (method == "STAR") {
    res <- RunSTAR(
      sample.dir = sample.dir, ref = ref, out.folder = out.folder, thread = localcores,
      star.path = st.path, star.paras = st.paras
    )
    if (is.null(res)) {
      all.txt <- list.files(path = out.folder, pattern = ".txt$", recursive = TRUE, full.names = TRUE)
      # read files
      count.list <- lapply(
        all.txt,
        function(x) {
          sample.count <- data.table::fread(file = x, select = c(1, count.col)) %>% as.data.frame()
          sample.count <- sample.count[!sample.count[[1]] %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"), ]
          colnames(sample.count) <- c("GeneName", gsub(pattern = "(.*).txt", replacement = "\\1", x = basename(x)))
          sample.count
        }
      )
      # create count matrix
      count.mat <- Reduce(f = function(x, y) {
        merge.mat <- merge(x, y, by = "GeneName", all = T)
      }, x = count.list)
      rownames(count.mat) <- count.mat$GeneName
      count.mat$GeneName <- NULL
      # loadding into DESeq2
      de.obj <- Loading2DESeq2(mat = count.mat, meta = meta.data, fmu = fmu)
      return(de.obj)
    } else {
      message("Some samples failed to run, skipping loading into DESeq2!")
      return(NULL)
    }
  }
}

#' Extract Run, Distinguish RNA-seq Type, Download Fastq, Perform Read Mapping, Load Output to R.
#'
#' @param gsm GSM number. Default: NULL (use \code{acce}).
#' @param acce GEO accession number. Default: NULL (use \code{gsm}).
#' \code{acce} and \code{gsm} cannot be both NULL.
#' @param force.type Force the RNA-seq type, used when failing to automatically identify the RNA-seq type. Available value: "10x", "Smart-seq2", "bulk".
#' If not NULL, skip automatic identification of RNA-seq type. Default: NULL.
#' @param star.ref Path of folder containing STAR reference,
#' used when bulk RNA-seq or Smart-seq2 scRNA-seq/mini-bulk RNA-seq. Default: NULL.
#' @param cellranger.ref Path of folder containing 10x-compatible transcriptome reference,
#' used when 10x Genomics scRNA-seq. Default: NULL.
#' @param out.folder Output folder. Default: NULL (current working directory).
#' @param timeout Maximum request time. Default: 36000000.
#' @param star.path Path to \code{STAR}. Default: NULL (conduct automatic detection).
#' @param cellranger.path Path to cellranger. Default: NULL (conduct automatic detection).
#' @param download.method Method to download fastq files, chosen from "download.file", "ascp" and "wget". Default: "wget".
#' @param ascp.path Path to ascp (/path/bin/ascp), please ensure that the relative path of asperaweb_id_dsa.openssh file
#' (/path/bin/ascp/../etc/asperaweb_id_dsa.openssh). Default: NULL (conduct automatic detection).
#' @param wget.path Path to wget. Default: NULL (conduct automatic detection).
#' @param star.paras Parameters for \code{STAR}.
#' Default: "--outBAMsortingThreadN 4 --twopassMode None".
#' @param cellranger.paras Parameters for \code{cellranger}.
#' Default: "--chemistry=auto --jobmode=local".
#' @param localcores Number of cores used, same as \code{localcores} for \code{cellranger} and \code{runThreadN} for \code{STAR}.
#' Default: 4.
#' @param localmem Set max GB the pipeline may request at one time. Default: 16.
#' @param count.col Column contains used count data (2: unstranded; 3: \code{stranded=yes}; 4: \code{stranded=reverse}),
#' use when bulk RNA-seq or Smart-seq2 scRNA-seq/mini-bulk RNA-seq. Default: 2.
#'
#' @return List of R objects.
#' @importFrom magrittr %>%
#' @importFrom GEOquery getGEO
#' @importFrom Biobase annotation experimentData pData phenoData notes sampleNames exprs
#' @importFrom parallel detectCores mclapply
#' @importFrom GEOfastq crawl_gsms
#' @importFrom curl curl_fetch_memory
#' @importFrom data.table fread
#' @importFrom tidyr separate_rows
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom utils download.file read.table
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom methods new
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # only bulk/Smart-seq2 RNA-seq
#' GSE127942.list <- DownloadFastq2R(
#'   acce = "GSE127942", star.ref = "/path/to/ref", out.folder = "/path/to/output",
#'   star.path = "/path/to/STAR", timeout = 3600000
#' )
#' # 10x Genomics scRNA-seq
#' GSE282929.list <- DownloadFastq2R(
#'   acce = "GSE282929", cellranger.ref = "/path/to/cellranger/ref",
#'   cellranger.path = "/path/to/cellranger",
#'   out.folder = "/path/to/output", timeout = 3600000
#' )
#' # mixture of 10x Genomics scRNA-seq and bulk RNA-seq
#' GSE305141.list <- DownloadFastq2R(
#'   acce = "GSE305141", star.ref = "/path/to/star/ref", cellranger.ref = "/path/to/cellranger/ref",
#'   star.path = "/path/to/STAR", cellranger.path = "/path/to/cellranger",
#'   out.folder = "/path/to/output", timeout = 3600000
#' )
#' # given GSM number
#' GSE127942.list <- DownloadFastq2R(
#'   gsm = c("GSM3656922", "GSM3656923"), star.ref = "/path/to/ref", out.folder = "/path/to/output",
#'   star.path = "/path/to/STAR", timeout = 3600000
#' )
#' }
DownloadFastq2R <- function(gsm = NULL, acce = NULL, force.type = NULL, star.ref = NULL, cellranger.ref = NULL, out.folder = NULL,
                            timeout = 36000000, star.path = NULL, cellranger.path = NULL,
                            download.method = c("wget", "download.file", "ascp"), ascp.path = NULL, wget.path = NULL,
                            star.paras = "--outBAMsortingThreadN 4 --twopassMode None",
                            cellranger.paras = "--chemistry=auto --jobmode=local",
                            localcores = 4, localmem = 16, count.col = 2) {
  # check parameter
  download.method <- match.arg(arg = download.method)
  # check ref
  if (is.null(star.ref) && is.null(cellranger.ref)) {
    stop("The star.ref and cellranger.ref are NULL, please provide at least one valid value!")
  }
  message("Step1: extract runs with GEO accession number or GSM number.")
  run.df <- suppressMessages(ExtractRun(gsm = gsm, acce = acce, timeout = 36000))
  if (is.null(force.type)) {
    message("Step2: distinguish bulk RNA-seq, 10x Genomics scRNA-seq, and Smart-seq2 automatically.")
    run.type.list <- DistinguishRNA(geo.runs = run.df)
  } else {
    message("Step2: the RNA-seq type is manually specified!")
    valid.force.type <- intersect(force.type, c("10x", "Smart-seq2", "bulk"))
    if (length(valid.force.type) == 0) {
      stop("There is no valid force.type, available values: ", paste(c("10x", "Smart-seq2", "bulk"), collapse = ", "))
    } else {
      if (length(valid.force.type) > 1) {
        message("Multiple force.type detected, use the first: ", valid.force.type[1])
        valid.force.type <- valid.force.type[1]
      }
      if (valid.force.type == "10x") {
        run.type.list$bulk.rna <- data.frame()
        run.type.list$scrna.ss2 <- data.frame()
        run.type.list$scrna.10x <- run.df
      } else if (valid.force.type == "Smart-seq2") {
        run.type.list$bulk.rna <- data.frame()
        run.type.list$scrna.10x <- data.frame()
        run.type.list$scrna.ss2 <- run.df
      } else if (valid.force.type == "bulk") {
        run.type.list$scrna.ss2 <- data.frame()
        run.type.list$scrna.10x <- data.frame()
        run.type.list$bulk.rna <- run.df
      }
    }
  }
  # check ref
  if ((nrow(run.type.list$bulk.rna) > 0) || (nrow(run.type.list$scrna.ss2) > 0)) {
    if (is.null(star.ref)) {
      stop("Detected bulk RNA-seq or Smart-seq2 scRNA-seq/mini-bulk RNA-seq datasets/runs, please provide star.ref!")
    }
  }
  if (nrow(run.type.list$scrna.10x) > 0) {
    if (is.null(cellranger.ref)) {
      stop("Detected 10x Genomics scRNA-seq datasets/runs, please provide cellranger.ref!")
    }
  }
  # download and mapping
  if (is.null(out.folder)) {
    out.folder <- getwd()
  }
  # create fastq output folder
  fq.out.folder <- file.path(out.folder, "fastq")
  mapping.folder <- file.path(out.folder, "mapping")
  dir.create(path = fq.out.folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = mapping.folder, showWarnings = FALSE, recursive = TRUE)
  # output list
  out.li <- list()
  if (nrow(run.type.list$bulk.rna) > 0) {
    bulk.fq.folder <- file.path(fq.out.folder, "bulkRNAseq")
    message("Step3: Download ", nrow(run.type.list$bulk.rna), " bulk RNA-seq datasets/runs to: ", bulk.fq.folder)
    bulk.down <- DownloadFastq(
      gsm.df = run.type.list$bulk.rna, out.folder = bulk.fq.folder, download.method = download.method,
      ascp.path = ascp.path, quiet = TRUE, wget.path = wget.path,
      timeout = timeout, parallel = FALSE, use.cores = NULL,
      format.10x = FALSE, remove.raw = FALSE
    )
    if (!is.null(bulk.down)) {
      warning("Error occured when downloading bulk RNA-seq datasets/runs, skip subsequent mapping and loading processes.")
      bulk.deobj <- NULL
    } else {
      message("Step4: read mapping and load to R (DESeq2).")
      bulk.mapping.folder <- file.path(mapping.folder, "bulkRNAseq")
      bulk.gsm.folders <- unique(file.path(bulk.fq.folder, run.type.list$bulk.rna$gsm_name))
      bulk.deobj <- Fastq2R(
        sample.dir = bulk.gsm.folders, ref = star.ref, method = "STAR",
        out.folder = bulk.mapping.folder, st.path = star.path, st.paras = star.paras,
        localcores = localcores, count.col = count.col, meta.data = NULL, fmu = NULL
      )
    }
    out.li$bulk.rna <- bulk.deobj
  }
  if (nrow(run.type.list$scrna.10x) > 0) {
    fq.10x.folder <- file.path(fq.out.folder, "10x")
    message("Step3: Download ", nrow(run.type.list$scrna.10x), " 10x Genomics scRNA-seq datasets/runs to: ", fq.10x.folder)
    sc.10x.down <- DownloadFastq(
      gsm.df = run.type.list$scrna.10x, out.folder = fq.10x.folder, download.method = download.method,
      ascp.path = ascp.path, quiet = TRUE, wget.path = wget.path,
      timeout = timeout, parallel = FALSE, use.cores = NULL,
      format.10x = TRUE, remove.raw = TRUE
    )
    if (!is.null(sc.10x.down)) {
      warning("Error occured when downloading 10x Genomics scRNA-seq datasets/runs, skip subsequent mapping and loading processes.")
      sc.10x.down <- NULL
    } else {
      message("Step4: read mapping and load to R (Seurat).")
      sc.10x.mapping.folder <- file.path(mapping.folder, "10x")
      sc.10x.gsm.folders <- unique(file.path(fq.10x.folder, run.type.list$scrna.10x$gsm_name))
      sc.10x.seu <- Fastq2R(
        sample.dir = sc.10x.gsm.folders, ref = cellranger.ref, method = "CellRanger",
        out.folder = sc.10x.mapping.folder, st.path = cellranger.path, st.paras = cellranger.paras,
        localcores = localcores, localmem = localmem, merge = FALSE
      )
    }
    out.li$scrna.10x <- sc.10x.seu
  }
  if (nrow(run.type.list$scrna.ss2) > 0) {
    ss2.fq.folder <- file.path(fq.out.folder, "smartseq2")
    message("Step3: Download ", nrow(run.type.list$scrna.ss2), " Smart-seq2 scRNA-seq/mini-bulk RNA-seq datasets/runs to: ", ss2.fq.folder)
    ss2.down <- DownloadFastq(
      gsm.df = run.type.list$scrna.ss2, out.folder = ss2.fq.folder, download.method = download.method,
      ascp.path = ascp.path, quiet = TRUE, wget.path = wget.path,
      timeout = timeout, parallel = FALSE, use.cores = NULL,
      format.10x = FALSE, remove.raw = FALSE
    )
    if (!is.null(ss2.down)) {
      warning("Error occured when downloading Smart-seq2 scRNA-seq/mini-bulk RNA-seq datasets/runs, skip subsequent mapping and loading processes.")
      ss2.deobj <- NULL
    } else {
      message("Step4: read mapping and load to R (DESeq2).")
      ss2.mapping.folder <- file.path(mapping.folder, "smartseq2")
      ss2.gsm.folders <- unique(file.path(ss2.fq.folder, run.type.list$scrna.ss2$gsm_name))
      ss2.deobj <- Fastq2R(
        sample.dir = ss2.gsm.folders, ref = star.ref, method = "STAR",
        out.folder = ss2.mapping.folder, st.path = star.path, st.paras = star.paras,
        localcores = localcores, count.col = 2, meta.data = NULL, fmu = NULL
      )
    }
    out.li$scrna.ss2 <- ss2.deobj
  }
  if ((nrow(run.type.list$bulk.rna)) == 0 && (nrow(run.type.list$scrna.10x)) == 0 && (nrow(run.type.list$scrna.ss2)) == 0) {
    message("Nothing to be done (no bulk RNA-seq, Smart-seq2 scRNA-seq/mini-bulk RNA-seq, 10x Genomics scRNA-seq datasets/runs detected).")
  }
  return(out.li)
}
