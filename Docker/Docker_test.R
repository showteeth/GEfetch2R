library(GEfetch2R)

# raw data ---------------
CheckAPI(database = c("GEO", "SRA/ENA"))

# extract metadata
GSE130636.runs <- ExtractRun(acce = "GSE130636", platform = "GPL20301")
# a small test
GSE130636.runs <- GSE130636.runs[GSE130636.runs$run %in% c("SRR9004346", "SRR9004351"), ]

# download with prefetch
GSE130636.down <- DownloadSRA(
  gsm.df = GSE130636.runs,
  prefetch.path = "/opt/sratoolkit.3.0.6-ubuntu64/bin/prefetch",
  out.folder = "/home/rstudio/download_sra/prefetch"
)

# download from ENA using download.file
GSE130636.down <- DownloadSRA(
  gsm.df = GSE130636.runs, download.method = "download.file",
  timeout = 360000, out.folder = "/home/rstudio/download_sra/download_file",
  rename = TRUE, parallel = TRUE, use.cores = 2
)

# download from ENA using ascp
GSE130636.down <- DownloadSRA(
  gsm.df = GSE130636.runs, download.method = "ascp",
  ascp.path = "/opt/conda/bin/ascp", max.rate = "300m",
  rename = TRUE, out.folder = "/home/rstudio/download_sra/ascp",
  parallel = TRUE, use.cores = 2
)

# download from ENA using wget
GSE130636.down <- DownloadSRA(
  gsm.df = GSE130636.runs, download.method = "wget",
  wget.path = "/usr/bin/wget", timeout = 360000, 
  out.folder = "/home/rstudio/download_sra/wget", rename = TRUE, 
  parallel = TRUE, use.cores = 2
)

# count matrix -----------
CheckAPI(database = c("GEO", "PanglaoDB", "UCSC Cell Browser"))

# * GEO --------
# set VROOM_CONNECTION_SIZE to avoid error: Error: The size of the connection buffer (786432) was not large enough
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
# extract metadata
GSE297431.meta <- ExtractGEOMeta(acce = "GSE297431")
GSE297431.meta[1:3, c("title", "geo_accession", "source_name_ch1", "description", "cell type")]
GSE297431.meta.supp <- ExtractGEOMeta(
  acce = "GSE297431", down.supp = TRUE,
  supp.idx = 2 # specify the index of used supplementary file
)
head(GSE297431.meta.supp)

# Smart-seq2
GSE94820.seu = ParseGEO(acce = "GSE94820", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)
GSE94820.seu
GSE241004.seu = ParseGEO(acce = "GSE241004", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)
GSE241004.seu
GSE182219.seu = ParseGEO(acce = "GSE182219", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)
GSE182219.seu
GSE195839.seu = ParseGEO(acce = "GSE195839", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)
GSE195839.seu

# 10x Genomics
GSE310026.seu = ParseGEO(acce = "GSE310026", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/rstudio/GEO", merge = FALSE, timeout = 36000)
GSE310026.seu
GSE270070.seu = ParseGEO(acce = "GSE270070", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/rstudio/GEO", merge = FALSE, timeout = 36000)
GSE270070.seu
GSE274844.seu = ParseGEO(acce = "GSE274844", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/rstudio/GEO", merge = FALSE, timeout = 36000)
GSE274844.seu
GSE293930.seu = ParseGEO(acce = "GSE293930", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/rstudio/GEO", merge = FALSE, timeout = 36000)
GSE293930.seu
GSE277137.seu = ParseGEO(acce = "GSE277137", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/rstudio/GEO", merge = FALSE, timeout = 36000)
GSE277137.seu

# * PanglaoDB -----------
# cell type composition
lung.composition <- ExtractPanglaoDBComposition(sra = "SRA570744")

# count matrix
lung.seu <- ParsePanglaoDB(sra = "SRA570744", srs = "SRS2253536")

StatDBAttribute(df = PanglaoDBMeta, filter = c("species", "protocol"), database = "PanglaoDB")
hsa.meta <- ExtractPanglaoDBMeta(
  species = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"),
  show.cell.type = TRUE, cell.num = c(1000, 2000)
)
hsa.composition <- ExtractPanglaoDBComposition(
  meta = hsa.meta
)
hsa.seu <- ParsePanglaoDB(hsa.meta[1:3, ], merge = TRUE)

# * UCSC ----------
# extract cell type composition
ut.sample.ct <- ExtractCBComposition(link = c(
  "https://cells.ucsc.edu/?ds=adult-ureter", # collection
  "https://cells.ucsc.edu/?ds=adult-testis" # dataset
))
ut.sample.ct[1:5, c("title", "CellType", "Num")]

# processed objects ----------
CheckAPI(database = c("GEO", "Zenodo", "CELLxGENE", "Human Cell Atlas"))

# * GEO ----------------
# rds
GSE286410.seu = ParseGEOProcessed(acce = "GSE286410", timeout = 360000, supp.idx = 3, out.folder = "/home/rstudio/GEO/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)
GSE286410.seu
# loom
ParseGEOProcessed(acce = "GSE254944", timeout = 360000, supp.idx = 1, out.folder = "/home/rstudio/GEO/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))
# RData
ParseGEOProcessed(acce = "GSE306745", timeout = 360000, supp.idx = 1, out.folder = "/home/rstudio/GEO/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))
GSE306745.list <- LoadRData(
  rdata = "/home/rstudio/GEO/Processed/GSE306745/GSE306745_Microglia2.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
# h5ad
ParseGEOProcessed(acce = "GSE290433", timeout = 360000, supp.idx = 2, out.folder = "/home/rstudio/GEO/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))

# * Zenodo -------------
# single doi
zebrafish.df <- ExtractZenodoMeta(doi = "10.5281/zenodo.7243603")
# vector dois
multi.dois <- ExtractZenodoMeta(doi = c("1111", "10.5281/zenodo.7243603", "10.5281/zenodo.7244441"))

# * cellxgene -------------
cellxgene.given.h5ad <- ParseCELLxGENE(
  link = c(
    "https://cellxgene.cziscience.com/collections/77f9d7e9-5675-49c3-abed-ce02f39eef1b", # collection
    "https://cellxgene.cziscience.com/e/e12eb8a9-5e8b-4b59-90c8-77d29a811c00.cxg/" # dataset
  ),
  timeout = 36000000,
  out.folder = "/home/rstudio/download_cellxgene"
)

# * HCA ---------------
hca.given.download <- ParseHCA(
  link = c(
    "https://explore.data.humancellatlas.org/projects/902dc043-7091-445c-9442-d72e163b9879",
    "https://explore.data.humancellatlas.org/projects/cdabcf0b-7602-4abf-9afb-3b410e545703"
  ), timeout = 36000000,
  out.folder = "/home/rstudio/download_hca"
)

# object conversion -------------
# * test data --------------------
library(Seurat)
library(scRNAseq)
pbmc_small
head(pbmc_small@meta.data)
seger <- scRNAseq::SegerstolpePancreasData()
seger

# * convert SeuratObject to other objects -------------------------
# SeuratObject to SingleCellExperiment
sce.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "SCE")
sce.obj
# SeuratObject to CellDataSet
cds.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", reduction = "tsne", to = "CellDataSet")
cds.obj
# SeuratObject to cell_data_set
cds3.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "cell_data_set")
cds3.obj
# SeuratObject to AnnData
Seu2AD(
  seu.obj = pbmc_small, method = "sceasy",
  out.folder = "/home/rstudio/conversion",
  assay = "RNA", slot = "counts", conda.path = "/opt/conda"
)
Seu2AD(seu.obj = pbmc_small, method = "SeuratDisk", out.folder = "/home/rstudio/conversion",
       assay = "RNA", save.scale = TRUE)
Seu2AD(seu.obj = pbmc_small, method = "scDIOR",
       out.folder = "/home/rstudio/conversion", assay = "RNA", save.scale = TRUE)
# SeuratObject to loom
ExportSeurat(
  seu.obj = pbmc_small, assay = "RNA", to = "loom",
  loom.file = "/home/rstudio/conversion/pbmc_small.loom"
)

# * convert other objects to SeuratObject ------------------------
# generate the pbmc3k.h5ad: https://github.com/showteeth/GEfetch2R/blob/main/man/benchmark/generate_pbmc3k_anndata.ipynb
# upload to the home directory
# or download: download.file("https://github.com/showteeth/GEfetch2R/raw/refs/heads/main/man/benchmark/pbmc3k.h5ad", destfile = "/home/rstudio/pbmc3k.h5ad")

# SingleCellExperiment to SeuratObject
seu.obj.sce <- ImportSeurat(
  obj = sce.obj, from = "SCE", count.assay = "counts",
  data.assay = "logcounts", assay = "RNA"
)
# CellDataSet to SeuratObject
seu.obj.cds <- ImportSeurat(
  obj = cds.obj, from = "CellDataSet",
  count.assay = "counts", assay = "RNA"
)
# cell_data_set to SeuratObject
seu.obj.cds3 <- ImportSeurat(
  obj = cds3.obj, from = "cell_data_set",
  count.assay = "counts", data.assay = "logcounts", assay = "RNA"
)
# AnnData to SeuratObject
ann.sceasy <- AD2Seu(
  anndata.file = "~/pbmc3k.h5ad", method = "sceasy",
  assay = "RNA", slot = "scale.data"
)
ann.seu <- AD2Seu(
  anndata.file = "~/pbmc3k.h5ad",
  method = "SeuratDisk", assay = "RNA", load.assays = c("RNA")
)
ann.scdior <- AD2Seu(
  anndata.file = "~/pbmc3k.h5ad",
  method = "scDIOR", assay = "RNA"
)
ann.schard <- AD2Seu(
  anndata.file = "~/pbmc3k.h5ad",
  method = "schard", assay = "RNA", use.raw = T
)
ann.seuscdior <- AD2Seu(
  anndata.file = "~/pbmc3k.h5ad",
  method = "SeuratDisk+scDIOR", assay = "RNA", load.assays = c("RNA")
)
# loom to SeuratObject
seu.obj.loom <- ImportSeurat(loom.file = "~/conversion/pbmc_small.loom", from = "loom")

# * conversion between SingleCellExperiment and AnnData ----------------
# SingleCellExperiment to AnnData
SCE2AD(
  sce.obj = seger, method = "zellkonverter",
  out.folder = "/home/rstudio/conversion", slot = "counts",
  conda.path = "/opt/conda"
)
SCE2AD(
  sce.obj = seger, method = "sceasy", out.folder = "/home/rstudio/conversion",
  slot = "counts", conda.path = "/opt/conda"
)
seger.scdior <- seger
library(SingleCellExperiment)
rowData(seger.scdior)$varm <- NULL
SCE2AD(sce.obj = seger.scdior, method = "scDIOR", out.folder = "/home/rstudio/conversion")

# AnnData to SingleCellExperiment
sce.zell <- AD2SCE(
  anndata.file = "/home/rstudio/pbmc3k.h5ad",
  method = "zellkonverter", slot = "scale.data",
  use.raw = TRUE, conda.path = "/opt/conda"
)
sce.scdior <- AD2SCE(
  anndata.file = "/home/rstudio/pbmc3k.h5ad",
  method = "scDIOR", assay = "RNA",
  use.raw = TRUE, conda.path = "/opt/conda"
)
sce.schard <- AD2SCE(
  anndata.file = "/home/rstudio/pbmc3k.h5ad",
  method = "schard", use.raw = TRUE
)

# * conversion between SingleCellExperiment and loom ------------------
# SingleCellExperiment to loom
SCELoom(
  from = "SingleCellExperiment", to = "loom", sce = seger,
  loom.file = "/home/rstudio/conversion/seger.loom"
)
# loom to SingleCellExperiment
seger.loom <- SCELoom(
  from = "loom", to = "SingleCellExperiment",
  loom.file = "/home/rstudio/conversion/seger.loom"
)




