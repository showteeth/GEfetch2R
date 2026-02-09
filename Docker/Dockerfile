# based on ubuntu
FROM rocker/rstudio:4.4.3
LABEL maintainer="songyb0519@gmail.com"

# install system requirements
USER root
ADD sources.list /etc/apt/
# Ubuntu 24.04
ADD ubuntu.sources /etc/apt/sources.list.d/
RUN apt-get update
RUN apt-get install -y -q libfontconfig1-dev libz-dev
RUN apt-get install -y -q pkg-config libxml2 libxml2-dev libfreetype-dev
RUN apt-get install -y -q libpng-dev libuuid1 libtiff-dev cmake
RUN apt-get install -y -q libcairo2-dev libxt-dev libhdf5-dev libbz2-dev
RUN apt-get install -y -q tk libgdal-dev libharfbuzz-dev libfribidi-dev libudunits2-dev
RUN apt-get install -y -q patch libglpk40
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

# set mirror
## CRAN package mirror
RUN echo "options('repos' = c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/'))" >> /usr/local/lib/R/etc/Rprofile.site
RUN echo "options('timeout' = 9999999)" >> /usr/local/lib/R/etc/Rprofile.site
COPY Renviron /etc/R/Renviron.site
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_version('BiocManager', upgrade = 'never')"
RUN R -e "devtools::install_version('SeuratObject', version = '4.1.4', upgrade = 'never')"
RUN R -e "devtools::install_version('Seurat', version = '4.4.0', upgrade = 'never')"
RUN R -e "devtools::install_github('alexvpickering/GEOfastq@572589716fd2747658eb80c70e523eed6c5688b6', upgrade = 'never')"
RUN R -e "devtools::install_github('cellgeni/sceasy@0cfc0e39da4b3ce4abf19d2171daa0e4d2acdd03', upgrade = 'never')"
RUN R -e "devtools::install_github('mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c', upgrade = 'never')"
RUN R -e "devtools::install_version('R.utils', upgrade = 'never')"
COPY Renviron ~/.Renviron
RUN R -e "devtools::install_github('satijalab/seurat-wrappers@b6d519b69e8a6364f7275d116aa8b0c1ee783364', upgrade = 'never')"
RUN R -e "BiocManager::install('SingleCellExperiment', update = FALSE)"
RUN R -e "devtools::install_github('theislab/zellkonverter@7b118653a471330b3734dcfee60c3537352ecb8d', upgrade = 'never')"
RUN R -e "devtools::install_version('blob', upgrade = 'never')"
RUN R -e "BiocManager::install('scRNAseq', update = FALSE)"
RUN R -e "devtools::install_version('grr', version = '0.9.5', upgrade = 'never')"
RUN R -e "devtools::install_version('Matrix.utils', version = '0.9.8', upgrade = 'never')"
RUN R -e "devtools::install_github('cole-trapnell-lab/leidenbase@430f22af6982cc7d2e6e77f6b0df47bc970dcbce', upgrade = 'never')"
RUN R -e "devtools::install_github('cole-trapnell-lab/monocle3@87f6e88760566065179ae5039d524da7d032cc15', upgrade = 'never')"
RUN R -e "devtools::install_version('qlcMatrix', version = '0.9.7', upgrade = 'never')"
RUN R -e "BiocManager::install('monocle', update = FALSE)"
RUN R -e "devtools::install_bioc('3.12/GEOquery', upgrade = 'never')"
RUN R -e "BiocManager::install('scater', update = FALSE)"
RUN R -e "BiocManager::install('LoomExperiment', update = FALSE)"
RUN R -e "BiocManager::install('DESeq2', update = FALSE)"
RUN R -e "BiocManager::install('edgeR', update = FALSE)"
RUN R -e "devtools::install_github('cellgeni/schard@45a869885887f5ac01d7eff839a7b2b923b06c2e', upgrade = 'never')"
RUN R -e "devtools::install_github('JiekaiLab/dior@2b1ea47b6661c8a10d9455f3baeeccb8f12be2f0', upgrade = 'never')"
RUN R -e "install.packages('tiledbsoma', repos = c('https://tiledb-inc.r-universe.dev', 'https://cloud.r-project.org'))"
RUN R -e "install.packages('cellxgene.census', repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org'))"
RUN R -e "devtools::install_version('rPanglaoDB', version = '0.2.1', upgrade = 'never')"
RUN R -e "devtools::install_version('readxl', version = '1.4.3', upgrade = 'never')"
RUN R -e "BiocManager::install('BiocStyle', update = FALSE)"

# create python environment
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -O /tmp/miniconda.sh \
	&& /bin/bash /tmp/miniconda.sh -bp /opt/conda \
	&& rm -rf /tmp/miniconda.sh
# install Python packages
RUN /opt/conda/bin/conda install -y python=3.7.3 && \
	/opt/conda/bin/pip install scanpy==1.9.3 && \
	/opt/conda/bin/pip install loompy
RUN /opt/conda/bin/conda install -c bioconda -y samtools && \
	/opt/conda/bin/pip install diopy
RUN /opt/conda/bin/conda install -c hcc -y aspera-cli=3.9.6
RUN /opt/conda/bin/conda install -c conda-forge -y ncurses
ENV PATH="${PATH}:/opt:/opt/conda/bin"

# command-line tools
RUN cd /opt && wget --quiet https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.6/sratoolkit.3.0.6-ubuntu64.tar.gz -O sratoolkit.3.0.6-ubuntu64.tar.gz && \
	tar -zxf sratoolkit.3.0.6-ubuntu64.tar.gz 
ENV PATH="${PATH}:/opt/sratoolkit.3.0.6-ubuntu64/bin"
RUN /opt/conda/bin/pip install --no-cache-dir -U git+https://github.com/rvalieris/parallel-fastq-dump.git && \
	/opt/conda/bin/conda clean --all --yes
RUN cd /opt && wget --quiet https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux -O bamtofastq_linux && chmod +x bamtofastq_linux

# add codes for large-scale evaluation (GEO)
RUN mkdir -p /home/rstudio/GEOBench
ADD GEO_evaluation_download.R /home/rstudio/GEOBench/
ADD GEO_evaluation_createObject.R /home/rstudio/GEOBench/
# add codes for Docker usability test
ADD Docker_test.R /home/rstudio/
ADD pbmc3k.h5ad /home/rstudio/

# install GEfetch2R
RUN R -e "devtools::install_github('showteeth/GEfetch2R', upgrade = 'never')"