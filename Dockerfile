FROM rocker/binder

## USER
USER root

## Update
RUN apt-get update

# Change workdirectory
ENV HOME /home
WORKDIR ${HOME}

## Install Conda
RUN apt-get install -y wget bzip2
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
RUN rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH /home/miniconda3/bin:$PATH

## Update
RUN conda update conda
RUN conda update --all

## Add chanel
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

## Install snakemake
RUN conda install -c bioconda -c conda-forge snakemake

## Install fastqc, bowtie2, htseq and samtools
RUN conda install -y fastqc bowtie2 htseq samtools

## Install Bioconductor package
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") ; BiocManager::install("DESeq2", version = "3.8")'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") ; BiocManager::install("edgeR", version = "3.8")'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") ; BiocManager::install("genefilter", version = "3.8")'

# Install CRAN package
RUN Rscript -e "install.packages('devtools', repos='https://cran.rstudio.com/', dependencies = TRUE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN Rscript -e "install.packages(c('shinydashboard','DT', 'FactoMineR', 'corrplot','plotly'), repos='https://cran.rstudio.com/', dependencies = TRUE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN Rscript -e "install.packages(c('shinyWidgets','colourpicker'), repos='https://cran.rstudio.com/', dependencies = TRUE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds


# Install Sarstools
RUN Rscript -e 'library(devtools) ; install_github("PF2-pasteur-fr/SARTools", build_vignettes=TRUE)'

# Install fastq dump

RUN wget --quiet "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-ubuntu64.tar.gz"
RUN tar -xzf sratoolkit.2.9.2-ubuntu64.tar.gz
RUN rm sratoolkit.2.9.2-ubuntu64.tar.gz
ENV PATH /home/sratoolkit.2.9.2-ubuntu64/bin:$PATH

## USER
USER rstudio
ENV HOME /home/rstudio
WORKDIR ${HOME}

CMD jupyter notebook --ip 0.0.0.0 --NotebookApp.token='' --NotebookApp.password=''
