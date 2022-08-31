FROM bioconductor/bioconductor_docker:RELEASE_3_15

LABEL \
    author = "David Currie" \
    description = "Temporary Dockerfile for nf-core/spatialtranscriptomics" \
    maintainer = "david.currie21@imperial.ac.uk"

## update system libraries
RUN apt-get update

RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("BioBase=2.56.0",\
    "BiocGenerics=0.42.0", "GeomxTools=3.0.1", ask = F))'

# install renv & packages
RUN Rscript -e 'install.packages("renv=0.15.5")'
RUN Rscript -e 'install.packages("cowplot=1.1.1")'
RUN Rscript -e 'install.packages("dplyr=1.0.9")'
RUN Rscript -e 'install.packages("GGally=2.1.2")'
RUN Rscript -e 'install.packages("ggforce=0.3.3")'
RUN Rscript -e 'install.packages("ggplot2=3.3.6")'
RUN Rscript -e 'install.packages("ggrepel=0.9.1")'
RUN Rscript -e 'install.packages("knitr=1.39")'
RUN Rscript -e 'install.packages("pheatmap=1.0.12")'
RUN Rscript -e 'install.packages("reshape2=1.4.4")'
RUN Rscript -e 'install.packages("Rtsne=0.16")'
RUN Rscript -e 'install.packages("scales=1.2.0")'
RUN Rscript -e 'install.packages("tidyr=1.2.0")'
RUN Rscript -e 'install.packages("umap=0.2.8.0")'