FROM bioconductor/bioconductor_docker:devel

LABEL \
    author = "David Currie" \
    description = "Temporary Dockerfile for nf-core/spatialtranscriptomics" \
    maintainer = "david.currie21@imperial.ac.uk"

## update system libraries
RUN apt-get update

RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("Biobase", "BiocGenerics", "GeomxTools", ask = F))'

# install renv & packages
RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'install.packages("cowplot")'
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("GGally")'
RUN Rscript -e 'install.packages("ggforce")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("ggrepel")'
RUN Rscript -e 'install.packages("knitr")'
RUN Rscript -e 'install.packages("pheatmap")'
RUN Rscript -e 'install.packages("reshape2")'
RUN Rscript -e 'install.packages("rtsne")'
RUN Rscript -e 'install.packages("scales")'
RUN Rscript -e 'install.packages("umap")'