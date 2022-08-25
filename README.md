# ![nf-core/spatialtranscriptomicsgeomx](docs/images/nf-core-spatialtranscriptomicsgeomx_logo_light.png#gh-light-mode-only) ![nf-core/spatialtranscriptomicsgeomx](docs/images/nf-core-spatialtranscriptomicsgeomx_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/spatialtranscriptomicsgeomx/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/spatialtranscriptomicsgeomx/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/spatialtranscriptomicsgeomx/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/spatialtranscriptomicsgeomx/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?logo=Amazon%20AWS)](https://nf-co.re/spatialtranscriptomicsgeomx/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/spatialtranscriptomicsgeomx)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23spatialtranscriptomicsgeomx-4A154B?logo=slack)](https://nfcore.slack.com/channels/spatialtranscriptomicsgeomx)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/spatialtranscriptomicsgeomx** is a bioinformatics pipeline to process and analyse spatial transcriptomics data from the GeoMx Digital Spatial Profiler.

This pipeline was built using R, [Nextflow](https://www.nextflow.io), and Docker/singularity containers.

## Pipeline summary

The pipeline performs the following processes:
0. Load the data
1. Overview of the data
2. Quality control and preprocessing
   2.1 Segement quality control
   2.2 Probe quality control
   2.3 Aggregate probe-level counts to gene-level
   2.4 Measuring the limit of quantification
   2.5 Filtering segments and genes with low detection rates
3. Normalisation
4. Unsupervised analysis
5. Differential expression
6. Visual differentially expressed genes

These steps map to script in [bin/](./bin/).

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.04.4`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/spatialtranscriptomicsgeomx -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. It is recommended that you use Docker to run this pipeline. You can download the [Docker image](https://hub.docker.com/r/currie32/st-geomx) from Docker Hub.

5. Start running your own analysis!

   ```console
   nextflow run main.nf -with-docker currie32/st-geomx
   ```

6. You can view the results from running the pipeline by opening the jupyter notebook: [summary_of_outputs.ipynb](./summary_of_outputs.ipynb)

## Documentation

The nf-core/spatialtranscriptomicsgeomx pipeline comes with documentation about the pipeline [usage](https://nf-co.re/spatialtranscriptomicsgeomx/usage), [parameters](https://nf-co.re/spatialtranscriptomicsgeomx/parameters) and [output](https://nf-co.re/spatialtranscriptomicsgeomx/output).

## Credits

nf-core/spatialtranscriptomicsgeomx was originally written by David Currie.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#spatialtranscriptomicsgeomx` channel](https://nfcore.slack.com/channels/spatialtranscriptomicsgeomx) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/spatialtranscriptomicsgeomx for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
