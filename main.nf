#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/spatialtranscriptomicsgeomx
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/spatialtranscriptomicsgeomx
    Website: https://nf-co.re/spatialtranscriptomicsgeomx
    Slack  : https://nfcore.slack.com/channels/spatialtranscriptomicsgeomx
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.path_to_data = "$baseDir"
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = 'results'

workflow {

    image_load_data = load_data(params.path_to_data)
    image_sample_overview = sample_overview(params.path_to_data, image_load_data)
    image_segment_qc = segment_qc(params.path_to_data, image_sample_overview)
    image_probe_qc = probe_qc(params.path_to_data, image_segment_qc)
    image_aggregate_counts = aggregate_counts(params.path_to_data, image_probe_qc)
    image_limit_of_quantification = limit_of_quantification(params.path_to_data, image_aggregate_counts)
    image_filtering = filtering(params.path_to_data, image_limit_of_quantification)
    image_normalisation = normalisation(params.path_to_data, image_filtering)
    image_unsupervised_analysis = unsupervised_analysis(params.path_to_data, image_normalisation)
    image_differential_expression = differential_expression(params.path_to_data, image_unsupervised_analysis)
    visualising_de_genes(params.path_to_data, image_differential_expression)
}

process load_data {
    input:
    path path_to_data

    output:
    stdout emit: image_load_data

    script:
    """
    2_load_data.R $path_to_data
    """
}

process sample_overview {
    input:
    path path_to_data
    file image_load_data

    output:
    stdout emit: image_sample_overview

    script:
    """
    3_sample_overview.R $path_to_data
    """
}

process segment_qc {
    input:
    path path_to_data
    file image_sample_overview

    output:
    stdout emit: image_sample_overview

    script:
    """
    4_1_segment_qc.R $path_to_data
    """
}

process probe_qc {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    4_2_probe_qc.R
    """
}

process aggregate_counts {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    4_3_aggregate_counts.R
    """
}

process limit_of_quantification {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    4_4_limit_of_quantification.R
    """
}

process filtering {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    4_5_filtering.R
    """
}

process normalisation {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    5_normalisation.R
    """
}

process unsupervised_analysis {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    6_unsupervised_analysis.R
    """
}

process differential_expression {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    7_differential_expression.R
    """
}

process visualising_de_genes {
    input:
    file file_in

    output:
    file file_out

    script:
    """
    8_visualising_de_genes.R
    """
}
