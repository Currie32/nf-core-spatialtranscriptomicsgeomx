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
    image_visualing_de_genes = visualising_de_genes(params.path_to_data, image_differential_expression)
}

process load_data {
    input:
    path path_to_data

    output:
    stdout emit: image_load_data

    script:
    """
    mkdir -p data
    mkdir -p image
    mkdir -p plots
    0_load_data.R $path_to_data
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
    1_sample_overview.R $path_to_data
    """
}

process segment_qc {
    input:
    path path_to_data
    file image_sample_overview

    output:
    stdout emit: image_segment_qc

    script:
    """
    2_1_segment_qc.R $path_to_data
    """
}

process probe_qc {
    input:
    path path_to_data
    file image_segment_qc

    output:
    stdout emit: image_probe_qc

    script:
    """
    2_2_probe_qc.R $path_to_data
    """
}

process aggregate_counts {
    input:
    path path_to_data
    file image_probe_qc

    output:
    stdout emit: image_aggregate_counts

    script:
    """
    2_3_aggregate_counts.R $path_to_data
    """
}

process limit_of_quantification {
    input:
    path path_to_data
    file image_aggregate_counts

    output:
    stdout emit: image_limit_of_quantification

    script:
    """
    2_4_limit_of_quantification.R $path_to_data
    """
}

process filtering {
    input:
    path path_to_data
    file image_limit_of_quantification

    output:
    stdout emit: image_filtering

    script:
    """
    2_5_filtering.R $path_to_data
    """
}

process normalisation {
    input:
    path path_to_data
    file image_filtering

    output:
    stdout emit: image_normalisation

    script:
    """
    3_normalisation.R $path_to_data
    """
}

process unsupervised_analysis {
    input:
    path path_to_data
    file image_normalisation

    output:
    stdout emit: image_unsupervised_analysis

    script:
    """
    4_unsupervised_analysis.R $path_to_data
    """
}

process differential_expression {
    input:
    path path_to_data
    file image_unsupervised_analysis

    output:
    stdout emit: image_differential_expression

    script:
    """
    5_differential_expression.R $path_to_data
    """
}

process visualising_de_genes {
    input:
    path path_to_data
    file image_differential_expression

    output:
    stdout emit: image_unsupervised_analysis

    script:
    """
    6_visualising_de_genes.R $path_to_data
    """
}
