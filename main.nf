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

workflow {

    2_load_data()
    3_sample_overview()
    4_1_segment_qc()
    4_2_probe_qc()
    4_3_aggregate_counts()
    4_4_limit_of_quantification()
    4_5_filtering()
    5_normalisation()
    6_unsupervised_analysis()
    7_differential_expression()
    8_visualising_de_genes()
}

process 2_load_data {

    script:
    """
    2_load_data.R
    """
}

process 3_sample_overview {

    script:
    """
    3_sample_overview.R
    """
}

process 4_1_segment_qc {

    script:
    """
    4_1_segment_qc.R
    """
}

process 4_2_probe_qc {

    script:
    """
    4_2_probe_qc.R
    """
}

process 4_3_aggregate_counts {

    script:
    """
    4_3_aggregate_counts.R
    """
}

process 4_4_limit_of_quantification {

    script:
    """
    4_4_limit_of_quantification.R
    """
}

process 4_5_filtering {

    script:
    """
    4_5_filtering.R
    """
}

process 5_normalisation {

    script:
    """
    5_normalisation.R
    """
}

process 6_unsupervised_analysis {

    script:
    """
    6_unsupervised_analysis.R
    """
}

process 7_differential_expression {

    script:
    """
    7_differential_expression.R
    """
}

process 8_visualising_de_genes {

    script:
    """
    8_visualising_de_genes.R
    """
}
