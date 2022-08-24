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

// Paramters for processes

// Parameters used by more than 1 processes
pathBase = "$baseDir"
slideNameLong1 = "disease"
slideNameShort1 = "d"
slideNameLong2 = "normal"
slideNameShort2 = "n"
class1 = "DKD"
class2 = "normal"
region1 = "tubule"
region2 = "glomerulus"

// load_data
pathInputData = "input_data"
pkcFiles = "Hsa_WTA_v1.0.pkc"
annotationFile = "kidney_AOI_Annotations_all_vignette.xlsx"

// segment_qc
minSegmentReads = 1000
percentTrimmed = 80
percentStitched = 80
percentAligned = 75
percentSaturation = 50
minNegativeCount = 1
maxNTCCount = 9000
minNuclei = 20
minArea = 1000

// probe_qc
minProbeRatio = 0.1
percentFailGrubbs = 20

// limit_of_quantification
minLOQ = 2
cutoffLOQ = 2

// filtering
geneDetectionRateThreshold = 0.1
genesOfInterest = "PDCD1,CD274,IFNG,CD8A,CD68,EPCAM,KRT18,NPHS1,NPHS2,CALB1,CLDN8"
detectionRateThreshold = 0.1

// visualising_de_genes
gene1 = "PDHA1"
gene2 = "ITGB1"


workflow {

    imageLoadData = load_data(
        pathBase,
        pathInputData,
        pkcFiles,
        annotationFile
    )
    imageSampleOverview = sample_overview(
        imageLoadData,
        pathBase,
        slideNameLong1,
        slideNameShort1,
        slideNameLong2,
        slideNameShort2
    )
    imageSegmentQC = segment_qc(
        imageSampleOverview,
        pathBase,
        minSegmentReads,
        percentTrimmed,
        percentStitched,
        percentAligned,
        percentSaturation,
        minNegativeCount,
        maxNTCCount,
        minNuclei,
        minArea,
        region1
    )
    imageProbeQC = probe_qc(
        imageSegmentQC,
        pathBase,
        minProbeRatio,
        percentFailGrubbs
    )
    imageAggregateCounts = aggregate_counts(
        imageProbeQC,
        pathBase
    )
    imageLimitOfQuantification = limit_of_quantification(
        imageAggregateCounts,
        pathBase,
        minLOQ,
        cutoffLOQ
    )
    imageFiltering = filtering(
        imageLimitOfQuantification,
        pathBase,
        geneDetectionRateThreshold,
        genesOfInterest,
        detectionRateThreshold,
        slideNameLong1,
        slideNameShort1,
        slideNameLong2,
        slideNameShort2
    )
    imageNormalisation = normalisation(
        imageFiltering,
        pathBase
    )
    imageUnsupervisedAnalysis = unsupervised_analysis(
        imageNormalisation,
        pathBase
    )
    imageDifferentialExpression = differential_expression(
        imageUnsupervisedAnalysis,
        pathBase,
        genesOfInterest,
        region1,
        region2,
        class1,
        class2
    )
    imageVisualingDEGenes = visualising_de_genes(
        imageDifferentialExpression,
        pathBase,
        class1,
        class2,
        region1,
        region2,
        gene1,
        gene2
    )
}

process load_data {
    input:
    path pathBase
    val pathInputData
    val pkcFiles
    val sannotationFile

    output:
    stdout emit: imageLoadData

    script:
    """
    mkdir -p data
    mkdir -p image
    mkdir -p plots
    0_load_data.R $pathBase $pathInputData $pkcFiles $annotationFile
    """
}

process sample_overview {
    input:
    file imageLoadData
    path pathBase
    val slideNameLong1
    val slideNameShort1
    val slideNameLong2
    val slideNameShort2

    output:
    stdout emit: imageSampleOverview

    script:
    """
    1_sample_overview.R $pathBase $slideNameLong1 $slideNameShort1 $slideNameLong2 $slideNameShort2
    """
}

process segment_qc {
    input:
    file imageSampleOverview
    path pathBase
    val minSegmentReads
    val percentTrimmed
    val percentStitched
    val percentAligned
    val percentSaturation
    val minNegativeCount
    val maxNTCCount
    val minNuclei
    val minArea
    val region1

    output:
    stdout emit: imageSegmentQC

    script:
    """
    2_1_segment_qc.R $pathBase $minSegmentReads $percentTrimmed $percentStitched $percentAligned \
    $percentSaturation $minNegativeCount $maxNTCCount $minNuclei $minArea $region1
    """
}

process probe_qc {
    input:
    file imageSegmentQC
    path pathBase
    val minProbeRatio
    val percentFailGrubbs

    output:
    stdout emit: imageProbeQC

    script:
    """
    2_2_probe_qc.R $pathBase $minProbeRatio $percentFailGrubbs
    """
}

process aggregate_counts {
    input:
    file imageProbeQC
    path pathBase

    output:
    stdout emit: imageAggregateCounts

    script:
    """
    2_3_aggregate_counts.R $pathBase
    """
}

process limit_of_quantification {
    input:
    file imageAggregateCounts
    path pathBase
    val minLOQ
    val cutoffLOQ

    output:
    stdout emit: imageLimitOfQuantification

    script:
    """
    2_4_limit_of_quantification.R $pathBase $minLOQ $cutoffLOQ
    """
}

process filtering {
    input:
    file imageLimitOfQuantification
    path pathBase
    val geneDetectionRateThreshold
    stdin genesOfInterest
    val detectionRateThreshold
    val slideNameLong1
    val slideNameShort1
    val slideNameLong2
    val slideNameShort2

    output:
    stdout emit: imageFiltering

    script:
    """
    2_5_filtering.R $pathBase $geneDetectionRateThreshold $genesOfInterest $detectionRateThreshold \
        $slideNameLong1 $slideNameShort1 $slideNameLong2 $slideNameShort2
    """
}

process normalisation {
    input:
    file imageFiltering
    path pathBase

    output:
    stdout emit: imageNormalisation

    script:
    """
    3_normalisation.R $pathBase
    """
}

process unsupervised_analysis {
    input:
    file imageNormalisation
    path pathBase

    output:
    stdout emit: imageUnsupervisedAnalysis

    script:
    """
    4_unsupervised_analysis.R $pathBase
    """
}

process differential_expression {
    input:
    file imageUnsupervisedAnalysis
    path pathBase
    stdin genesOfInterest
    val region1
    val region2
    val class1
    val class2
    
    output:
    stdout emit: imageDifferentialExpression

    script:
    """
    5_differential_expression.R $pathBase $genesOfInterest $region1 $region2 $class1 $class2
    """
}

process visualising_de_genes {
    input:
    file imageDifferentialExpression
    path pathBase
    val class1
    val class2
    val region1
    val region2
    val gene1
    val gene2

    output:
    stdout emit: imageUnsupervisedAnalysis

    script:
    """
    6_visualising_de_genes.R $pathBase $class1 $class2 $region1 $region2 $gene1 $gene2
    """
}
