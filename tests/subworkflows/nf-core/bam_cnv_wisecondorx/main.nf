#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM_CNV_WISECONDORX } from '../../../../subworkflows/nf-core/bam_cnv_wisecondorx/main'

workflow test_bam_cnv_wisecondorx {

    input = [
        [id: 'test'],
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true),
    ]

    fasta = [
        [id: 'fasta'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ]

    fai = [
        [id: 'fai'],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true),
    ]

    reference = [
        [id: 'wc_reference'],
        [],
    ]

    blacklist = [
        [id: 'blacklist'],
        file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
    ]

    BAM_CNV_WISECONDORX(
        input,
        fasta,
        fai,
        reference,
        blacklist,
    )
}
