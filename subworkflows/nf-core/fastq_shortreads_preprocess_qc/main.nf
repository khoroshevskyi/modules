// TODO
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/umitools/extract/main'

include { PRINSEQPLUSPLUS  } from '../../../modules/nf-core/prinseqplusplus/main'
include { BBMAP_CLUMPIFY   } from '../../../modules/nf-core/bbmap/clumpify/main'

include { CAT_FASTQ        } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_SHORTREADS_PREPROCESS_QC {

    take:
    ch_reads         // channel: [ val(meta), [ fastq ] ]
    // TODO
    skip_umi_extract // boolean
    umi_discard_read // integer: 0, 1 or 2

    main:

    ch_versions = Channel.empty()
    // TODO

    // Pre - statistics
    // TODO
    // ...

    // Pre-processing
    // TODO
    // ...

    // Barcoding
    umi_reads = ch_reads
    umi_log = Channel.empty()
    if (!skip_umi_extract) {
        UMITOOLS_EXTRACT( ch_reads )
        umi_reads = UMITOOLS_EXTRACT.out.reads
        umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            UMITOOLS_EXTRACT.out.reads
                .map { meta, reads ->
                    meta.single_end ? [meta, reads] : [meta + ['single_end': true], reads[umi_discard_read % 2]]
                }
                .set { umi_reads }
        }
    }

    // Adapter removal and merging
    // TODO
    // ...

    // Complexity filtering
    // TODO
    // if (!skip_complexity_filtering) {
    //     PRINSEQPLUSPLUS( ... )
    //     ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
    // }

    // Deduplication
    // TODO
    // if (!skip_deduplication) {
    //     BBMAP_CLUMPIFY( ... )
    //     ch_versions = ch_versions.mix(BBMAP_CLUMPIFY.out.versions.first())
    // }

    // Host decontamination
    // TODO
    // ...

    // Final concatenation
    // TODO
    // if (!skip_final_concatenation) {
    //     CAT_FASTQ( ... )
    //     ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    // }

    // Post - statistics
    // TODO
    // ...


    emit:
    // TODO
    umi_log  = umi_log     // channel: [ val(meta), [ log ] ]

    versions = ch_versions // channel: [ versions.yml ]
}
