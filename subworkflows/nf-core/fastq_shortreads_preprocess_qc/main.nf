// TODO - swap modules with subworkflows when ready
// statistics
include { FASTQC as FASTQC_PRE              } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST             } from '../../../modules/nf-core/fastqc/main'
include { SEQFU_STATS as SEQFU_STATS_PRE    } from '../../../modules/nf-core/seqfu/stats/main'
include { SEQFU_STATS as SEQFU_STATS_POST   } from '../../../modules/nf-core/seqfu/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_PRE  } from '../../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_POST } from '../../../modules/nf-core/seqkit/stats/main'
include { SEQTK_COMP as SEQTK_COMP_PRE      } from '../../../modules/nf-core/seqtk/comp/main'
include { SEQTK_COMP as SEQTK_COMP_POST     } from '../../../modules/nf-core/seqtk/comp/main'
// preprocessing
// TODO include { SEQFU_CHECK                       } from '../../../modules/nf-core/seqfu/check/main'
include { SEQKIT_SANA                       } from '../../../modules/nf-core/seqkit/sana/main'
include { SEQKIT_PAIR                       } from '../../../modules/nf-core/seqkit/pair/main'
include { SEQKIT_SEQ                        } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_REPLACE                    } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQKIT_RMDUP                      } from '../../../modules/nf-core/seqkit/rmdup/main'
// barcoding
include { UMITOOLS_EXTRACT                  } from '../../../modules/nf-core/umitools/extract/main'
// adapter removal and merging
include { FASTP                             } from '../../../modules/nf-core/fastp/main'
include { TRIMGALORE                        } from '../../../modules/nf-core/trimgalore/main'
include { TRIMMOMATIC                       } from '../../../modules/nf-core/trimmomatic/main'
include { CUTADAPT                          } from '../../../modules/nf-core/cutadapt/main'
include { BBMAP_BBDUK                       } from '../../../modules/nf-core/bbmap/bbduk/main'
include { ADAPTERREMOVAL                    } from '../../../modules/nf-core/adapterremoval/main'
include { NGMERGE                           } from '../../../modules/nf-core/ngmerge/main'
include { LEEHOM                            } from '../../../modules/nf-core/leehom/main'
// complexity filtering
include { PRINSEQPLUSPLUS                   } from '../../../modules/nf-core/prinseqplusplus/main'
// deduplication
include { BBMAP_CLUMPIFY                    } from '../../../modules/nf-core/bbmap/clumpify/main'
// host decontamination
include { DEACON_INDEX                      } from '../../../modules/nf-core/deacon/index/main'
// TODO include { DEACON_FILTER                     } from '../../../modules/nf-core/deacon/filter/main'
include { HOSTILE_FETCH                     } from '../../../modules/nf-core/hostile/fetch/main'
include { HOSTILE_CLEAN                     } from '../../../modules/nf-core/hostile/clean/main'
// final concatenation
include { CAT_FASTQ                         } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_SHORTREADS_PREPROCESS_QC {

    take:
    ch_reads          // channel: [ val(meta), [ fastq ] ]
    // statistics
    skip_fastqc       // boolean
    skip_seqfu_stats  // boolean
    skip_seqkit_stats // boolean
    skip_seqtk_comp   // boolean
    // preprocessing
    // TODO skip_seqfu_check // boolean
    skip_seqkit_sana_pair // boolean
    skip_seqkit_seq       // boolean
    skip_seqkit_replace   // boolean
    skip_seqkit_rmdup     // boolean
    // barcoding
    skip_umitools_extract // boolean
    umi_discard_read      // integer: 0, 1 or 2
    // adapter removal and merging
    skip_fastp          // boolean
    skip_trimgalore     // boolean
    skip_trimmomatic    // boolean
    skip_cutadapt       // boolean
    skip_bbmap_bbduk    // boolean
    skip_adapterremoval // boolean
    skip_ngmerge        // boolean
    skip_leehom         // boolean
    // complexity filtering
    skip_prinseqplusplus // boolean
    // deduplication
    skip_bbmap_clumpify // boolean
    // host decontamination
    // TODO skip_deacon // boolean
    skip_hostile // boolean
    // final concatenation
    skip_cat_fastq // boolean

    main:

    ch_versions = Channel.empty()
    // TODO

    // pre - statistics
    // TODO
    // ...

    // preprocessing
    // TODO
    // ...

    // barcoding
    umi_reads = ch_reads
    umi_log = Channel.empty()
    if (!skip_umitools_extract) {
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

    // adapter removal and merging
    // TODO
    // ...

    // complexity filtering
    // TODO
    // if (!skip_complexity_filtering) {
    //     PRINSEQPLUSPLUS( ... )
    //     ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())
    // }

    // deduplication
    // TODO
    // if (!skip_deduplication) {
    //     BBMAP_CLUMPIFY( ... )
    //     ch_versions = ch_versions.mix(BBMAP_CLUMPIFY.out.versions.first())
    // }

    // host decontamination
    // TODO
    // ...

    // final concatenation
    // TODO
    // if (!skip_final_concatenation) {
    //     CAT_FASTQ( ... )
    //     ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    // }

    // post - statistics
    // TODO
    // ...


    emit:
    // TODO
    umi_log  = umi_log     // channel: [ val(meta), [ log ] ]

    versions = ch_versions // channel: [ versions.yml ]
}
