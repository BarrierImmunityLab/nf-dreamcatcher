nextflow.enable.types = true

include { DreamcatcherInput ; PreprocessedReads } from '../../lib/types.nf'
include { MAKE_UMI_FASTQ } from '../../modules/local/preprocess_reads/make_umi_fastq'
include { GET_BARCODE_CORRECTIONS } from '../../modules/local/preprocess_reads/get_barcode_corrections'
include { CORRECT_BARCODES } from '../../modules/local/preprocess_reads/correct_barcodes'
include { MAKE_BULK_FASTQ } from '../../modules/local/preprocess_reads/make_bulk_fastq'

workflow PREPROCESS_READS {
    take:
    samples: Channel<DreamcatcherInput>

    main:
    // 10x: UMI-tools formatted single-end reads with corrected barcodes
    ch_10x = samples.filter { r -> r.instrument == '10x' }
    ch_umi_fastq = MAKE_UMI_FASTQ(ch_10x)
    ch_barcode_corrections = GET_BARCODE_CORRECTIONS(ch_10x.join(ch_umi_fastq, by: 'sample_id'))
    ch_corrected = CORRECT_BARCODES(ch_umi_fastq.join(ch_barcode_corrections, by: 'sample_id')).map { r -> record(sample_id: r.sample_id, reads_r1: r.reads_r1, reads_r2: null, preprocess_log: r.preprocess_log) }

    // Bulk: single-end or paired-end reads
    ch_bulk = MAKE_BULK_FASTQ(samples.filter { r -> r.instrument != '10x' })

    emit:
    reads: Channel<PreprocessedReads> = ch_corrected.mix(ch_bulk)
}
