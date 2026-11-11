nextflow.enable.types = true

include { perlVersion ; pigzVersion } from '../../../../lib/versions'

process CORRECT_BARCODES {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        uncorrected_fastq: Path,
        barcodes_vs_whitelist_tsv: Path,
        make_umi_fastq_log: Path,
        get_barcode_corrections_log: Path
    )

    output:
    record(
        sample_id: sample_id,
        reads_r1: file('Unmapped_filt.R1.fastq.gz'),
        preprocess_log: file('preprocess_reads.log'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'pigz', version: eval(pigzVersion())) >> 'versions'

    script:
    """
    correct_barcodes.pl ${uncorrected_fastq} ${barcodes_vs_whitelist_tsv} Unmapped_filt.R1.fastq > correct_barcodes.log
    pigz -p${task.cpus} -f Unmapped_filt.R1.fastq
    cat ${make_umi_fastq_log} ${get_barcode_corrections_log} correct_barcodes.log > preprocess_reads.log
    """

    stub:
    """
    echo | gzip > Unmapped_filt.R1.fastq.gz
    touch preprocess_reads.log
    """
}
