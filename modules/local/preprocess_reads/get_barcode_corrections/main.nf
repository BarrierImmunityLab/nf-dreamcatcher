nextflow.enable.types = true

include { perlVersion ; toolVersion } from '../../../../lib/versions'

process GET_BARCODE_CORRECTIONS {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        unmapped_reads_with_logs: Iterable<Path>,
        uncorrected_fastq: Path
    )

    stage:
    stageAs unmapped_reads_with_logs, 'reads/*'

    output:
    record(
        sample_id: sample_id,
        barcodes_vs_whitelist_tsv: file('barcodes_vs_whitelist.tsv'),
        barcodes_to_correct_fa: file('barcodes_to_correct.fa'),
        whitelist_fa: file('whitelist.fa'),
        barcode_bowtie2_log: file('barcode_bowtie2.log'),
        get_barcode_corrections_log: file('get_barcode_corrections.log'),
    )

    topic:
    record(process: task.process, tool: 'bowtie2', version: eval(toolVersion("bowtie2 --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'

    script:
    """
    get_barcode_corrections.sh reads ${task.cpus} > get_barcode_corrections.log
    """

    stub:
    """
    touch barcodes_vs_whitelist.tsv barcodes_to_correct.fa whitelist.fa barcode_bowtie2.log get_barcode_corrections.log
    """
}
