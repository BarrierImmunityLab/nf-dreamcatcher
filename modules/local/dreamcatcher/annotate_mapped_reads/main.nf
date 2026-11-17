nextflow.enable.types = true

include { perlVersion ; samtoolsVersion } from '../../../../lib/versions'

process ANNOTATE_MAPPED_READS {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        fcounts_bam: Path,
        human_bam: Path
    )
    prefix: String

    output:
    record(
        sample_id: sample_id,
        read_to_gene_tsv: file("${prefix}.read_to_gene_w_mismatches.tsv"),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'samtools', version: eval(samtoolsVersion())) >> 'versions'

    script:
    """
    export CPUS=${task.cpus}
    annotate_mapped_reads.pl ${fcounts_bam} ${human_bam} > ${prefix}.read_to_gene_w_mismatches.tsv
    """

    stub:
    """
    touch ${prefix}.read_to_gene_w_mismatches.tsv
    """
}
