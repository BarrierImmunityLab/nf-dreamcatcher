nextflow.enable.types = true

include { toolVersion } from '../../../../lib/versions'

process FEATURECOUNTS {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-featurecounts:1.0.0"

    input:
    record(
        sample_id: String,
        bam: Path,
        gff: Path,
        paired: Boolean
    )
    prefix: String

    output:
    record(
        sample_id: sample_id,
        fcounts_tsv: file("${prefix}.fcounts.tsv"),
        fcounts_summary: file("${prefix}.fcounts.tsv.summary"),
        fcounts_bam: file("${prefix}_fcounts.bam"),
        fcounts_log: file("${prefix}.fcounts.log"),
    )

    topic:
    record(process: task.process, tool: 'subread', version: eval(toolVersion("featureCounts -v 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'

    script:
    def paired_args = paired ? '-p --countReadPairs' : ''
    """
    featureCounts -T ${task.cpus} ${paired_args} -M -O --fraction -t gene -g locus_tag -s 0 -a ${gff} \\
        -o ${prefix}.fcounts.tsv -R BAM ${bam} &> ${prefix}.fcounts.log
    mv ${bam}.featureCounts.bam ${prefix}_fcounts.bam
    """

    stub:
    """
    touch ${prefix}.fcounts.tsv ${prefix}.fcounts.tsv.summary ${prefix}_fcounts.bam ${prefix}.fcounts.log
    """
}
