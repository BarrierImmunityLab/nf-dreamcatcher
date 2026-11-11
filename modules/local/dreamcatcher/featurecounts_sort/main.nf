nextflow.enable.types = true

include { samtoolsVersion ; toolVersion } from '../../../../lib/versions'

process FEATURECOUNTS_SORT {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-mapping:1.0.0"

    input:
    record(
        sample_id: String,
        bam: Path,
        gff: Path,
        paired: Boolean
    )
    level: String

    output:
    record(
        sample_id: sample_id,
        counts_tsv: file(level == 'gene' ? 'top.fcounts.tsv' : 'top.strain_fcounts.tsv'),
        counts_summary: file(level == 'gene' ? 'top.fcounts.tsv.summary' : 'top.strain_fcounts.tsv.summary'),
        counts_log: file("top.${level}_fcounts.log"),
        sorted_bam: file("top_filtered_${level}_sorted.bam"),
        sorted_bam_bai: file("top_filtered_${level}_sorted.bam.bai"),
    )

    topic:
    record(process: task.process, tool: 'subread', version: eval(toolVersion("featureCounts -v 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'
    record(process: task.process, tool: 'samtools', version: eval(samtoolsVersion())) >> 'versions'

    script:
    if (!(level in ['gene', 'strain'])) {
        error("Unknown counting level: ${level}")
    }
    def attribute = level == 'gene' ? 'locus_tag' : 'strain_id'
    def counts = level == 'gene' ? 'top.fcounts.tsv' : 'top.strain_fcounts.tsv'
    def paired_args = paired ? '-p --countReadPairs' : ''
    """
    featureCounts -T ${task.cpus} ${paired_args} -M -O --fraction -t gene -g ${attribute} -s 0 -a ${gff} \\
        -o ${counts} -R BAM ${bam} &> top.${level}_fcounts.log
    samtools sort -@${task.cpus} ${bam}.featureCounts.bam > top_filtered_${level}_sorted.bam 2> /dev/null
    samtools index top_filtered_${level}_sorted.bam
    """

    stub:
    """
    touch top.fcounts.tsv top.fcounts.tsv.summary top.strain_fcounts.tsv top.strain_fcounts.tsv.summary top.${level}_fcounts.log top_filtered_${level}_sorted.bam top_filtered_${level}_sorted.bam.bai
    """
}
