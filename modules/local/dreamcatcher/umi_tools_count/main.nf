nextflow.enable.types = true

include { toolVersion } from '../../../../lib/versions'

process UMI_TOOLS_COUNT {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        sorted_bam: Path,
        sorted_bam_bai: Path
    )
    level: String

    output:
    record(
        sample_id: sample_id,
        umi_counts_tsv_gz: file("${level}.umi_counts.tsv.gz"),
        umi_tools_count_log: file("${level}.umi_tools_count.log"),
    )

    topic:
    record(process: task.process, tool: 'umi_tools', version: eval(toolVersion("umi_tools --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'

    script:
    """
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ${sorted_bam} -S ${level}.umi_counts.tsv.gz &> ${level}.umi_tools_count.log
    """

    stub:
    """
    echo | gzip > ${level}.umi_counts.tsv.gz
    touch ${level}.umi_tools_count.log
    """
}
