nextflow.enable.types = true

include { samtoolsVersion ; toolVersion } from '../../../../lib/versions'

process HISAT2_MAP {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-mapping:1.0.0"

    input:
    record(
        sample_id: String,
        reads_r1: Path,
        reads_r2: Path?,
        index_files: Iterable<Path>
    )
    index_name: String
    align_args: String
    bam_name: String
    log_name: String

    stage:
    stageAs index_files, '*'

    output:
    record(
        sample_id: sample_id,
        bam: file(bam_name),
        hisat2_log: file(log_name),
    )

    topic:
    record(process: task.process, tool: 'hisat2', version: eval(toolVersion("hisat2 --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'
    record(process: task.process, tool: 'samtools', version: eval(samtoolsVersion())) >> 'versions'

    script:
    def reads = reads_r2 ? "-1 ${reads_r1} -2 ${reads_r2}" : "-U ${reads_r1}"
    """
    set -o pipefail
    hisat2 -t -p${task.cpus} ${align_args} ${reads} -x ${index_name} 2> ${log_name} \\
        | samtools sort -@${task.cpus} --verbosity 0 -O BAM - > ${bam_name}
    """

    stub:
    """
    touch ${bam_name} ${log_name}
    """
}
