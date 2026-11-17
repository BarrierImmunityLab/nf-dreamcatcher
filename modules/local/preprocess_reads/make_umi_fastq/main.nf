nextflow.enable.types = true

include { perlVersion ; pigzVersion ; toolVersion } from '../../../../lib/versions'

process MAKE_UMI_FASTQ {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        unmapped_reads_with_logs: Iterable<Path>
    )

    stage:
    stageAs unmapped_reads_with_logs, 'reads/*'

    output:
    record(
        sample_id: sample_id,
        uncorrected_fastq: file('Unmapped_filt_uncorr.R1.fastq'),
        make_umi_fastq_log: file('make_umi_fastq.log'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'pigz', version: eval(pigzVersion())) >> 'versions'
    record(process: task.process, tool: 'pbzip2', version: eval(toolVersion("pbzip2 --version 2>&1 | grep -Eo 'v[0-9.]+' | head -n1 | sed 's/v//'"))) >> 'versions'

    script:
    """
    export CPUS=${task.cpus}
    echo -e "\\tpreprocess_unmapped_reads.sh: preparing unmapped reads in single cell mode, STARsolo raw unmapped reads are used." > make_umi_fastq.log
    make_umi_fastq.pl reads Unmapped_filt_uncorr.R1.fastq >> make_umi_fastq.log
    """

    stub:
    """
    touch Unmapped_filt_uncorr.R1.fastq make_umi_fastq.log
    """
}
