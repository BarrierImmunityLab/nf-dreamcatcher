nextflow.enable.types = true

include { perlVersion ; pigzVersion ; toolVersion } from '../../../../lib/versions'

process MAKE_BULK_FASTQ {
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
        reads_r1: file('Unmapped_filt.R1.fastq.gz'),
        reads_r2: file('Unmapped_filt.R2.fastq.gz', optional: true),
        preprocess_log: file('preprocess_reads.log'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'pigz', version: eval(pigzVersion())) >> 'versions'
    record(process: task.process, tool: 'pbzip2', version: eval(toolVersion("pbzip2 --version 2>&1 | grep -Eo 'v[0-9.]+' | head -n1 | sed 's/v//'"))) >> 'versions'

    script:
    """
    export CPUS=${task.cpus}
    echo -e "\\tpreprocess_unmapped_reads.sh: preparing unmapped reads in bulk mode.." > preprocess_reads.log
    make_bulk_fastq.pl reads >> preprocess_reads.log
    pigz -p${task.cpus} -f Unmapped_filt.R1.fastq
    if [[ -s Unmapped_filt.R2.fastq ]]
    then
        pigz -p${task.cpus} -f Unmapped_filt.R2.fastq
    fi
    """

    stub:
    """
    echo | gzip > Unmapped_filt.R1.fastq.gz
    if [ "${unmapped_reads_with_logs.any { f -> f.name.startsWith('Unmapped.out.mate2') }}" = "true" ]; then
        echo | gzip > Unmapped_filt.R2.fastq.gz
    fi
    touch preprocess_reads.log
    """
}
