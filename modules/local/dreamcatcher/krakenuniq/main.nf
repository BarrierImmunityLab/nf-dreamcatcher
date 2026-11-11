nextflow.enable.types = true

include { pigzVersion ; toolVersion } from '../../../../lib/versions'

process KRAKENUNIQ {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-krakenuniq:1.0.0"

    input:
    record(
        sample_id: String,
        reads_r1: Path,
        reads_r2: Path?
    )
    kuniq_db: Path

    output:
    record(
        sample_id: sample_id,
        kuniq_output_gz: file('kuniq.output.txt.gz'),
        kuniq_report_txt: file('kuniq.report.txt'),
        kuniq_log: file('kuniq.log'),
    )

    topic:
    record(process: task.process, tool: 'krakenuniq', version: eval(toolVersion("krakenuniq --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'
    record(process: task.process, tool: 'pigz', version: eval(pigzVersion())) >> 'versions'

    script:
    def preload_gb = Math.floor(task.memory.toGiga() * 0.85).intValue()
    def reads = reads_r2 ? "--paired ${reads_r1} ${reads_r2}" : "${reads_r1}"
    """
    krakenuniq --threads ${task.cpus} --db ${kuniq_db} --preload-size ${preload_gb}G --output kuniq.output.txt \\
        --report-file kuniq.report.txt ${reads} &> kuniq.log
    pigz -p${task.cpus} -f kuniq.output.txt

    if (( \$(wc -l < kuniq.report.txt) < 10 ))
    then
        echo "ERROR: KrakenUniq has produced an empty report file!" >&2
        exit 1
    fi
    """

    stub:
    """
    echo | gzip > kuniq.output.txt.gz
    touch kuniq.report.txt kuniq.log
    """
}
