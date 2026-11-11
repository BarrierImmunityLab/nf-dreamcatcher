nextflow.enable.types = true

include { toolVersion } from '../../../../lib/versions'

process HISAT2_BUILD {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-hisat2:1.0.0"

    input:
    record(
        sample_id: String,
        fna: Path
    )
    prefix: String

    output:
    record(
        sample_id: sample_id,
        index_files: files("${prefix}_strains.*.ht2*"),
        build_log: file("${prefix}.hisat2_build.log"),
    )

    topic:
    record(process: task.process, tool: 'hisat2', version: eval(toolVersion("hisat2 --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1"))) >> 'versions'

    script:
    """
    hisat2-build ${fna} ${prefix}_strains &> ${prefix}.hisat2_build.log
    """

    stub:
    """
    touch ${prefix}_strains.1.ht2 ${prefix}.hisat2_build.log
    """
}
