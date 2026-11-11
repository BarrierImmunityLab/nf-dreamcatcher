nextflow.enable.types = true

include { toolVersion } from '../../../../lib/versions'

process MAKE_COMBINED_REFERENCE {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        strains_list: Path
    )
    ncbi_dir: Path
    prefix: String

    output:
    record(
        sample_id: sample_id,
        fna: file("${prefix}_strains.fna"),
        gff: file("${prefix}_strains.gff"),
        mkref_log: file("${prefix}.mkref.log"),
    )

    topic:
    record(process: task.process, tool: 'gzip', version: eval(toolVersion('gzip --version 2>&1 | head -n1'))) >> 'versions'

    script:
    """
    make_combined_reference.sh ${ncbi_dir} ${strains_list} ${prefix}_strains &> ${prefix}.mkref.log
    """

    stub:
    """
    touch ${prefix}_strains.fna ${prefix}_strains.gff ${prefix}.mkref.log
    """
}
