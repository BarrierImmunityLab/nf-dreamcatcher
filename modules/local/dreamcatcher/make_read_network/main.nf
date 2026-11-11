nextflow.enable.types = true

include { perlVersion } from '../../../../lib/versions'

process MAKE_READ_NETWORK {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        reads_tsv: Path,
        filtered_annotated_fcounts_tsv: Path
    )
    output_name: String

    output:
    record(
        sample_id: sample_id,
        nodes_and_edges_tsv: file(output_name),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'

    script:
    """
    make_read_network.pl ${reads_tsv} ${filtered_annotated_fcounts_tsv} > ${output_name}
    """

    stub:
    """
    touch ${output_name}
    """
}
