nextflow.enable.types = true

include { perlVersion } from '../../../../lib/versions'

process EXTRACT_GENE_BIOTYPES {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        gff: Path
    )

    output:
    record(
        sample_id: sample_id,
        gene_biotypes_tsv: file('gene_biotypes.tsv'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'

    script:
    """
    extract_gene_biotypes.pl ${gff} > gene_biotypes.tsv
    """

    stub:
    """
    touch gene_biotypes.tsv
    """
}
