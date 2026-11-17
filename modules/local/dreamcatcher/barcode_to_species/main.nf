nextflow.enable.types = true

include { perlVersion } from '../../../../lib/versions'

process BARCODE_TO_SPECIES {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        strain_umi_counts_tsv_gz: Path
    )
    gcf_species_genus: Path

    output:
    record(
        sample_id: sample_id,
        barcode_to_species_tsv: file('barcode_to_species.tsv'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'

    script:
    """
    barcode_to_species.pl ${strain_umi_counts_tsv_gz} ${gcf_species_genus} > barcode_to_species.tsv
    """

    stub:
    """
    touch barcode_to_species.tsv
    """
}
