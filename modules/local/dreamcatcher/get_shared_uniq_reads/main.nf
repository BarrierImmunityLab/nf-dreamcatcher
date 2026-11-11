nextflow.enable.types = true

include { perlVersion ; toolVersion } from '../../../../lib/versions'

process GET_SHARED_UNIQ_READS {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        detected_read_to_gene_tsv: Path,
        filtered_gene_list: Path,
        filtered_annotated_fcounts_tsv: Path
    )

    output:
    record(
        sample_id: sample_id,
        filtered_read_type_strain_tsv: file('filtered.read_type_strain.tsv'),
        reads_by_strain_uniqness_tsv: file('reads_by_strain_uniqness.tsv'),
        filtered_uniq_strain_reads_list: file('filtered.uniq_stain_reads.list'),
        filtered_multi_strain_reads_list: file('filtered.multi_strain_reads.list'),
        rrna_uniq_counts_tsv: file('rRNA.uniq_counts.tsv'),
        protein_uniq_counts_tsv: file('protein.uniq_counts.tsv'),
        rrna_shared_reads_tsv: file('rRNA.shared_reads.tsv'),
        protein_shared_reads_tsv: file('protein.shared_reads.tsv'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'sort', version: eval(toolVersion('sort --version 2>&1 | head -n1'))) >> 'versions'

    script:
    def sort_ram_gb = Math.floor(task.memory.toGiga() * 0.85).intValue()
    """
    export RAM=${sort_ram_gb}
    export CPUS=${task.cpus}
    get_shared_uniq_reads.sh
    """

    stub:
    """
    touch filtered.read_type_strain.tsv reads_by_strain_uniqness.tsv filtered.uniq_stain_reads.list filtered.multi_strain_reads.list rRNA.uniq_counts.tsv protein.uniq_counts.tsv rRNA.shared_reads.tsv protein.shared_reads.tsv
    """
}
