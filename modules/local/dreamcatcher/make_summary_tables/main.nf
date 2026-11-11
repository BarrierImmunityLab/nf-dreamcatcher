nextflow.enable.types = true

include { perlVersion ; toolVersion } from '../../../../lib/versions'

process MAKE_SUMMARY_TABLES {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        fcounts_tsv: Path,
        read_to_gene_tsv: Path,
        gene_biotypes_tsv: Path
    )
    seq2taxid: Path
    prefix: String

    output:
    record(
        sample_id: sample_id,
        annotated_fcounts_tsv: file("${prefix}.annotated_fcounts.tsv"),
        summary_tsv: file("${prefix}.summary.tsv"),
        protein_genes_list: file("${prefix}.protein_genes.list"),
        rrna_genes_list: file("${prefix}.rRNA_genes.list"),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'awk', version: eval(toolVersion('awk --version 2>&1 | head -n1'))) >> 'versions'

    script:
    """
    make_summary_tables.pl ${fcounts_tsv} ${gene_biotypes_tsv} ${read_to_gene_tsv} ${seq2taxid} ${prefix}
    awk '{if (\$9 != "rRNA") {print \$1}}' ${prefix}.annotated_fcounts.tsv > ${prefix}.protein_genes.list
    awk '{if (\$9 == "rRNA") {print \$1}}' ${prefix}.annotated_fcounts.tsv > ${prefix}.rRNA_genes.list
    """

    stub:
    """
    touch ${prefix}.annotated_fcounts.tsv ${prefix}.summary.tsv ${prefix}.protein_genes.list ${prefix}.rRNA_genes.list
    """
}
