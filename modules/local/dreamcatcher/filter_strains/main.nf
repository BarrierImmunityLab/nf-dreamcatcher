nextflow.enable.types = true

include { perlVersion ; toolVersion } from '../../../../lib/versions'

process FILTER_STRAINS {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        detected_summary_tsv: Path,
        detected_annotated_fcounts_tsv: Path
    )

    output:
    record(
        sample_id: sample_id,
        filtered_summary_tsv: file('filtered.summary.tsv'),
        filtered_gene_list: file('filtered.gene.list'),
        filtered_annotated_fcounts_tsv: file('filtered.annotated_fcounts.tsv'),
        filtered_protein_genes_list: file('filtered.protein_genes.list'),
        filtered_rrna_genes_list: file('filtered.rRNA_genes.list'),
        filtered_strains_list: file('filtered_strains.list'),
        blacklisted_genes_list: file('blacklisted_genes.list'),
        filter_strains_log: file('filter_strains.log'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'awk', version: eval(toolVersion('awk --version 2>&1 | head -n1'))) >> 'versions'

    script:
    """
    cp ${moduleDir}/resources/blacklisted_genes.list .
    filter_strains.pl ${detected_summary_tsv} ${detected_annotated_fcounts_tsv} blacklisted_genes.list > filter_strains.log
    grep -F -f filtered.gene.list ${detected_annotated_fcounts_tsv} > filtered.annotated_fcounts.tsv || true
    awk '{if (\$9 != "rRNA") {print \$1}}' filtered.annotated_fcounts.tsv > filtered.protein_genes.list
    awk '{if (\$9 == "rRNA") {print \$1}}' filtered.annotated_fcounts.tsv > filtered.rRNA_genes.list
    cut -f1 filtered.summary.tsv > filtered_strains.list
    """

    stub:
    """
    if [[ "${sample_id}" == *_NO_FILTERED ]]
    then
        touch filtered.summary.tsv
    else
        echo GCF_000000001.1 > filtered.summary.tsv
    fi
    touch filtered.gene.list filtered.annotated_fcounts.tsv filtered.protein_genes.list filtered.rRNA_genes.list blacklisted_genes.list filter_strains.log
    cut -f1 filtered.summary.tsv > filtered_strains.list
    """
}
