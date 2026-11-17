nextflow.enable.types = true

include { perlVersion ; toolVersion } from '../../../../lib/versions'

process FILTER_TOP_GFF {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        top_annotated_fcounts_tsv: Path,
        top_gff: Path,
        blacklisted_genes_list: Path
    )
    seq2taxid: Path

    output:
    record(
        sample_id: sample_id,
        top_removed_genes_list: file('top_removed.genes.list'),
        top_bad_genes_list: file('top_bad.genes.list'),
        top_filtered_gff: file('top_strains_filtered.gff'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'
    record(process: task.process, tool: 'awk', version: eval(toolVersion('awk --version 2>&1 | head -n1'))) >> 'versions'

    script:
    """
    awk -F '\\t' '(\$9 == "rRNA" && \$10 > 3.4) || (\$9 != "rRNA" && \$10 > 2.3) || \$11 > 0.1' ${top_annotated_fcounts_tsv} | cut -f1 > top_removed.genes.list
    cat top_removed.genes.list ${blacklisted_genes_list} | sort | uniq > top_bad.genes.list
    filter_top_gff.pl ${top_gff} top_bad.genes.list ${seq2taxid} > top_strains_filtered.gff
    """

    stub:
    """
    touch top_removed.genes.list top_bad.genes.list top_strains_filtered.gff
    """
}
