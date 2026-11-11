nextflow.enable.types = true

include { rPackageVersion ; rVersion } from '../../../../lib/versions'

process PARSE_READ_NETWORK {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-r:1.0.0"

    input:
    record(
        sample_id: String,
        filtered_summary_tsv: Path,
        rrna_uniq_counts_tsv: Path,
        protein_uniq_counts_tsv: Path,
        nodes_and_edges_tsv: Path,
        rrna_nodes_and_edges_tsv: Path,
        protein_nodes_and_edges_tsv: Path
    )
    gcf_species_genus: Path

    output:
    record(
        sample_id: sample_id,
        filtered_cluster_tsv: file('filtered.cluster.tsv'),
        top_cluster_tsv: file('top.cluster.tsv'),
        top_strains_list: file('top_strains.list'),
        parse_read_network_log: file('parse_read_network.log'),
    )

    topic:
    record(process: task.process, tool: 'r', version: eval(rVersion())) >> 'versions'
    record(process: task.process, tool: 'igraph', version: eval(rPackageVersion('igraph'))) >> 'versions'
    record(process: task.process, tool: 'dplyr', version: eval(rPackageVersion('dplyr'))) >> 'versions'

    script:
    """
    parse_read_network.R ${gcf_species_genus} &> parse_read_network.log
    sort -u -t\$'\\t' -k15,15 --merge filtered.cluster.tsv > top.cluster.tsv
    cut -f1 top.cluster.tsv | grep -v RefSeq > top_strains.list
    """

    stub:
    """
    touch filtered.cluster.tsv top.cluster.tsv parse_read_network.log
    echo GCF_000000001.1 > top_strains.list
    """
}
