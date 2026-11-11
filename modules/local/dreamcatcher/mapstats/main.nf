nextflow.enable.types = true

include { perlVersion ; samtoolsVersion } from '../../../../lib/versions'

process MAPSTATS {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        preprocess_log: Path,
        kuniq_report_txt: Path,
        detected_strains_list: Path,
        detected_hisat2_log: Path?,
        detected_protein_genes_list: Path?,
        detected_rrna_genes_list: Path?,
        detected_read_to_gene_tsv: Path?,
        detected_fcounts_bam: Path?,
        detected_annotated_fcounts_tsv: Path?,
        filtered_strains_list: Path?,
        filtered_read_type_strain_tsv: Path?,
        filtered_protein_genes_list: Path?,
        filtered_rrna_genes_list: Path?,
        top_strains_list: Path?,
        top_hisat2_log: Path?,
        top_protein_genes_list: Path?,
        top_rrna_genes_list: Path?,
        top_sorted_bam: Path?,
        filtered_annotated_fcounts_tsv: Path?
    )

    stage:
    stageAs preprocess_log, '0_preprocessed_reads/*'
    stageAs kuniq_report_txt, '1_krakenuniq/*'
    stageAs detected_strains_list, '2_detected_strains/*'
    stageAs detected_hisat2_log, '2_detected_strains/*'
    stageAs detected_protein_genes_list, '2_detected_strains/*'
    stageAs detected_rrna_genes_list, '2_detected_strains/*'
    stageAs detected_read_to_gene_tsv, '2_detected_strains/*'
    stageAs detected_fcounts_bam, '2_detected_strains/*'
    stageAs detected_annotated_fcounts_tsv, '2_detected_strains/*'
    stageAs filtered_strains_list, '3_filtered_strains/*'
    stageAs filtered_read_type_strain_tsv, '3_filtered_strains/*'
    stageAs filtered_protein_genes_list, '3_filtered_strains/*'
    stageAs filtered_rrna_genes_list, '3_filtered_strains/*'
    stageAs filtered_annotated_fcounts_tsv, '3_filtered_strains/*'
    stageAs top_strains_list, '5_top_strains/*'
    stageAs top_hisat2_log, '5_top_strains/*'
    stageAs top_protein_genes_list, '5_top_strains/*'
    stageAs top_rrna_genes_list, '5_top_strains/*'
    stageAs top_sorted_bam, '5_top_strains/*'

    output:
    record(
        sample_id: sample_id,
        mapstats_tsv: file('mapstats.tsv'),
        filtered_chr_list: file('3_filtered_strains/filtered.chr.list', optional: true),
    )

    topic:
    record(process: task.process, tool: 'samtools', version: eval(samtoolsVersion())) >> 'versions'
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'

    script:
    """
    export CPUS=${task.cpus}
    calculate_mapstats.sh ${sample_id} > mapstats.tsv
    """

    stub:
    """
    printf 'Sample\\tTot_reads\\n${sample_id}\\t1\\n' > mapstats.tsv
    """
}
