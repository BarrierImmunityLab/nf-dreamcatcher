nextflow.enable.types = true

include { DreamcatcherInput ; DreamcatcherResult } from '../../lib/types.nf'
include { PREPROCESS_READS } from './preprocess_reads'
include { KRAKENUNIQ } from '../../modules/local/dreamcatcher/krakenuniq'
include { KUNIQ_TO_GCF } from '../../modules/local/dreamcatcher/kuniq_to_gcf'
include { MAKE_COMBINED_REFERENCE as MAKE_COMBINED_REFERENCE_DETECTED } from '../../modules/local/dreamcatcher/make_combined_reference'
include { MAKE_COMBINED_REFERENCE as MAKE_COMBINED_REFERENCE_TOP } from '../../modules/local/dreamcatcher/make_combined_reference'
include { HISAT2_BUILD as HISAT2_BUILD_DETECTED } from '../../modules/local/dreamcatcher/hisat2_build'
include { HISAT2_BUILD as HISAT2_BUILD_TOP } from '../../modules/local/dreamcatcher/hisat2_build'
include { HISAT2_MAP as HISAT2_MAP_HUMAN } from '../../modules/local/dreamcatcher/hisat2_map'
include { HISAT2_MAP as HISAT2_MAP_DETECTED } from '../../modules/local/dreamcatcher/hisat2_map'
include { HISAT2_MAP as HISAT2_MAP_TOP } from '../../modules/local/dreamcatcher/hisat2_map'
include { FEATURECOUNTS as FEATURECOUNTS_DETECTED } from '../../modules/local/dreamcatcher/featurecounts'
include { FEATURECOUNTS as FEATURECOUNTS_TOP } from '../../modules/local/dreamcatcher/featurecounts'
include { FEATURECOUNTS_SORT as FEATURECOUNTS_SORT_GENE } from '../../modules/local/dreamcatcher/featurecounts_sort'
include { FEATURECOUNTS_SORT as FEATURECOUNTS_SORT_STRAIN } from '../../modules/local/dreamcatcher/featurecounts_sort'
include { ANNOTATE_MAPPED_READS as ANNOTATE_MAPPED_READS_DETECTED } from '../../modules/local/dreamcatcher/annotate_mapped_reads'
include { ANNOTATE_MAPPED_READS as ANNOTATE_MAPPED_READS_TOP } from '../../modules/local/dreamcatcher/annotate_mapped_reads'
include { EXTRACT_GENE_BIOTYPES } from '../../modules/local/dreamcatcher/extract_gene_biotypes'
include { MAKE_SUMMARY_TABLES as MAKE_SUMMARY_TABLES_DETECTED } from '../../modules/local/dreamcatcher/make_summary_tables'
include { MAKE_SUMMARY_TABLES as MAKE_SUMMARY_TABLES_TOP_INITIAL } from '../../modules/local/dreamcatcher/make_summary_tables'
include { MAKE_SUMMARY_TABLES as MAKE_SUMMARY_TABLES_TOP } from '../../modules/local/dreamcatcher/make_summary_tables'
include { FILTER_STRAINS } from '../../modules/local/dreamcatcher/filter_strains'
include { GET_SHARED_UNIQ_READS } from '../../modules/local/dreamcatcher/get_shared_uniq_reads'
include { MAKE_READ_NETWORK as MAKE_READ_NETWORK_ALL } from '../../modules/local/dreamcatcher/make_read_network'
include { MAKE_READ_NETWORK as MAKE_READ_NETWORK_RRNA } from '../../modules/local/dreamcatcher/make_read_network'
include { MAKE_READ_NETWORK as MAKE_READ_NETWORK_PROTEIN } from '../../modules/local/dreamcatcher/make_read_network'
include { PARSE_READ_NETWORK } from '../../modules/local/dreamcatcher/parse_read_network'
include { FILTER_TOP_GFF } from '../../modules/local/dreamcatcher/filter_top_gff'
include { UMI_TOOLS_COUNT as UMI_TOOLS_COUNT_GENE } from '../../modules/local/dreamcatcher/umi_tools_count'
include { UMI_TOOLS_COUNT as UMI_TOOLS_COUNT_STRAIN } from '../../modules/local/dreamcatcher/umi_tools_count'
include { BARCODE_TO_SPECIES } from '../../modules/local/dreamcatcher/barcode_to_species'
include { MAPSTATS } from '../../modules/local/dreamcatcher/mapstats'

workflow DREAMCATCHER {
    take:
    samples: Channel<DreamcatcherInput>

    main:
    ncbi_dir = file("${params.dreamcatcher_refdata_root}/NCBI_datasets/bacteria_ref/ncbi_dataset/data")
    kuniq_db = file("${params.dreamcatcher_refdata_root}/KrakenUniq_select")
    seq2taxid = file("${params.dreamcatcher_refdata_root}/KrakenUniq_select/seqid2taxid.map.refseq")
    assembly_summary = file("${params.dreamcatcher_refdata_root}/KrakenUniq_select/assembly_summary_bacteria.txt")
    gcf_species_genus = file("${params.dreamcatcher_refdata_root}/KrakenUniq_select/gcf_species_genus.tsv")
    human_index_files = files("${params.dreamcatcher_refdata_root}/decoy_reference/T2T_mito_phix.*.ht2*")

    // STEP 1: Extract and preprocess the unmapped reads for further analysis..
    // 6 distinct cases: bulk/single-cell mode, paired-end/single-end reads. Only 10x is processed in single cell mode, and 10x reads are always coerced to UMI-Tools formatted single-end reads
    // Bulk reads are homopolymer-filtered and discarded if too short (<40 bp). Single end reads are homopolymer-filtered and barcodes are corrected vs the whitelist in case of STARsolo unmapped reads
    ch_reads = PREPROCESS_READS(samples).map { r -> record(sample_id: r.sample_id, reads_r1: r.reads_r1, reads_r2: r.reads_r2, paired: r.reads_r2 != null, preprocess_log: r.preprocess_log) }
    ch_paired = ch_reads.map { r -> record(sample_id: r.sample_id, paired: r.paired) }
    ch_10x = samples.filter { r -> r.instrument == '10x' }.map { r -> record(sample_id: r.sample_id) }

    // STEP 2: Running KrakenUniq ..
    // KrakenUniq quits with no error and no output in some cases, still not sure why. Just re-run the job (see the retries in conf/pbspro.config).
    ch_kuniq = KRAKENUNIQ(ch_reads, kuniq_db)

    // STEP 3: Extracting unique bacterial RefSeq accessions ..
    // This finds the best RefSeq genome for each detected species, using KrakenUniq output. If no strains passed the filtering, we quit, but we DO generate a mapstats.tsv (mostly filled with zeroes).
    ch_detected_strains = KUNIQ_TO_GCF(ch_kuniq, gcf_species_genus, assembly_summary)
    ch_with_detected_strains = ch_detected_strains.filter { r -> r.detected_strains_list.size() > 0 }

    // STEP 4: Generating a combined genome and GFF file for all bacterial species ..
    ch_detected_ref = MAKE_COMBINED_REFERENCE_DETECTED(
        ch_with_detected_strains.map { r -> record(sample_id: r.sample_id, strains_list: r.detected_strains_list) },
        ncbi_dir,
        'detected',
    )
    // STEP 5: Making combined Hisat2 index ..
    ch_detected_index = HISAT2_BUILD_DETECTED(ch_detected_ref, 'detected')
    // STEP 6: Mapping original unmapped reads to the combined reference using Hisat2 ..
    // Notice the differential use of "-k/--max-seeds" and the relaxed soft-clipping penalty.
    ch_detected_map = HISAT2_MAP_DETECTED(
        ch_reads.join(ch_detected_index, by: 'sample_id'),
        'detected_strains',
        '-k 5000 --sp 1,1 --no-spliced-alignment --no-unal',
        'detected_strains.bam',
        'detected.hisat2_map.log',
    )
    // STEP 7: Map same corrected reads to hisat2 index of T2T assembly.
    ch_human_bam = HISAT2_MAP_HUMAN(
        ch_reads.map { r -> record(sample_id: r.sample_id, reads_r1: r.reads_r1, reads_r2: r.reads_r2, index_files: human_index_files) },
        'T2T_mito_phix',
        '-k 1 --max-seeds 1000 --sp 1,1 --no-unal',
        'human_remap.bam',
        'hisat2_human.log',
    ).map { r -> record(sample_id: r.sample_id, human_bam: r.bam) }

    // STEP 8: Running featureCounts with accurate multimapper counting; annotate each read with gene, and number of bacterial and human mismatches ..
    ch_detected_fcounts = FEATURECOUNTS_DETECTED(
        ch_detected_map.join(ch_detected_ref.map { r -> record(sample_id: r.sample_id, gff: r.gff) }, by: 'sample_id').join(ch_paired, by: 'sample_id'),
        'detected',
    )
    // STEP 9: Annotate reads with bacterial gene, and mismatch rate vs bacterial and human genomes ..
    ch_detected_read_to_gene = ANNOTATE_MAPPED_READS_DETECTED(ch_detected_fcounts.join(ch_human_bam, by: 'sample_id'), 'detected')
    // STEP 10: Making a gene biotype table ..
    ch_gene_biotypes = EXTRACT_GENE_BIOTYPES(ch_detected_ref.map { r -> record(sample_id: r.sample_id, gff: r.gff) })
    // STEP 11: Filtering/annotating the featureCounts table and generating summary per RefSeq accession ..
    // Annotate the table and generate a summary per bacterial assembly (== accession)
    ch_detected_tables = MAKE_SUMMARY_TABLES_DETECTED(
        ch_detected_fcounts.join(ch_detected_read_to_gene, by: 'sample_id').join(ch_gene_biotypes, by: 'sample_id'),
        seq2taxid,
        'detected',
    ).map { r -> record(sample_id: r.sample_id, detected_annotated_fcounts_tsv: r.annotated_fcounts_tsv, detected_summary_tsv: r.summary_tsv, detected_protein_genes_list: r.protein_genes_list, detected_rrna_genes_list: r.rrna_genes_list) }

    // STEP 12: Filtering the accessions: keeping strains that pass our filtering criteria..
    // Filter the annotated table based on the number of rRNA and non-rRNA genes etc. If no strains passed the filtering, we quit, but we DO generate a mapstats.tsv.
    ch_filtered = FILTER_STRAINS(ch_detected_tables)
    ch_with_filtered_strains = ch_filtered.filter { r -> r.filtered_summary_tsv.size() > 0 }

    // STEP 13: Make the network file using the overlapping reads ..
    // Use fast(er) sort, since files can be big. Note that we ONLY use the filtered genes to estimate the overlap etc
    ch_shared_reads = GET_SHARED_UNIQ_READS(
        ch_with_filtered_strains.join(ch_detected_read_to_gene.map { r -> record(sample_id: r.sample_id, detected_read_to_gene_tsv: r.read_to_gene_tsv) }, by: 'sample_id')
    )
    ch_network_inputs = ch_shared_reads.join(ch_with_filtered_strains, by: 'sample_id')
    ch_network_all = MAKE_READ_NETWORK_ALL(
        ch_network_inputs.map { r -> record(sample_id: r.sample_id, reads_tsv: r.filtered_read_type_strain_tsv, filtered_annotated_fcounts_tsv: r.filtered_annotated_fcounts_tsv) },
        'nodes_and_edges.tsv',
    )
    ch_network_rrna = MAKE_READ_NETWORK_RRNA(
        ch_network_inputs.map { r -> record(sample_id: r.sample_id, reads_tsv: r.rrna_shared_reads_tsv, filtered_annotated_fcounts_tsv: r.filtered_annotated_fcounts_tsv) },
        'rRNA.nodes_and_edges.tsv',
    )
    ch_network_protein = MAKE_READ_NETWORK_PROTEIN(
        ch_network_inputs.map { r -> record(sample_id: r.sample_id, reads_tsv: r.protein_shared_reads_tsv, filtered_annotated_fcounts_tsv: r.filtered_annotated_fcounts_tsv) },
        'protein.nodes_and_edges.tsv',
    )
    // STEP 14: Parse the network and get the definitive list of strains ..
    // Cluster the read network into a more representative set of strains!
    ch_clusters = PARSE_READ_NETWORK(
        ch_with_filtered_strains.map { r -> record(sample_id: r.sample_id, filtered_summary_tsv: r.filtered_summary_tsv) }.join(ch_shared_reads.map { r -> record(sample_id: r.sample_id, rrna_uniq_counts_tsv: r.rrna_uniq_counts_tsv, protein_uniq_counts_tsv: r.protein_uniq_counts_tsv) }, by: 'sample_id').join(ch_network_all, by: 'sample_id').join(ch_network_rrna.map { r -> record(sample_id: r.sample_id, rrna_nodes_and_edges_tsv: r.nodes_and_edges_tsv) }, by: 'sample_id').join(ch_network_protein.map { r -> record(sample_id: r.sample_id, protein_nodes_and_edges_tsv: r.nodes_and_edges_tsv) }, by: 'sample_id'),
        gcf_species_genus,
    )

    // STEP 15: Making the new (top strain) bacterial reference and annotation, and building Hisat2 reference ..
    ch_top_ref = MAKE_COMBINED_REFERENCE_TOP(
        ch_clusters.map { r -> record(sample_id: r.sample_id, strains_list: r.top_strains_list) },
        ncbi_dir,
        'top',
    )
    ch_top_index = HISAT2_BUILD_TOP(ch_top_ref, 'top')
    // STEP 16: Mapping filtered unmapped reads to the top strain reference using Hisat2 ..
    ch_top_map = HISAT2_MAP_TOP(
        ch_reads.join(ch_top_index, by: 'sample_id'),
        'top_strains',
        '-k 100 --sp 1,1 --no-spliced-alignment --no-unal',
        'top_strains.bam',
        'top.hisat2_map.log',
    )
    // STEP 17: Running featureCounts on top strains with accurate multimapper counting, output BAM
    ch_top_fcounts = FEATURECOUNTS_TOP(
        ch_top_map.join(ch_top_ref.map { r -> record(sample_id: r.sample_id, gff: r.gff) }, by: 'sample_id').join(ch_paired, by: 'sample_id'),
        'top',
    )
    // STEP 18: Annotate reads with bacterial gene, and mismatch rate vs bacterial and human genomes ..
    ch_top_read_to_gene = ANNOTATE_MAPPED_READS_TOP(ch_top_fcounts.join(ch_human_bam, by: 'sample_id'), 'top')
    // STEP 19: Make annotated fcounts and summary table for top strains ..
    ch_top_initial_tables = MAKE_SUMMARY_TABLES_TOP_INITIAL(
        ch_top_fcounts.join(ch_top_read_to_gene, by: 'sample_id').join(ch_gene_biotypes, by: 'sample_id'),
        seq2taxid,
        'top',
    )

    // STEP 20: Make an annotation file of top strains without the bad genes, and with GCF (strain RefSeq) ID ..
    // Top strains need to be filtered again, in order for all counts to be correct
    ch_top_filtered_gff = FILTER_TOP_GFF(
        ch_top_initial_tables.map { r -> record(sample_id: r.sample_id, top_annotated_fcounts_tsv: r.annotated_fcounts_tsv) }.join(ch_top_ref.map { r -> record(sample_id: r.sample_id, top_gff: r.gff) }, by: 'sample_id').join(ch_with_filtered_strains.map { r -> record(sample_id: r.sample_id, blacklisted_genes_list: r.blacklisted_genes_list) }, by: 'sample_id'),
        seq2taxid,
    )
    // STEP 21: Redo the counting with a new GFF ..
    // gene: counts per locus_tag. strain: counting needed for UMI-tools; we will count reads with same UMI, barcode and strain once. This avoids overcounting due to rRNA operons
    ch_top_recount_inputs = ch_top_map
        .join(ch_top_filtered_gff.map { r -> record(sample_id: r.sample_id, gff: r.top_filtered_gff) }, by: 'sample_id')
        .join(ch_paired, by: 'sample_id')
    ch_top_gene_counts = FEATURECOUNTS_SORT_GENE(ch_top_recount_inputs, 'gene')
    ch_top_strain_counts = FEATURECOUNTS_SORT_STRAIN(ch_top_recount_inputs, 'strain')
    // STEP 22: REMAKE annotated fcounts and summary table for top strains - now without human/mismatched/blacklisted genes..
    ch_top_tables = MAKE_SUMMARY_TABLES_TOP(
        ch_top_gene_counts.map { r -> record(sample_id: r.sample_id, fcounts_tsv: r.counts_tsv) }.join(ch_top_read_to_gene, by: 'sample_id').join(ch_gene_biotypes, by: 'sample_id'),
        seq2taxid,
        'top',
    )

    // STEP 23: Running final quantification using UMI-tools .. (single-cell mode only)
    // One BAM has actual bacterial genes, another has RefSeq strain ID for "gene" (XT)
    UMI_TOOLS_COUNT_GENE(ch_top_gene_counts.join(ch_10x, by: 'sample_id'), 'gene')
    ch_umi_strain = UMI_TOOLS_COUNT_STRAIN(ch_top_strain_counts.join(ch_10x, by: 'sample_id'), 'strain')
    // STEP 24: Generate barcode-to-species table from per-gene UMI-tools .. (single-cell mode only)
    ch_barcode_to_species = BARCODE_TO_SPECIES(
        ch_umi_strain.map { r -> record(sample_id: r.sample_id, strain_umi_counts_tsv_gz: r.umi_counts_tsv_gz) },
        gcf_species_genus,
    )

    // STEP 25: Calculating the mapping stats ..
    // Only for samples that ended legitimately (no detected strains, no filtered strains, or all steps done)
    ch_mapstats_inputs = ch_reads
        .map { r -> record(sample_id: r.sample_id, preprocess_log: r.preprocess_log) }
        .join(ch_kuniq.map { r -> record(sample_id: r.sample_id, kuniq_report_txt: r.kuniq_report_txt) }, by: 'sample_id')
        .join(ch_detected_strains.map { r -> record(sample_id: r.sample_id, detected_strains_list: r.detected_strains_list) }, by: 'sample_id')
        .join(ch_detected_map.map { r -> record(sample_id: r.sample_id, detected_hisat2_log: r.hisat2_log) }, by: 'sample_id', remainder: true)
        .join(ch_detected_tables.map { r -> record(sample_id: r.sample_id, detected_protein_genes_list: r.detected_protein_genes_list, detected_rrna_genes_list: r.detected_rrna_genes_list, detected_annotated_fcounts_tsv: r.detected_annotated_fcounts_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_detected_read_to_gene.map { r -> record(sample_id: r.sample_id, detected_read_to_gene_tsv: r.read_to_gene_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_detected_fcounts.map { r -> record(sample_id: r.sample_id, detected_fcounts_bam: r.fcounts_bam) }, by: 'sample_id', remainder: true)
        .join(ch_filtered.map { r -> record(sample_id: r.sample_id, filtered_strains_list: r.filtered_strains_list, filtered_protein_genes_list: r.filtered_protein_genes_list, filtered_rrna_genes_list: r.filtered_rrna_genes_list, filtered_annotated_fcounts_tsv: r.filtered_annotated_fcounts_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_shared_reads.map { r -> record(sample_id: r.sample_id, filtered_read_type_strain_tsv: r.filtered_read_type_strain_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_clusters.map { r -> record(sample_id: r.sample_id, top_strains_list: r.top_strains_list) }, by: 'sample_id', remainder: true)
        .join(ch_top_map.map { r -> record(sample_id: r.sample_id, top_hisat2_log: r.hisat2_log) }, by: 'sample_id', remainder: true)
        .join(ch_top_tables.map { r -> record(sample_id: r.sample_id, top_protein_genes_list: r.protein_genes_list, top_rrna_genes_list: r.rrna_genes_list) }, by: 'sample_id', remainder: true)
        .join(ch_top_gene_counts.map { r -> record(sample_id: r.sample_id, top_sorted_bam: r.sorted_bam) }, by: 'sample_id', remainder: true)
        .filter { r ->
            def no_detected_strains = r.detected_strains_list.size() == 0
            def no_filtered_strains = r.filtered_strains_list != null && r.filtered_strains_list.size() == 0
            def all_steps_done = r.top_protein_genes_list != null && r.top_sorted_bam != null
            no_detected_strains || no_filtered_strains || all_steps_done
        }
    ch_mapstats = MAPSTATS(ch_mapstats_inputs)

    // Results
    ch_results = ch_mapstats
        .join(ch_kuniq.map { r -> record(sample_id: r.sample_id, kuniq_report_txt: r.kuniq_report_txt) }, by: 'sample_id')
        .join(ch_detected_tables.map { r -> record(sample_id: r.sample_id, detected_summary_tsv: r.detected_summary_tsv, detected_annotated_fcounts_tsv: r.detected_annotated_fcounts_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_filtered.map { r -> record(sample_id: r.sample_id, filtered_summary_tsv: r.filtered_summary_tsv, filtered_annotated_fcounts_tsv: r.filtered_annotated_fcounts_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_clusters.map { r -> record(sample_id: r.sample_id, filtered_cluster_tsv: r.filtered_cluster_tsv, top_cluster_tsv: r.top_cluster_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_top_tables.map { r -> record(sample_id: r.sample_id, top_summary_tsv: r.summary_tsv, top_annotated_fcounts_tsv: r.annotated_fcounts_tsv) }, by: 'sample_id', remainder: true)
        .join(ch_barcode_to_species, by: 'sample_id', remainder: true)
        .join(samples.map { r -> record(sample_id: r.sample_id, instrument: r.instrument) }, by: 'sample_id')
        .filter { r -> r.instrument != '10x' || r.top_summary_tsv == null || r.barcode_to_species_tsv != null }

    emit:
    results: Channel<DreamcatcherResult> = ch_results
}
