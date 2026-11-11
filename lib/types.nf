nextflow.enable.types = true

record DreamcatcherInput {
    sample_id: String
    instrument: String
    unmapped_reads_with_logs: Iterable<Path>
}

def toDreamcatcherInput(row: Map<String, String>) {
    return record(
        sample_id: row.sample_id,
        instrument: row.instrument,
        unmapped_reads_with_logs: files("${row.reads_dir}/{Unmapped.out.mate[12]*,Log.out,Log.final.out}"),
    )
}

record DreamcatcherResult {
    sample_id: String
    mapstats_tsv: Path
    kuniq_report_txt: Path
    detected_summary_tsv: Path?
    detected_annotated_fcounts_tsv: Path?
    filtered_summary_tsv: Path?
    filtered_annotated_fcounts_tsv: Path?
    filtered_cluster_tsv: Path?
    top_summary_tsv: Path?
    top_cluster_tsv: Path?
    top_annotated_fcounts_tsv: Path?
    barcode_to_species_tsv: Path?
}

record PreprocessedReads {
    sample_id: String
    reads_r1: Path
    reads_r2: Path?
    preprocess_log: Path
}
