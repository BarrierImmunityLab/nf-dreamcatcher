nextflow.enable.types = true

include { toDreamcatcherInput } from './lib/types.nf'
include { DREAMCATCHER } from './subworkflows/local/dreamcatcher'

workflow {

    main:
    ch_samples = channel.fromPath(params.input)
        .flatMap { csv -> csv.splitCsv(header: true) }
        .map { row -> toDreamcatcherInput(row) }
    ch_results = DREAMCATCHER(ch_samples)

    // Tool versions
    ch_versions = channel.topic('versions')
        .unique()
        .map { r -> "${r.process},${r.tool},\"${r.version}\"" }
        .mix(channel.of("pipeline,nextflow,\"${nextflow.version}\""))
        .collectFile(name: 'versions.csv', seed: 'process,tool,version', newLine: true, sort: true)

    publish:
    results = ch_results
    versions = ch_versions
}

output {
    results {
        path { r ->
            r.mapstats_tsv >> "${r.sample_id}/"
            r.kuniq_report_txt >> "${r.sample_id}/"
            r.detected_summary_tsv >> "${r.sample_id}/"
            r.detected_annotated_fcounts_tsv >> "${r.sample_id}/"
            r.filtered_summary_tsv >> "${r.sample_id}/"
            r.filtered_annotated_fcounts_tsv >> "${r.sample_id}/"
            r.filtered_cluster_tsv >> "${r.sample_id}/"
            r.top_summary_tsv >> "${r.sample_id}/"
            r.top_cluster_tsv >> "${r.sample_id}/"
            r.top_annotated_fcounts_tsv >> "${r.sample_id}/"
            r.barcode_to_species_tsv >> "${r.sample_id}/"
        }
    }
    versions {
        path '.'
    }
}
