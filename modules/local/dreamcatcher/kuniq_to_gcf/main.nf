nextflow.enable.types = true

include { perlVersion } from '../../../../lib/versions'

process KUNIQ_TO_GCF {
    tag sample_id
    label "use_singularity"
    container "ghcr.io/barrierimmunitylab/dreamcatcher-utils:1.0.0"

    input:
    record(
        sample_id: String,
        kuniq_report_txt: Path
    )
    gcf_species_genus: Path
    assembly_summary: Path

    output:
    record(
        sample_id: sample_id,
        detected_strains_list: file('detected_strains.list'),
        kuniq_to_gcf_log: file('kuniq_to_GCF.log'),
    )

    topic:
    record(process: task.process, tool: 'perl', version: eval(perlVersion())) >> 'versions'

    script:
    """
    kuniq_to_GCF.pl ${kuniq_report_txt} ${gcf_species_genus} ${assembly_summary} > detected_strains.list 2> kuniq_to_GCF.log
    """

    stub:
    """
    if [[ "${sample_id}" == *_NO_STRAINS ]]
    then
        touch detected_strains.list
    else
        echo GCF_000000001.1 > detected_strains.list
    fi
    touch kuniq_to_GCF.log
    """
}
