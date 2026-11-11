# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow (DSL2, static types) port of the `dreamcatcher` bash pipeline, which finds and quantifies bacterial reads in the host-unmapped reads of bulk and 10x RNA-seq. The original ran every tool in one container and one job, sized for KrakenUniq's memory and hisat2's cores together. Here every step is a process with its own resources, and the processes share four images grouped by the cost of re-running them (`-resume` is invalidated by an image change): `krakenuniq`, `mapping` (HISAT2, featureCounts), `r`, and `utils` for the quick steps (perl, samtools, Bowtie2, UMI-tools, whitelists). It runs standalone (`main.nf`) and is meant to be included by `atlas-ingest-pipeline` as the `DREAMCATCHER` subworkflow. See `README.md` for the algorithm and the mapping from the original steps to processes.

The reference data (~1 TB) and the CPUs/memory it needs are only on the HPC, so full runs are tested there by hand. Locally only `-stub-run` and the untracked `tests/replay` can run.

## Commands

```sh
# Lint/type-check every script (no other build step)
nextflow lint .

# Stub-run the whole pipeline: every process has a stub: block, local executor, no containers or reference data.
# tests/stub/ holds a 10x, a paired-end bulk and a single-end bulk sample plus empty reference fixtures.
nextflow run main.nf -profile stub -stub-run

# Real run on the QIMRB cluster
nextflow run main.nf -profile qimrb,pbspro,trace --input samplesheets/samples.csv -resume   # see run_nextflow.sh

# Replay processes on the intermediate files of an example output (needs docker with the images tagged
# ghcr.io/barrierimmunitylab/dreamcatcher-<name>:1.0.0, e.g. with docker tag), then compare with the example.
# tests/replay is git-ignored: local testing during the rewrite only
nextflow run tests/replay -c tests/replay/nextflow.config --example_dir example-outputs/dreamcatcher-bulk/SRX24417731 \
    --dreamcatcher_refdata_root <dir containing KrakenUniq_select/>
tests/replay/compare.sh local_dev/replay/results/SRX24417731 example-outputs/dreamcatcher-bulk/SRX24417731

# Build and test an image locally (CI of the containers repo does the same)
docker build --platform linux/amd64 -t local/dreamcatcher-utils:1.0.0 images/dreamcatcher-utils/1.0.0
docker run --rm --platform linux/amd64 -v "$PWD/images/dreamcatcher-utils/1.0.0:/t:ro" local/dreamcatcher-utils:1.0.0 bash /t/test.sh
```

## Code style

Same conventions as `atlas-ingest-pipeline`: `nextflow.enable.types = true` at the top of every script; data flows through channels as records (`record(sample_id: ..., ...)`), process inputs use the destructured form `record(sample_id: String, bam: Path)`, each process has one `record(...)` output, call sites assign the result and read fields by name, channels are combined with `join(..., by: 'sample_id')`, and record types shared across scripts live in `lib/types.nf`. Run `nextflow lint .` after every change. Prefer fewer comments; in workflow bodies a one-line comment names each channel-plus-process block, and a channel chain is grouped with the process call that consumes it.

Field names must be unique across records that are joined, so give output fields descriptive names (`detected_summary_tsv`, not `summary`) or project them with `map` before joining. A process called from several places is included under aliases (`include { HISAT2_MAP as HISAT2_MAP_TOP }`) and is parameterised through extra inputs (`prefix: String`), never through `ext` in the config.

## Architecture

- `main.nf` reads the samplesheet (`sample_id,instrument,reads_dir`), calls `DREAMCATCHER`, and publishes the key files per sample plus `versions.csv` with workflow outputs (`publish:` / `output {}`). `outputDir` must be set in the config and in each profile, since config expressions are evaluated before profile params apply.
- `subworkflows/local/dreamcatcher.nf` is the pipeline; `preprocess_reads.nf` is the 10x/bulk preprocessing it calls. Their processes are in `modules/local/dreamcatcher/<process>/` and `modules/local/preprocess_reads/<process>/`, one folder per process with `main.nf` and `resources/usr/bin/` for the scripts only that process uses. `nextflow.enable.moduleBinaries = true` puts `resources/usr/bin` on `PATH`, so a script is called by name; a script shared by two processes would have to be duplicated, so keep each script with the one process that calls it. `blacklisted_genes.list` is in `FILTER_STRAINS/resources/` and passed on to `FILTER_TOP_GFF` as an output.
- The scripts are the original perl/bash/R scripts. Edits were limited to removing the `$CMD` container wrapper, the BAM-input branches, and hard-coded paths. They keep their original file names, which Nextflow stages so that scripts reading fixed names from the current directory (`parse_read_network.R`, `calculate_mapstats.sh`) keep working. `MAPSTATS` stages its inputs into the original numbered folders (`stageAs ..., '2_detected_strains/*'`) for the same reason.
- Early exits: no detected strains, or none passing the filter, are `filter { ... .size() > 0 }` branches. `MAPSTATS` takes all later outputs as optional inputs (`Path?`, joined with `remainder: true`) and always runs, so those samples finish with a zero-filled `mapstats.tsv`.
- Resources live in `conf/pbspro.config` by process name (aliases get their own names). KrakenUniq's `--preload-size` is 85% of `task.memory`; the memory-bound sorts of `GET_SHARED_UNIQ_READS` use the same formula. Every process has `label "use_singularity"` (the `qimrb` profile loads the module for it) and `container "ghcr.io/barrierimmunitylab/dreamcatcher-<name>:<version>"` (the full URL, not a parameter).
- Every process declares the versions of the tools it runs in a `topic:` section (helpers in `lib/versions.nf`); the probes run inside the container.
- Reference data is passed as `Path` inputs and read through symlinks, so the Singularity bind mounts (`-B /working` in `qimrb`) must cover both the work directory and the reference data. The cell barcode whitelists are bundled in the `dreamcatcher-utils` image (`/opt/whitelists`, `$WHITELISTS_DIR`) and looked up by the file name recorded in STARsolo's `Log.out`, not by its absolute path.
- Images are in `images/<name>/<version>/{Dockerfile,test.sh}`, in the layout of `barrier-immunity-lab-containers` (they have not been copied there yet). Never edit a published version folder; make a new one and change the version in the processes.

## What "identical output" means

Compared with the original pipeline's outputs, files match line for line except: (1) several perl scripts (`kuniq_to_GCF.pl`, `make_summary_tables.pl`, `filter_strains.pl`, `make_read_network.pl`) print in hash order, which is different on every run, also in the original, so those files match after sorting; (2) `parse_read_network.R` uses Leiden clustering, which is random. `tests/replay/compare.sh` allows (1).
