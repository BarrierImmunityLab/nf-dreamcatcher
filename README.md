# Dreamcatcher

<img src="https://github.com/apredeus/dreamcatcher/blob/main/img/dreamcatcher.png">

Dreamcatcher is an advanced pipeline for accurate identification of _bacterial_ reads in human RNA-seq data. It has two modes: 

  - bulk, in which it identifies bacterial species present in the sample and quantifies their abundance & gene expression; 
  - single cell, which includes several extra steps which, in addition to the "bulk" output, result in a (barcode) x (bacterial species) matrix, which could be conveniently added to the single cell/spatial metadata. 

Unfortunately, only 10x protocols are currently supported on a single-cell level. You can run Smart-seq(1/2/3) or other plate-based data in a bulk mode though. 

## Motivation

Main motivation behind `dreamcatcher` is as follows: most metagenome profiling methods can be roughtly divided into k-mer based and mapping-based. Mapping-based methods are usually considered to be too slow to run against a modern bacterial database, containing 10's of thousands of genomes. At the same time, k-mer based methods produce overwhelming number of false positives, which come from two sources: (1) spurious matches to foreign sequences in the reference genomes; (2) huge sampling bias towards certain bacterial species that are human pathogens. Both problems are exacerbated in case of RNA-seq, since most (~90%) of the detected reads come from bacterial rRNA operons. 

At the same time, careful case-by-case consideration of bacteria of interest has shown that it's possible to distinguish spurious and proper alignments when using mapped reads and accurate genome annotation. Furthermore, when only minimal filtering is applied, k-mer-based algorithms do not suffer from substantial false negative rates; simply put, if the bug is present in the reference and in the sample, they would detect it. We decided to use a combination of these methods, and further leverage gene annotation to reliably define the bacterial species truly present in the sample. 

## Algorithm

`dreamcatcher` leverages RefSeq genomes with matched curated genome annotations (see section on the Reference preparation below). Our current reference follows regular Kraken/Kraken2/KrakenUniq approach. This allows us to put a comprehensive RefSeq-based 

  - reads are mapped using STAR, STARsolo, Cell Ranger, or a similar tool. 
  - in case of 10x single cell/Visium protocol, unmapped reads are converted into a "UMI-tools" style of formatting, with barcode and UMI sequences added to the read name; 
  - KrakenUniq is used to classify the unmapped reads, and assemblies with at least minimal number of unique k-mers are selected (usually 200-1000 strains). These are the _detected strains_;
  - Hisat2 is used to map the same unmapped reads to this selection of genomes, and featureCounts is used to find the expressed genes; 
  - expression table is then annotated according to the gene type (rRNA or not), mean number of mismatches per 100 bp, etc;
  - following this, unmapped reads are mapped to T2T human assembly as well as PhiX; 
  - genes that have too many mismatches, that map to human assembly or PhiX, or genes from the blacklist are removed; 
  - strains which do not express a minimal number of genes are removed as contaminants/artifacts. What remains is referred to as _filtered strains_;
  - following this, a graph of shared reads is generated for the filtered strains; 
  - the graph is clustered iteratively, based on shared read topology and taxonomy;
  - this results in selection of strains named _top strains_; this is the main result of the "bulk" part of the pipeline;
  - in case of single-cell mode, the unmapped reads are then mapped to the genomes of _top strains_, and the resulting BAM file is counted using `featureCounts`; 
  - using `UMI-tools` allows us to get UMI counts per barcode, per _top strains_ gene; 
  - finally, we sum the expression up per species and generate a final (barcode) x (species) table. 

## Output 

In the original bash pipeline all files of a sample were written to one directory with the following layout. The Nextflow pipeline keeps every one of these files as a declared output of the process that makes it (see `modules/local/dreamcatcher/*/main.nf`), so they are all available in the Nextflow work directory, but only the key files listed further down are published. The numbered folders are:

  - `0_preprocessed_reads` - with read preprocessing logs and barcode corrections (if performed); 
  - `1_krakenuniq` - the output of `krakenuniq` unmapped read classification; 
  - `2_detected_strains` - list of the detected strains, as well as `hisat2` indexes, BAM files, annotated `featureCounts` outputs, and per-read annotations;
  - `3_filtered_strains` - list of the filtered strains, as well as annotated per-gene and per-read outputs; 
  - `4_read_networks` - all the files that are the input and the output of the read graph collapsing algorithm; 
  - `5_top_strains` - list of the top strains, as well as `hisat2` indexes, BAM files, annotated `featureCounts` outputs; 
  - `6_umitools` - results of `umi-tools` processing, in case of 10x single cell or Visium experiments. Not present for bulk samples.  

The pipeline publishes these key files for each sample (`<outdir>/<sample_id>/`):
  
  - `kuniq.report.txt` - KrakenUniq report file that contains metagenomic classification for all reads unmapped to host, and can be visualised using [Pavian](https://github.com/fbreitwieser/pavian); 
  - `detected.annotated_fcounts.tsv`, `detected.summary.tsv` - per-gene and per-strain metrics for the detected strains, derived by mapping reads to a set of strains (`detected_strains.list`) with > 3 unique k-mers in KrakenUniq;  
  - `filtered.annotated_fcounts.tsv`, `filtered.summary.tsv` - per-gene and per-strain metrics for the filtered strains. Details of the filtering process could be found in `3_filtered_strains/filter_strains.log`; 
  - `filtered.cluster.tsv`, `top.cluster.tsv` - results of read network clustering (input and output strains);
  - `mapstats.tsv` - detailed statistics of reads that were mapped/assigned to a gene in detected, filtered, and top strains; 
  - `barcode_to_species.tsv` - corrected barcode by identified bacterial species matrix. Not present for bulk samples.

Two file types that don't have a header are `*.summary.tsv` and `*.annotated_fcounts.tsv`. Given below is the description of each file. 

Files `detected.annotated_fcounts.tsv`, `filtered.annotated_fcounts.tsv`, and `top.annotated_fcounts.tsv` all are made using `featureCounts` output, and extra annotation for each gene was added to them. There are 14 columns overall, and following table explains each one using an example output: 

<div align="center">

| RefSeq locus tag | Contig | Gene start | Gene end | Strand | Gene length | Split read count | Raw read count | Gene biotype | Mismatch rate (per 100bp) | Fraction T2T-human/phiX | RefSeq strain ID | Strain name | Strain taxid | 
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
P9N54_RS23305 | NZ_CP121294.1 | 4776927 | 4778102 | - | 1176 | 0.142853 | 1 | protein_coding | 0.662 | 0.000 | GCF_029674705.1 | Escherichia coli O155:H21 strain=NWU_1 | 3038394 |
P9N54_RS23360 | NZ_CP121294.1 | 4788469 | 4790003 | - | 1535 | 20.634720 | 6231 | rRNA | 0.516 | 0.000 | GCF_029674705.1 | Escherichia coli O155:H21 strain=NWU_1 | 3038394 |
NX106_RS23395 | NZ_JAMSJK010000099.1 | 1 | 321 | + | 321 | 0.125000 | 1 | rRNA | 2.532 | 1.000 | GCF_024733345.1 | Escherichia ruysiae strain=C61-1 | 2608867 |
NX106_RS23520 | NZ_JAMSJK010000126.1 | 194 | 280 | + | 87 | 0.003418 | 4 | rRNA | 0.625 | 0.000 | GCF_024733345.1 | Escherichia ruysiae strain=C61-1 | 2608867 |

</div>

Files `detected.summary.tsv`, `filtered.summary.tsv`, and `top.summary.tsv` summarise all the strains in the respective groups based on the detected genes. There are 13 columns that look as described below. Note that "protein" gene in the summary table is anything that is not rRNA in RefSeq `gene_biotype`; thus, noncoding RNAs, tRNAs, and other miscellaneus gene types are sometimes (rarely) included here.  

<div align="center">

| RefSeq strain ID | Strain name | Strain taxid | Detected rRNA genes | Split read sum - rRNA | Raw read sum - rRNA | Mean rRNA mismatch rate (per 100bp) | Mean rRNA fraction T2T-human/phiX | Detected protein genes | Split read sum - protein | Raw read sum - protein | Mean protein mismatch rate (per 100bp) | Mean protein fraction T2T-human/phiX |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
GCF_014083985.1 | Streptomyces sp. SCUT-3 | 2684469 | 14 | 0.061646 | 194 | 0.000 | 0.000 | 0 | 0.000000 | 0 | 0.000 | 0.000 |
GCF_002861525.1 | Actinomyces urogenitalis strain=UMB0319 | 103621 | 2 | 0.349945 | 54 | 0.039 | 0.000 | 3 | 2.653839 | 10 | 0.403 | 0.000 | 
GCF_014640315.1 | Hymenobacter glacieicola strain=CGMCC 1.12990 | 1562124 | 1 | 0.000397 | 4 | 0.000 | 0.000 | 2 | 1.022217 | 3 | 1.686 | 0.333 |
GCF_029674705.1 | Escherichia coli O155:H21 strain=NWU_1 | 3038394 | 3 | 59.922806 | 10997 | 0.346 | 0.000 | 24 | 21.642791 | 68 | 0.796 | 0.000 |
 
</div>

Mean number of mismatches and human/phiX contamination is used for filtering. Thus, no strains in the `filtered` or `top` cohort should have genes with mean mismatches over 3.4 (for rRNA)/2.3 (for protein), or with more than 0.100 human/phiX overlap. 

## Reference preparation 

Dreamcatcher reference takes a while to make. The list of used RefSeq strains, together with the assembly metadata, can be found in this repository (`/data/assembly_summary_bacteria.txt.gz`). Two of the biggest items in the referense are (1) download of RefSeq genomes and annotations; (2) creation of the KrakenUniq database; (3) creation of T2T+phiX decoy `hisat2` index. 

To get RefSeq genomes and annotations in GFF3 format, you would need to install [NCBI datasets and dataformat utilities](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/). Following this, take the strain list (`/data/refseq_strain.list` in this repo) and run

```bash
datasets download genome accession --include genome,gff3,gtf --dehydrated --inputfile refseq_strain.list
unzip ncbi_dataset.zip
datasets rehydrate --gzip  --directory .
```
This will create a directory with subdirectories named `GCF_007970465.1` etc; each individual subdirectory should contain the gzipped genome file (ending with `genomic.fna.gz`) and GFF3 file (`genomic.gff.gz`). Make sure all strains have these two files downloaded successfully before moving on. The current database has 51,800 bacterial species, and takes up about 200 Gb of space.  

To create `krakenuniq` database, we have used the following commands. To re-create the exact database used in the paper, the strain list above needs to be used. The command below used 32 cores and limited the database t o 480 Gb for sensitivity. The final folder needs to be cleaned up and takes up about 600 Gb. 

```bash
krakenuniq-download --threads 32 --db . 'refseq/vertebrate_mammalian/Any/species_taxid=9606'
krakenuniq-download --threads 32 --db . refseq/viral/Any viral-neighbors
krakenuniq-download --threads 32 --db . --dust refseq/archaea
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Contig/refseq_category='representative genome'
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Contig/refseq_category='reference genome'
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Scaffold/refseq_category='representative genome'
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Scaffold/refseq_category='reference genome'

krakenuniq-build --jellyfish-bin /software/cellgeni/miniconda3/bin/jellyfish --max-db-size 480 --db . --kmer-len 31 --threads 32 --taxids-for-genomes --taxids-for-sequences
```
This takes a long time to create (approximately 2 weeks on 32 cores). Please create a free [Globus](https://www.globus.org/) account and let me know the email you've registered with if you need to download the pre-made reference. 

Additional files required by `dreamcatcher` are: 

  - `seqid2taxid.map.refseq` - full mapping of all sequences available in our KrakenUniq database. These are generated by `krakenuniq` when making the reference; 
  - `assembly_summary_bacteria.txt` - GenBank metadata of all used RefSeq bacterial assemblies, including assembly level and strain info; 
  - `gcf_species_genus.tsv` - mapping of all bacterial RefSeq IDs to NCBI species/genus names and TaxIDs. 

All of these smaller files are available as archives in `/data` subdirectory of this repo. 

## Pipeline organization 

The pipeline is written in [Nextflow](https://www.nextflow.io/) (DSL2, static types), with `perl` used extensively for text file parsing. Each step is its own process with its own resource request: for example, `KRAKENUNIQ` asks for a lot of memory but little else, while the `hisat2` steps ask for many cores but little memory. There are four container images, grouped by how costly it is to re-run the processes that use them: `dreamcatcher-krakenuniq`, `dreamcatcher-mapping` (HISAT2 and featureCounts), `dreamcatcher-r`, and `dreamcatcher-utils` for the quick steps. The images are defined in `images/` (one folder per image, `images/<name>/<version>/Dockerfile`) and are published as `ghcr.io/barrierimmunitylab/<name>:<version>`.

```
PREPROCESS_READS ─┬─ KRAKENUNIQ ─ KUNIQ_TO_GCF ─ detected strains: MAKE_COMBINED_REFERENCE ─ HISAT2_BUILD ─ HISAT2_MAP ─ FEATURECOUNTS ─ ANNOTATE_MAPPED_READS ─ MAKE_SUMMARY_TABLES ─┐
                  └─ HISAT2_MAP (T2T + phiX) ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
─ FILTER_STRAINS ─ GET_SHARED_UNIQ_READS ─ MAKE_READ_NETWORK (x3) ─ PARSE_READ_NETWORK ─ top strains: (same mapping steps, then FILTER_TOP_GFF and a recount) ─ UMI_TOOLS_COUNT (10x only) ─ BARCODE_TO_SPECIES ─ MAPSTATS
```

The numbered steps of the original `dreamcatcher` script correspond to processes as follows.

| Original step | Process(es) |
|:-|:-|
| 1. Extract and preprocess unmapped reads | `MAKE_UMI_FASTQ`, `GET_BARCODE_CORRECTIONS`, `CORRECT_BARCODES` (10x) or `MAKE_BULK_FASTQ` (bulk) |
| 2. KrakenUniq | `KRAKENUNIQ` |
| 3. Detected strains | `KUNIQ_TO_GCF` |
| 4-5. Combined reference and `hisat2` index | `MAKE_COMBINED_REFERENCE_DETECTED`, `HISAT2_BUILD_DETECTED` |
| 6-8. Mapping to the detected strains and to T2T+phiX, counting | `HISAT2_MAP_DETECTED`, `HISAT2_MAP_HUMAN`, `FEATURECOUNTS_DETECTED` |
| 9-11. Per-read, per-gene and per-strain tables | `ANNOTATE_MAPPED_READS_DETECTED`, `EXTRACT_GENE_BIOTYPES`, `MAKE_SUMMARY_TABLES_DETECTED` |
| 12. Strain filtering | `FILTER_STRAINS` |
| 13. Read network | `GET_SHARED_UNIQ_READS`, `MAKE_READ_NETWORK_ALL/RRNA/PROTEIN` |
| 14. Clustering | `PARSE_READ_NETWORK` |
| 15-19. Top strains: reference, mapping, counting, tables | `MAKE_COMBINED_REFERENCE_TOP`, `HISAT2_BUILD_TOP`, `HISAT2_MAP_TOP`, `FEATURECOUNTS_TOP`, `ANNOTATE_MAPPED_READS_TOP`, `MAKE_SUMMARY_TABLES_TOP_INITIAL` |
| 20-22. Remove bad genes, recount | `FILTER_TOP_GFF`, `FEATURECOUNTS_SORT_GENE`, `FEATURECOUNTS_SORT_STRAIN`, `MAKE_SUMMARY_TABLES_TOP` |
| 23-24. UMI counting, barcode-to-species table (10x) | `UMI_TOOLS_COUNT_GENE`, `UMI_TOOLS_COUNT_STRAIN`, `BARCODE_TO_SPECIES` |
| 25-26. Mapping statistics | `MAPSTATS` |

Differences from the original script:

  - input is the STAR/STARsolo output directory only (`Unmapped.out.mate1/2` with `Log.out` and `Log.final.out`); the `--bam` and `--skip-remap` options are gone (use `-resume` instead);
  - the human/phiX remapping does not depend on KrakenUniq, so it runs at the same time as KrakenUniq;
  - if no strain is detected, or none passes the filtering, the sample finishes successfully with a zero-filled `mapstats.tsv` (the original wrote it and then exited with an error);
  - the cell barcode whitelists are bundled in the `dreamcatcher-utils` image (from the same cellgeni/STARsolo commit as the `cellgeni-starsolo` image) and the whitelist is looked up there by the file name recorded in the STARsolo `Log.out`, since the recorded path is not visible inside this container; there is no `whitelists_dir` reference data.

### Running the analysis

Create a samplesheet with one row per sample. `instrument` is `10x` for single cell/spatial data; any other value processes the sample as bulk. `reads_dir` is the STAR/STARsolo output directory.

```
sample_id,instrument,reads_dir
SRX20761959,10x,/path/to/starsolo/SRX20761959
SRX24417731,bulk,/path/to/star/SRX24417731
```

and run

```bash
nextflow run main.nf -profile qimrb,pbspro,trace --input samples.csv -resume
```

The profiles `qimrb` (cluster paths, Singularity bind mounts), `pbspro` (executor and the per-process resources in `conf/pbspro.config`) and `trace` (execution trace with the peak memory of every task) are those of the [atlas-ingest-pipeline](https://github.com/BarrierImmunityLab/atlas-ingest-pipeline). Edit `conf/pbspro.config` after looking at the trace. The reference data root (`--dreamcatcher_refdata_root`) must contain what is described in the Reference preparation section:

```
<dreamcatcher_refdata_root>/NCBI_datasets/bacteria_ref/ncbi_dataset/data/   RefSeq genomes and annotations
<dreamcatcher_refdata_root>/KrakenUniq_select/                              KrakenUniq database, seqid2taxid.map.refseq, assembly_summary_bacteria.txt, gcf_species_genus.tsv
<dreamcatcher_refdata_root>/decoy_reference/T2T_mito_phix.*.ht2             T2T + phiX hisat2 index
```

Results are published to `<outputDir>/<sample_id>/`, together with `<outputDir>/versions.csv`.

### Running it from another pipeline

`subworkflows/local/dreamcatcher.nf` takes a channel of records `(sample_id, instrument, unmapped_reads_with_logs)` and emits the key files of every sample as a record (see `DreamcatcherResult` in `lib/types.nf`). It uses the parameters `dreamcatcher_refdata_root`. The calling pipeline needs `nextflow.enable.moduleBinaries = true` (the scripts of every process are in its `resources/usr/bin` folder), the Singularity bind mounts and `module load singularity` as in the `qimrb` profile, and the resources in `conf/pbspro.config`. The scripts and `blacklisted_genes.list` are read from the pipeline directory inside the containers, so that directory (like the work directory and the reference data) must be visible there: under `/working` for the `qimrb` bind mount, or under `$HOME`, which Singularity binds by default.

### Testing

  - `nextflow run main.nf -profile stub -stub-run` checks the wiring of all processes without any tools or reference data;
  - to compare a full run with an example output, copy the Nextflow work directory of the run and compare the files of the process directories by name with the numbered folders of the example.
