module load singularity/3.11.4

nextflow run main.nf \
    -profile qimrb,pbspro,trace \
    --input samplesheets/samples.csv \
    -resume
