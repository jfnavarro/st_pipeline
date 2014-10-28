*** Copyright Spatial Transcriptomics, 2014. ***


This folder contains test data for the pipeline.
Data was created from an old (MiSeq F1) original dataset for Mouse GRCM 38 as follows:

- A repeat-masked version of the shortest mouse chromosome, 19, was obtained from Ensembl.
- A Bowtie2-index, 'chromosome19', was built for the chromosome.
- All reads that were mapped against chromosome 19 *for the original dataset* were retained as new input.
- The HTSeq original annotation file was stripped of all other content than chromosome 19, and modified
  to play with the new Bowtie index.
- The chip file (130307_Design3_27mer.txt) was unaltered, as was the contaminant genome (R45S5_R5S1).


Example of a run (started with this folder as working directory):

st_pipeline_run.py ./input/miseqF1/testdata_R1.fastq ./input/miseqF1/testdata_R2.fastq --ids ./config/idfiles/130307_Design3_27mer.txt --ref-map ./config/genomes/mouse_grcm38/chromosome19 --ref-annotation ./config/annotations/mouse_grcm38_chromosome19.gtf --allowed-missed 5 --allowed-kimer 7 --min-length-qual-trimming 28 --mapping-fw-trimming 51 --mapping-rv-trimming 5 --length-id 27 --min-quality-trimming 20 --verbose --no-clean-up --bowtie-threads 3 --error-id 0 --start-id 0 --htseq-mode intersection-nonempty --contaminant-bowtie2-index ./config/contaminant_genomes/R45S5_R5S1/rnagenome --expName testdata --output-folder ./output --temp-folder ./tmp




