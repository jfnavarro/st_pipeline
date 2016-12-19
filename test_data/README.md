This dataset was generated
from the publicly available raw FASTQ files
of the Mouse Olfatory Bulb Replicate number 7
from the publication http://science.sciencemag.org/content/353/6294/78


The data was analysed with the ST Pipeline 1.1.6
and with a STAR genome index generated from
the Mus Musculus Ensembl annotation version 86. 
The annotation file used Mus Musculus GenCode 25 vM11.
A contaminant genome STAR index was used generated
from the Ensembl non conding RNA Mus musculus version 86.
The IDs file used to demultiplex is the 1000L5.
 
The following settings were used:

st_pipeline_run.py \
  --output-folder OUTPUT \
  --ids id.txt \
  --ref-map path_to_genome_index \
  --ref-annotation path_to_annotation_file.gtf \
  --expName SOME_NAME \
  --remove-polyA 10 \
  --remove-polyT 10 \
  --remove-polyG 10 \
  --remove-polyC 10 \
  --htseq-no-ambiguous \
  --verbose \
  --mapping-threads 16 \
  --log-file OUTPUT_log.txt \
  --two-pass-mode \
  --umi-filter \
  --filter-AT-content 90 \
  --filter-GC-content 90 \
  --contaminant-index path_to_contaminant_index \
  --min-length-qual-trimming 30 \
  --disable-clipping \
  R1.fastq R2.fastq
