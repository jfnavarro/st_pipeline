The ST Pipeline has been optimized for speed, robustness and it is very easy
to use with many parameters to adjust all the settings.
The ST Pipeline is fully parallel and has constant memory use.
The ST Pipeline allows to skip any of the steps and multiple customization options.
The ST Pipeline allows to use the genome or the transcriptome as reference.

The following files/parameters are commonly required :
- FASTQ files (Read 1 containing the spatial information and the UMI and read 2 containing the genomic sequence)
- A genome index generated with STAR
- An annotation file in GTF or GFF3 format (optional when using a transcriptome)
- The file containing the barcodes and array coordinates (look at the folder "ids" to use it as a reference).
Basically this file contains 3 columns (BARCODE, X and Y), so if you provide this
file with barcodes identifying cells (for example), the ST pipeline can be used for single cell data.
This file is also optional if the data is not barcoded (for example RNA-Seq data).
- A name for the dataset

The ST pipeline has multiple parameters mostly related to trimming, mapping and annotation
but generally the default values are good enough. You can see a full
description of the parameters typing "st_pipeline_run.py --help" after you have installed the ST pipeline.

The input FASTQ files can be given in gzip/bzip format as well.

Basically what the ST pipeline does (default mode) is :
- Quality trimming (read 1 and read 2):
  - Remove low quality bases
  - Sanity check (reads same length, reads order, etc..)
  - Check quality UMI
  - Remove artifacts (PolyT, PolyA, PolyG, PolyN and PolyC) of user defined length
  - Check for AT and GC content
  - Discard reads with a minimum number of bases of that failed any of the checks above
- Contamimant filter e.x. rRNA genome (Optional)
- Mapping with STAR (only read 2)
- Demultiplexing with [Taggd](https://github.com/jfnavarro/taggd) (only read 1)
- Keep reads (read 2) that contain a valid barcode and are correctly mapped
- Annotate the reads with htseq-count (slightly modified version)
- Group annotated reads by barcode (spot position), gene and genomic location (with an offset) to get a read count
- In the grouping/counting only unique molecules (UMIs) are kept.

You can see a graphical more detailed description of the workflow in the documents `workflow.pdf` and `workflow_extended.pdf`

The output dataset is a matrix of counts (genes as columns, spots as rows) in TSV format.
The ST pipeline will also output a log file with useful stats and information.
