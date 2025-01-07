Introduction
------------

The ST Pipeline contains the tools and scripts needed to process 
and analyze the raw files generated with the Spatial Transcriptomics 
or Visium in FASTQ format to generate datasets for down-stream analysis. 
The ST pipeline can also be used to process single cell RNA-seq data as 
long as a file with barcodes identifying each cell is provided.
The ST Pipeline can also process RNA-Seq datasets generated with 
or without UMIs. 

The ST Pipeline has been optimized for speed, robustness and 
it is very easy to use with many parameters to adjust all the settings.
The ST Pipeline is fully parallel and has constant memory use. 
The ST Pipeline allows to skip any of the steps and to use the 
genome or the transcriptome as reference. 

The following files/parameters are required:

- FASTQ files (Read 1 containing the spatial information and the UMI 
  and read 2 containing the genomic sequence) 
- A genome index generated with STAR 
- An annotation file in GTF or GFF format (optional)
- The file containing the barcodes and array coordinates 
   (look at the folder "ids" and chose the correct one). 
   Basically this file contains 3 columns (BARCODE, X and Y), 
   so if you provide this file with barcodes identinfying cells (for example), 
   the ST pipeline can be used for single cell data.
   This file is optional too. 
- A name for the dataset

The ST pipeline has multiple parameters mostly related to trimming, 
mapping and annotation but generally the default values are good enough. 
You can see a full description of the parameters 
typing "st_pipeline_run.py --help" after you have installed the ST pipeline.

The input FASTQ files can be given in gzip/bzip format as well. 

Basically what the ST pipeline does is:

- Quality trimming (read 1 and read 2):
    - Remove low quality bases
    - Sanity check (reads same length, reads order, etc..)
    - Check quality UMI (if provided)
    - Remove artifacts (PolyT, PolyA, PolyG, PolyN and PolyC) of user defined length
    - Check for AT and GC content
    - Discard reads with a minimum number of bases of that failed any of the checks above
- Contamimant filter e.x. rRNA genome (Optional)
- Mapping with STAR (only read 2)
- Demultiplexing with [Taggd](https://github.com/jfnavarro/taggd) (only read 1)
- Keep reads (read 2) that contain a valid barcode and are correctly mapped
- Annotate the reads with htseq-count (optional)
- Group annotated reads by barcode(spot position) and gene to get a read count
- In the grouping/counting only unique molecules (UMIs) are kept. 

You can see a graphical more detailed description of the workflow in the documents workflow.pdf and workflow_extended.pdf

The output will be a matrix of counts (genes as columns, spots as rows)
and a log file with useful information and stats.