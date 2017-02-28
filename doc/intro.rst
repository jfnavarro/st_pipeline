Introduction
------------

ST Pipeline contains the tools and scripts needed to process and analyze the raw
files generated with the Spatial Transcriptomics method in FASTQ format.
The ST pipeline can also be used to process single cell data as long as
a file with barcodes identifying each cell is provided.

The following files/parameters are required :

    * FASTQ files (Read 1 containing the BARCODE and the UMI and read
      2 containing the genomic sequence)
    * A genome index generated with STAR
    * An annotation file in GTF or GFF format
    * The file containing the barcodes and array coordinates (look at the folder
      "ids" and chose the correct one.
    * A name for the dataset

The ST pipeline has multiple parameters mostly related to trimming, mapping and
annotation but generally the default values are good enough. You can see a full
description of the parameters typing "st_pipeline_run.py --help" after you have
installed the ST pipeline.

The raw data can be given in gzip format as well.

Basically what the ST pipeline does is :

    * Quality trimming (read 1 and read 2) :
        * Remove low quality bases
        * Sanity check (reads same length, reads order, etc..)
        * Check quality UMI
        * Remove artifacts (PolyT, PolyA, PolyG and PolyC)
        * Check for AT and GC content
        * Discard reads with a minimum number of bases of that failed any of the
          checks above
    * Contamimant filter e.x. rRNA genome (Optional)
    * Mapping with STAR (only read 2)
    * Demultiplexing with taggd
      (https://github.com/SpatialTranscriptomicsResearch/taggd) (only read 1)
    * Keep reads (read 2) that contain a valid barcode and are correctly mapped
    * Annotate the reads with htseq-count
    * Group annotated reads by barcode(spot position) and gene to get a read count
    * In the grouping/counting only unique molecules (UMIs) are kept.

You can see a graphical more detailed description of the workflow in the
documents workflow.pdf and workflow_extended.pdf

The output will be a data frame file with the counts (genes as columns, spots as
rows), a BED file containing the transcripts (Read name, coordinate, gene,
etc..), and a JSON file with useful stats. The ST pipeline will also output a log
file with useful information.
