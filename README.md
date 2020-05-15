# Spatial Transcriptomics Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![PyPI version](https://badge.fury.io/py/stpipeline.svg)](https://badge.fury.io/py/stpipeline)
[![Build Status](https://travis-ci.org/jfnavarro/st_pipeline.svg?branch=master)](https://travis-ci.org/jfnavarro/st_pipeline)

The ST Pipeline contains the tools and scripts needed to process and analyze the raw files generated with the Spatial Transcriptomics and Visium in FASTQ format to generate datasets for down-stream analysis. 
The ST pipeline can also be used to process single cell RNA-seq data as long as a file with barcodes identifying each cell is provided (same template as the files in the folder "ids").

The ST Pipeline has been optimized for speed, robustness and it is very easy to use with many parameters to adjust all the settings.
The ST Pipeline is fully parallel and has constant memory use. 
The ST Pipeline allows to skip any of the steps and to use the genome or the transcriptome as reference. 

The following files/parameters are required :
- FASTQ files (Read 1 containing the spatial information and the UMI and read 2 containing the genomic sequence) 
- A genome index generated with STAR 
- An annotation file in GTF or GFF3 format (optional when using a transcriptome)
- The file containing the barcodes and array coordinates (look at the folder "ids" and chose the correct one). Basically this file contains 3 columns (BARCODE, X and Y), so if you provide this file with barcodes identinfying cells (for example), the ST pipeline can be used for single cell data. This file is also optional if the data is not barcode (for example RNA-Seq data).
- A name for the dataset

The ST pipeline has multiple parameters mostly related to trimming, mapping and annotation but generally the default values are good enough. You can see a full description of the parameters typing "st_pipeline_run.py --help" after you have installed the ST pipeline.

The input FASTQ files can be given in gzip/bzip format as well. 

Basically what the ST pipeline does is :
- Quality trimming (read 1 and read 2) :
	- Remove low quality bases
	- Sanity check (reads same length, reads order, etc..)
	- Check quality UMI
	- Remove artifacts (PolyT, PolyA, PolyG, PolyN and PolyC) of user defined length
	- Check for AT and GC content
	- Discard reads with a minimum number of bases of that failed any of the checks above
- Contamimant filter e.x. rRNA genome (Optional)
- Mapping with STAR (only read 2)
- Demultiplexing with [Taggd](https://github.com/SpatialTranscriptomicsResearch/taggd) (only read 1)
- Keep reads (read 2) that contain a valid barcode and are correctly mapped
- Annotate the reads with htseq-count (slightly modified version)
- Group annotated reads by barcode(spot position), gene and genomic location (with an offset) to get a read count
- In the grouping/counting only unique molecules (UMIs) are kept. 

You can see a graphical more detailed description of the workflow in the documents workflow.pdf and workflow_extended.pdf

The output will be a matrix of counts (genes as columns, spots as rows),
a BED file containing the transcripts (Read name, coordinate, gene, etc..), and a JSON
file with useful stats.
The ST pipeline will also output a log file with useful information.

**Installation**

We recommend you install a virtual environment like Pyenv or Anaconda before you install the pipeline. 

The ST Pipeline works with python 3.6 or bigger.

You can install the ST Pipeline using PyPy:

    pip install stpipeline
 
Alternatively, you can build the ST Pipeline yourself:

First clone the repository 

    git clone <stpipeline repository> 
    
or download a tar/zip from the releases section and unzip it

    unzip stpipeline_release.zip
    
Access the cloned ST Pipeline folder or the folder where the tar/zip file has been decompressed. 

    cd stpipeline

To install the pipeline type then

    python setup.py build
    python setup.py install

To run a test type (you need internet connection to run the tests)

    python setup.py test

To see the different options type 

    st_pipeline_run.py --help
    
**Requirements**

The ST Pipeline requires STAR installed in the system (minimum version 2.5.4 if you use a ST Pipeline version >= 1.6.0):
https://github.com/alexdobin/STAR

If you use anaconda you can install STAR with

    conda install -c bioconda star
    
The ST Pipeline requires samtools installed in the system
If you use anaconda you can install Samtools with

    conda install -c bioconda samtools openssl=1.0

The ST Pipeline recommends a computer with at least 32GB of RAM (depending on the size of the genome) and 8 cpu cores. 

**Dependencies** 

The ST Pipeline depends on some Python packages that will
be automatically installed during the installation process. 
You can see them in the file dependencies.txt

**Example**

An example run would be

	st_pipeline_run.py --expName test --ids ids_file.txt --ref-map path_to_index --log-file log_file.txt --output-folder /home/me/results --ref-annotation annotation_file.gtf file1.fastq file2.fastq 

**Emsembl ids**

If you used an Ensembl annotation file and you would like change
the ouput file so it contains gene Ids/names instead of Ensembl ids. 
You can use this tool that comes with the ST Pipeline

	convertEnsemblToNames.py --annotation path_to_annotation_file --output st_data_updated.tsv st_data.tsv
	
**Merge demultiplexed FASTQ files**

If you used different indexes to sequence and need to merge the files
you can use the script merge_fastq.py

	merge_fastq.py --run-path path_to_run_folder --out-path path_to_output --identifiers S1 S2 S3 S4
	
Where identifiers will be strings that identify each demultiplexed sample. 

**Filter out genes by gene type**

If you want to remove from the dataset (matrix in TSV) genes corresponding
to certain gene types (For instance to keep only protein_coding). You can do
so with the script filter_gene_type_matrix.py

	filter_gene_type_matrix.py --counts-matrix stdata.tsv --gene-types-keep protein-coding --outfile new_stdata.tsv --annotation path_to_annotation_file
	
You may include the parameter --ensembl-ids if your gene names are represented as gene ids instead.

**Remove spots from dataset**

If you want to remove spots from a dataset (matrix in TSV) for instance
to keep only spots inside the tissue. You can do so with the script adjust_matrix_coordinates.py

	adjust_matrix_coordinates.py --counts-matrix stadata.tsv --outfile new_stdata.tsv --coordinates-file coordinates.txt
	
Where coordinates.txt will be a tab delimited file with 6 columns:

orig_x orig_y new_x new_y new_pixel_x new_pixel_y

Only spots whose coordinates in the file will be kept and then optionally you
can update the coordinates in the matrix choosing for the new array or pixel coordinates.

**Quality stats**

The ST Pipeline generate useful statistical information in the LOG file but if you
want to obtain more detail information about the quality of the data, you can run the following script:

	st_qa.py --input-data stdata.tsv 
	
If you want to perform quality stats on multiple samples you can run:

	multi_qa.py --counts-table-files stdata1.tsv stadata2.tsv stdata3.tsv ... 
	
Multi_qa.py generates violing plots, correlation plots/tables and more useful information and 
it allows to log the counts for the correlation.
	
**Documentation**

You can see a more detailed documentation in the folder "doc_out".

**Example data**

You can see a real dataset obtained from the public data from
the following publication (http://science.sciencemag.org/content/353/6294/78)
in the folder called "data".

**License**

The ST pipeline is open source under the MIT license which means that you can use it, change it and re-distribute but you must always refer to our license (see LICENSE and AUTHORS).

**Reference**

If you use the ST Pipeline, please refer its publication: 

ST Pipeline: An automated pipeline for spatial mapping of unique transcripts
Oxford BioInformatics
10.1093/bioinformatics/btx211

**Contact**

For questions, bugs, feedback, etc.. you can contact 
Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


