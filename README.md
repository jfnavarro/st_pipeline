# Spatial Transcriptomics Pipeline

[![Build Status](https://travis-ci.org/jfnavarro/st_pipeline.svg?branch=master)](https://travis-ci.org/jfnavarro/st_pipeline)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Coverage Status](https://coveralls.io/repos/github/jfnavarro/st_pipeline/badge.svg?branch=master)](https://coveralls.io/github/jfnavarro/st_pipeline?branch=master)
[![PyPI version](https://badge.fury.io/py/stpipeline.svg)](https://badge.fury.io/py/stpipeline)

The ST Pipeline contains the tools and scripts needed to process and analyze the raw files generated with the Spatial Transcriptomics method in FASTQ format to generated datasets for down-stream analysis. 
The ST pipeline can also be used to process single cell data as long as a file with barcodes identifying each cell is provided.

The ST Pipeline has been optimized for speed, robustness and it is very easy to use with many parameters to adjust all the settings.

The following files/parameters are required :
- FASTQ files (Read 1 containing the spatial information and the UMI and read 2 containing the genomic sequence) 
- A genome index generated with STAR 
- An annotation file in GTF or GFF format
- The file containing the barcodes and array coordinates (look at the folder "ids" and chose the correct one). Basically this file contains 3 columns (BARCODE, X and Y), so if you provide this file with barcodes identinfying cells (for example), the ST pipeline can be used for single cell data.
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
- Group annotated reads by barcode(spot position) and gene to get a read count
- In the grouping/counting only unique molecules (UMIs) are kept. 

You can see a graphical more detailed description of the workflow in the documents workflow.pdf and workflow_extended.pdf

The output will be a matrix of counts (genes as columns, spots as rows),
a BED file containing the transcripts (Read name, coordinate, gene, etc..), and a JSON
file with useful stats.
The ST pipeline will also output a log file with useful information.

**Installation**

We recommend you install a virtual environment like Pyenv or Anaconda before you install the pipeline. 
The ST Pipeline works with python 2.7.

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
    
Alternatively, you can install the ST Pipeline using PyPy:

    pip install stpipeline
    
**Example**

An example run would be

	st_pipeline_run.py --ids ids_file.txt --ref-map path_to_index --log-file log_file.txt --output-folder /home/me/results --ref-annotation annotation_file.gtf file1.fastq file2.fastq 

**Emsembl ids**

If you use an Ensembl annotation file and you would like change
the ouput file so it contains gene Ids/names instead of Ensembl ids. 
You can use this tool that comes with the ST Pipeline

	convertEnsemblToNames.py --names-map map.txt --output st_data_updated.tsv st_data.tsv
	
Where map.txt is a tab delimited file with two columns:

ENSEMBL_ID	GENE_NAME

And st_data.tsv is the output from the ST Pipeline.

**Quality stats**

The ST Pipeline generate useful statistical information in the LOG file but if you
want to obtain more detail information about the quality of the data, you can run the following script:

	st_qa.py --input-data stdata.tsv
	
**Documentation**

You can see a more detailed documentation in the folder "doc_out".

**Example data**

You can see a real dataset obtained from the public data from
the following publication (http://science.sciencemag.org/content/353/6294/78)
in the folder called "data".

**License**

The ST pipeline is open source under the MIT license which means that you can use it, change it and re-distribute but you must always refer to our license (see LICENSE and AUTHORS).

**Reference**

If you use the ST Pipeline, please refer to it by including this:

TODO: add reference

**Contact**

For questions, bugs, feedback, etc.. you can contact 
Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

**Dependencies** 

The ST Pipeline depends on some Python packages that will
be automatically installed during the installation process. 
You can see them in the file dependencies.txt

**Requirements**

The ST Pipeline requires to have installed
in the system the aligner STAR (minimum version 2.5.0) :
https://github.com/alexdobin/STAR

The ST Pipeline requieres between
32GB and 64GB of RAM depending
on the size of the data. 
It is recommended to run it
in a system with at least 8 cpu cores. 

