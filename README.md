# Spatial Transcriptomics (ST) Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-310/)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-311/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-312/)
[![PyPI version](https://badge.fury.io/py/stpipeline.svg)](https://badge.fury.io/py/stpipeline)
[![Build Status](https://github.com/jfnavarro/st_pipeline/actions/workflows/dev.yml/badge.svg)](https://github.com/jfnavarro/st_pipeline/actions/workflows/dev)

The ST Pipeline provides the tools, algorithms and scripts needed to process and analyze the raw
data generated with Spatial Transcriptomics or Visium in FASTQ format to generate datasets for down-stream analysis.

The ST Pipeline can also be used to process single cell/nuclei RNA-seq data as long as a
file with molecular `barcodes` identifying each cell is provided (same template as the files in the folder "ids").

The ST Pipeline can also be used to process bulk RNA-seq data, in this case the barcodes file is not required.

The ST Pipeline has been optimized for speed, robustness and it is very easy to use with many parameters to adjust all the settings.
The ST Pipeline is fully parallel and it has constant memory use.
The ST Pipeline allows to skip any of the main steps and provides multiple customization options.
The ST Pipeline allows to use either the genome or the transcriptome as reference.

Basically what the ST pipeline does (default mode) is:

- Quality trimming step (read 1 and read 2):
  - Remove low quality bases
  - Sanity check (reads same length, reads order, etc..)
  - Check quality UMI
  - Remove artifacts (PolyT, PolyA, PolyG, PolyN and PolyC) of user defined length
  - Check for AT and GC content
  - Discard reads with a minimum number of bases of that failed any of the checks above
- Contamimant filter step (e.x. rRNA genome) (Optional)
- Mapping with [STAR](https://github.com/alexdobin/STAR) step (only read 2) (Optional)
- Demultiplexing with [Taggd](https://github.com/jfnavarro/taggd) step (only read 1) (Optional)
- Keep reads (read 2) that contain a valid barcode and are correctly mapped
- Annotate the reads to the reference (Optional)
- Group annotated reads by barcode (spot position), gene and genomic location (with an offset) to get a read count
- In the grouping/counting only unique molecules (UMIs) are kept (Optional)

You can see a graphical more detailed description of the workflow in the documents `workflow.pdf` and `workflow_extended.pdf`

The output dataset is a matrix of counts (genes as columns, spots as rows) in TSV format.
The ST pipeline will also output a log file with useful stats and information.

## Installation

For users see [install](docs/installation.md)

For developers [contributing](CONTRIBUTING.md)

## Usage

See [usage](docs/usage.md)

## Authors

See [authors](AUTHORS.md)

## License

The ST pipeline is open source under the MIT license which means that you can use it,
change it and re-distribute but you must always refer to our license (see LICENSE).

## Credits

If you use the ST Pipeline, please refer its publication:
ST Pipeline: An automated pipeline for spatial mapping of unique transcripts
Oxford BioInformatics
10.1093/bioinformatics/btx211

## Example dataset

You can see a real dataset obtained from the public data from
the following publication (http://science.sciencemag.org/content/353/6294/78)
in the folder called "data".

## Contact

For questions, bugs, feedback, etc.. you can contact:

Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
