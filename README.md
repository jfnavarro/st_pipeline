ST Pipeline contains the tools and scripts needed to process and analyze the sequenced
reads files generated in the Spatial Transcriptomics group in FASTQ format. 

The required input files are two fastq files, a genome reference, an annotation reference and the file
containing the IDs (from the array-chip-plate) plus the optional configuration parameters. 

The output will be two JSON files and a BED file, one containing
the features (id,cordinates,gene) and another one containing the raw reads and the id.

We recommend you install a virtual environment like Pyenv or Anaconda. 

The ST Pipeline has the following dependencies :

HTSeq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License, version 3, can be found here: http://www.gnu.org/licenses/gpl-3.0-standalone.html

STAR
Spliced Transcripts Alignment to a Reference
© Alexander Dobin, 2009-2015

STAR is a free open source software distributed under GPLv3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For details, see <http://www.gnu.org/licenses/>.

To install the pipeline type 

- python setup.py
- python setup.py install

To run a test type
- python setup.py test

To see the different options type -h or --help

