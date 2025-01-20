# Changelog

## Version 2.0.0
* Refactor code to modern Python (black, mypy, isort)
* Added Github Actions
* Added pre-commit hooks
* Changed the build configuration to Poetry
* Added a Docker container
* Changed tests to pytest
* Updated versions of dependencies
* Performed code optimizations and clean ups
* Added more tests for almost full coveragge
* Bumped taggd to 0.4.0
* Changed documentation
* Added --demultiplex-chunk-size option

## Version 1.8.2
* Added annotation (htseq) feature type as parameter

## Version 1.8.1
* Fix a small bug when having the UMI before the BARCODE

## Version 1.8.0
* Improved the unit-tests
* Make input data a positional arguments in extra scripts

## Version 1.7.9
* --ref-map is not a required parameter now when --disable_mapping is used
* Create dataset and compute saturation steps are only performed if the annotated
file is present
* Fixed bug in merge_fastq.py
* Fixed bugs in st_qa.py

## Version 1.7.8
* Updated docs
* Added flags to skip (trimming, mapping and annotation)
* Fixed a bug in the compute saturation option
* Restored the unittests
* Improvements when dealing with ambiguous and no_feature genes

## Version 1.7.6
* Ported to Python 3

## Version 1.7.3
* Fixed a small bug when soring the output of taggd
* Improved the st_qa.py script

## Version 1.7.2
* Removed parallel code for un-necessary parts
* Added option to provide saturation points

## Version 1.6.5
* Improved the st_qa.py script
* Few small improvements in the annotation step
* Added support to use a transcriptome (--transcriptome)

## Version 1.6.2
* Fixed a error in the parameters

## Version 1.6.1
* Added option to disable the barcode demultiplexing step
* Added option to disable the UMI filtering step

## Version 1.6.0
* Made the parsing of unique UMIs gene by gene and parallel
* AdjacentBi is now the default method for UMI counting
* Made the trimming function output the trimmed R2 in BAM format with the barcode and UMI
* Made the mapping function works with a BAM file as input (latest STAR release)
* Made the annotation function parallel
* Made the quality step parallel
* Improvements in speed and memory (constant memory use)
* The STAR genome loading strategy can be now set
* Added an affinity based method to cluster UMIs
* Added option to set the STAR BAM sort memory limit

## Version 1.5.1
* Fixed a bug that will output the matrix of counts inverted

## Version 1.5.0
* st_qa.py generate expression heatmap plots
* Fixed a minor bug in the computation of the saturation curves
* adjust_matrix_coordinates now does not update the coordinates by default
* adjust_matrix_coordiantes works with the latest ST Spot detector format
* small updates in the fastq merging script
* relaxed a bit the restriction checks for some parameters

## Version 1.4.5
* Added extra scripts:
    - merge_bcl : merge BCL files based in patterns
    - filter_gene_type_matrix : filter gene in output data based on Ensembl gene types
* Bumped Pysam and HTSeq

## Version 1.4.1
* Small update to make the PIP installation more robust

## Version 1.4.0
* Small update to make the PIP installation more robust

## Version 1.3.5
* Optimized the counting of UMIs by strand, start-pos and offset

## Version 1.3.4
* Fixed a typo in one of the parameter that caused the pipeline to not run

## Version 1.3.3
* Disabled spliced alignments by default

## Version 1.3.2
* Optimized the mapping and annotation steps

## Version 1.3.1
* Homopolymer miss-matches is a parameter now
* Removing now also PolyNs (parameter)

## Version 1.3.0
* Added more methods to cluster UMIs
* Optimized the UMI counting algorithm
* Optimized the memory use

## Version 1.2.6
* Take into account soft-clipped bases when computing start/end positions

## Version 1.2.5
* Changed the limit range of some parameters

## Version 1.2.4
* Fixed small bugs
* Small improvements in st_qa.py and convertEnsemblToNames.py

## Version 1.2.3
* Bumped TaggD version
* Added more stats to the dataset output
* Added scripts to compute stats
* Added new option for TaggD

## Version 1.2.2
* Fixed bugs in convertEnsemblToNames
* Added some parameters for TaggD demultiplexing
* Bumped version of TaggD

## Version 1.2.1
* Made homopolymers filters enabled by default
* Added a test dataset to the docs

## Version 1.2.0
* Fixed a small bug in the deletion of the tmp folder

## Version 1.1.7
* Make sure to remove tmp files even if an error happens

## Version 1.1.6
* Fixed bug that would leave some files in /tmp
* Allowed mis-matches when removing adaptors is now 2

## Version 1.1.5
* Removed some un-necessary parameters

## Version 1.1.1
* Simplified the two pass mode

## Version 1.1.0
* Added flag to discard reads mapping to anti-sense strand
* Parameters for GC content filter instead of using the same value as AT content filter
* Fixed a small bug in the logging of some parameters

## Version 1.0.4
* When removing adaptors (homopolymers streches) allow to up to 3 missmatches
* Added GC content filter (same % as AT content)

## Version 1.0.3
* Fixed a minor bug in the counting of UMIs or - strand

## Version 1.0.2
* If no temp folder is given a new unique one is created on top of the execution folder
* integrate createDataset.py into the code of the pipeline
* Adjusted some parameters names and descriptions (no UMI is default)
* Added sliding window when counting unique molecules
* Added support to bzip

## Version 1.0.1
* Fixed small bug in the parsing of the umi quality parameter

## Version 1.0.0
* Added option to check for UMI quality
* Optimized the UMI template check code
* Optimized how the unique molecules are counted
* Better stats for the quality filter step
* Updated convertEnsemblToNames script
* Updated stringdocs

## Version 0.9.9
* Small bug fixes

## Version 0.9.6
* Fixed a bug with the non ambiguous option
* Fix a bug in the saturation computation

## Version 0.9.5
* When a R2 is trimmed its correspondant R1 is trimmed as well

## Version 0.9.4
* Fixed a stupid bug in the compute saturation option

## Version 0.9.3
* Changed the rRNA filter so the BAM output does not need to be sorted

## Version 0.9.2
* Fixed a bug in the parsing of parameters

## Version 0.9.1
* Fixed a small bug with the location of discarded files

## Version 0.9.0
* Replaced JSON for data frame in the output format
* Replaced python gzip for system call (faster)
* Changed the logic of how the filenames are stored and handled

## Version 0.8.9
* Improved the error messages and error handling

## Version 0.8.8
* Removed barcodes IDs from the output file

## Version 0.8.7
* Updated comments, manual and license
* Small improvements

## Version 0.8.5
* Fixed a bug in the computation of saturation curves

## Version 0.8.4
* Added a normal hash with INT keys to increase speed and reduce memory
* Using the gene_id for annotation again

## Version 0.8.3
* Added parameter for strandness in annotation (yes by default)
* Simplified a bit the quality trimming step (do not account for user input trimmed bases)

## Version 0.8.2
* Added stats for annotated reads
* Replaced shelve dict for sqldict
* Fixed some small bugs in the annotation

## Version 0.8.1
* Removed the pair mode keep option
* Removed un-neccessary pair mode and mapped checks
  after alignment

## Version 0.8.0
* Added option to do the STAR 2 pass mode
* Removed option to run pipeline without IDs
* Speed improvements
* Perform demultiplex after mapping
* No attaching the barcode to reverse reads
* Removing some parameters
* Some improvements in stDataPlotter
* Option to use BAM format
* Removed annotation filtering step
* Removed forward trimming parameters
* Output gene names even with ENSEMBL

## Version 0.7.7
* Small memory improvements
* Updates in plotting script

## Version 0.7.6
* End coordinates now contain the whole read length
* Make annotation strand aware (reverse)
* Updated to STAR 2.5

## Version 0.7.5
* Fixed a small bug

## Version 0.7.4
* Added some memory improvements

## Version 0.7.3
* Added parameters for inverse trimming
* Memory and speed optimizations in createDatasets
* Added option for low_memory use

## Version 0.7.2
* Added unique genes to saturation points
* Added option to keep non-annotated reads

## Version 0.7.1
* Fixed some small bugs

## Version 0.7.0
* Fixed a bug in the saturation points
* Removed counttrie as option for clustering
* Updated and improved CTTS scripts
* Updated datfa plotter color list

## Version 0.6.9
* Fixed a bug in the saturation points

## Version 0.6.8
* Improved speed and memory in createDatasets
* Changed saturation points to fixed values that grow exp
* Improved speed in computation of saturation points
* Small bug fixes
* Upgraded json2Scatter with many improvements
* Rename json2scatter to stDataPlotter

## Version 0.6.7
* Fixed a bug in the hierarchical clustering
* Added the input parameter to qa_stats
* Append experiment name to output files
* Added option to compute saturation points
* Added tool to plot stdata and clusters with aligned image

## Version 0.6.6
* Fixed a bug in the hierarchical clustering
* Fixed a bug in the printed stats

## Version 0.6.5
* Fixed a bug in retrieving the version of the software
* Added time stamps in different steps
* Added a UMI template quality filter

## Version 0.6.4
* Fixed a bug in counttrie clustering method
* Improved sorting of molecular barcodes prior clustering
* Added hiearachical clustering option

## Version 0.6.3
* Removed reads.json
* Added qa_stats.json to the output
* Restored old versioning system
* Removed hadoop related stuff
* Added support for gziped input files

## Version 0.6.2
* Improved the log a bit
* Added parameters for max,min intron size and max gap size

## Version 0.6.1
* Fixed some bugs in the prefix trie

## Version 0.5.9
* Added an option to find molecular barcodes clusters using a prefix trie

## Version 0.5.8
* Fixed a bug in the function to retrieve the pipeline version

## Version 0.5.7
* Fixed a bug with --disable-multimap option

## Version 0.5.6
* Fixed a typo in a parameter
* Fixed a bug that caused some parameters to not work

## Version 0.5.5
* Added some extra debugging info in createDatasets
* Output the read name in the BED output file
* Changed --allowed-kimera for --allowed-kmer
* Added version as parameter and log message

## Version 0.5.4
* Added parameter to disable soft clipping in mapping
* Disable softclipping in rRNA filter
* Make sure that discarded reads after rRNA filter are replaced by Ns
* Improved stats info a bit

## Version 0.5.3
* Bumped Taggd to 0.2.2

## Version 0.5.2
* Fixed a bug in the rRNA filter that would cause to not discard
rRNA mapped reads

## Version 0.5.1
* Added check when UMI is the same as barcode
* Added more stats
* Added percentiles distributiosn stats for createDAtaset
* Added support for BAM and SAM (not functional now)
* Added option to disable multiple aligned reads
* Fixed a bug in the bed file

## Version 0.5.0
* Added AT content filter in quality trimming
* Added min mapped length filter after mapping
* Make sure one of the multiple aligned reads is set as not multiple
aligned so it can be annotated
* Discard the other multiple aligned reads after mapping
* Disable sorting
* Restored back to use gene_id as column for annotation

## Version 0.4.9
* Changed naming convention
* Added support for normal RNA analysis

## Version 0.4.8
* Improved STAR configuration
* Added mapping post processing to filter out and adjust reversed reads
* Changed to use gene_name for annotation
* Fixed some bugs and some improvements
* Fixed bugs in the trimming

## Version 0.4.7
* Improved stats
* Fixed a bug that would remove original input files
* Added a script to convert ENSEMBL ids to gene names

## Version 0.4.6
* Fixed a bug that would not compute the number of discarded reads
when using molecular barcodes

## Version 0.4.5
* Fixed a bug in the barcodes JSON output

## Version 0.4.4
* Fixed a bug in the molecular barcodes algorithm
* Fixed a bug that would keep the original fastq reads in the system
* Update taggd version

## Version 0.4.3
* Small improvements with error checking and log in the mapping
* Fixed a bug that would remove the file after filtering annoted reads
* Make the sorting by name instead by position due to a bug in htseq-count

## Version 0.4.2
* Fixed a bug in the capture of parameters

## Version 0.4.1
* Improved the logs
* Fixed few bugs

## Version 0.4.0
* Added back taggd
* Added BED file to output
* Added STAR
* Optimized workflow
* do rRNA filter first
* Optimized annotation
* Optimized trimming
* Output reads do not contain duplicates

## Version 0.3.9
* Allowing molecular barcodes to be before the barcodes

## Version 0.3.8
* Added back findIndexes

## Version 0.3.7
* Removed cutadapt dependency

## Version 0.3.6
* Fixed a bug in the installation

## Version 0.3.5
* Added options to remove PolyC fix bugs in adaptors removal

## Version 0.3.4
* Added test for STAR and STAR binary to dependencies
* Added TAGGD and removed findIndexes
* Improved install script
* Added options to remove adaptors (PolyA, PolyT and PolyG)
* Exchanged Bowtie as primary mapper with STAR.

## Version 0.3.3
* Added option to keep files with discarded reads/barcodes
* Internal refactoring and optimization

## Version 0.3.2
* Outputted reads JSON now only has the portion of the read that was used to map
* Cutadapt is integrated but only using the quality trimming for now
* Internal refactoring and optimizations

## Version 0.3.1
* Added small unit-test for molecular barcodes
* Added more molecular barcodes algorithms (using a naive one for now)
* Fixed small issues in JSON parsing libraries

## Version 0.3.0
* Rewrite createDatasets.py
* Clean up repository and deprecated files
* Change the unit-test library and structure
* Refactor the unit-test (use pipeline API instead of command line calls)
* Ensure unit-test remove tmp files when failing
* Add better error handling
* Add unit-test for Molecular Barcodes
* Add Molecular Barcodes functionality
* General refactor and clean up
* Add invoke options (clean, build, install)
* Fix an important bug in createDatasets that caused incorrect
  computation of reads counts

## Version 0.2.5
* Improved installers
* Small bug fixes
* Added basic uni-test to do a run of the pipeline

## Version 0.2.4
* Some optimizations and bug fixes

## Version 0.2.3
* Fixed a error with new version of HTSeq-count that will discard more reads

## Version 0.2.2
* Added extra parameters
* Fixed some typos
* Fixed a bug that caused to remove some bases from the barcode ID in the rw reads

## Version 0.2.1
* code refactored and modularized
* add argparse for parameters parsing
* add API for Amazon EMR and terminal version
* better error handling
* optimized code
* new version of FindIndexes
* remove dependencies
* added proper installers and documentation
