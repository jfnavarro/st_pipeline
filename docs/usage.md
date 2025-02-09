# Usage

## Input format

The following files/parameters are commonly required:

- FASTQ files (Read 1 containing the barcode and the UMI and read 2 containing the genomic sequence)
- A genome index generated with STAR
- An annotation file in GTF or GFF3 format (optional when using a transcriptome)
- A file containing the barcodes and array coordinates (look at the folder "ids" to use it as a reference).
Basically this file contains 3 columns (BARCODE, X and Y), so if you provide this
file with barcodes identifying cells (for example), the ST pipeline can be used for single cell data.
This file is also optional if the data is not barcoded (for example RNA-Seq data).
- A name for the dataset

The ST pipeline has multiple parameters mostly related to trimming, mapping, demultiplexing and annotation
but generally the default values are good enough.

The input FASTQ files can be given in gzip/bzip format as well.

Note that the minimum read length is dependant on the type of kit used, and
should be adjusted accordingly, i.e. a 150bp kit should have a different
minimum read length than a 75bp kit.

Soft clipping is also not recommended when using the 75bp kit, due to the shorter length.

The UMI filter can be used for array batches 1000L6 and earlier. It is
not recommended to use it for array batches 1000L7 and newer as the UMI in
these arrays is fully randomised.

The output of the ST Pipeline is a counts matrix (TSV) and a log file with status and useful information.

You can see a full description of the parameters typing `st_pipeline_run --help` after you have installed the ST pipeline.

```console
st_pipeline_run [options] fastq_file_fw fastq_file_rv

  fastq_file_fw
    Read_1 containing the spatial barcodes and UMIs for each sequence.

  fastq_file_rv
    Read_2 containing the gene sequence corresponding to the sequence in
    Read_1.

  -h, --help            show this help message and exit
  --ids [FILE]          Path to the file containing the map of barcodes to the
                        array coordinates
  --ref-map [FOLDER]    Path to the folder with the STAR index for the genome
                        that you want to use as reference
  --ref-annotation [FILE]
                        Path to the reference annotation file (GTF or GFF
                        format is required) to be used to annotated the mapped
                        reads
  --expName [STRING]    Name of the dataset (The output files will prepend
                        this name)
  --contaminant-index [FOLDER]
                        Path to the folder with a STAR index with a
                        contaminant genome reference. Reads will be filtered
                        using the specified genome and mapping reads will be
                        discarded
  --no-clean-up         Do not remove temporary/intermediary files (useful for
                        debugging)
  --verbose             Show extra information on the log file
  --threads [INT]       Number of threads to use (default: 4)
  --bin-path [FOLDER]   Path to folder where binary executables are present
                        (system path by default)
  --log-file [STR]      Name of the file that we want to use to store the logs
                        (default output to screen)
  --output-folder [FOLDER]
                        Path of the output folder
  --temp-folder [FOLDER]
                        Path of the location for temporary files
  --keep-discarded-files
                        Keep files with discarded reads in every step
  --qual-64             Use phred-64 quality instead of phred-33(default) in
                        the quality trimming step
  --min-length-qual-trimming [INT]
                        Minimum length of the reads after trimming, shorter
                        reads will be discarded (default: 20)
  --min-quality-trimming [INT]
                        Minimum phred quality a base must have in order to be
                        kept in the quality trimming step (default: 20)
  --remove-polyA [INT]  Remove PolyA stretches of the given length from R2
                        (Use 0 to disable it) (default: 10)
  --remove-polyT [INT]  Remove PolyT stretches of the given length from R2
                        (Use 0 to disable it) (default: 10)
  --remove-polyG [INT]  Remove PolyG stretches of the given length from R2
                        (Use 0 to disable it) (default: 10)
  --remove-polyC [INT]  Remove PolyC stretches of the given length from R2
                        (Use 0 to disable it) (default: 10)
  --remove-polyN [INT]  Remove PolyN stretches of the given length from R2
                        (Use 0 to disable it) (default: 10)
  --homopolymer-mismatches [INT]
                        Number of mismatches allowed when removing
                        homopolymers (A, T, G, C or N) (default: 0)
  --filter-AT-content [INT%]
                        Discards reads whose number of A and T bases in total
                        are more or equal than the percentage given as input
                        (0-100) (default: 90)
  --filter-GC-content [INT%]
                        Discards reads whose number of G and C bases in total
                        are more or equal the percentage given as input
                        (0-100) (default: 90)
  --mapping-rv-trimming [INT]
                        Number of bases to trim in the reverse reads (R2) for
                        the mapping step (5' end) (default: 0)
  --inverse-mapping-rv-trimming [INT]
                        Number of bases to trim in the reverse reads (R2) for
                        the mapping step (3' end) (default: 0)
  --disable-multimap    If activated, multiple aligned reads obtained during
                        mapping will be all discarded. Otherwise the highest
                        scored one will be kept
  --disable-clipping    If activated, disable soft-clipping (local alignment)
                        in the mapping step
  --min-intron-size [INT]
                        Minimum allowed intron size when searching for splice
                        variants with STAR Splices alignments are disabled by
                        default (=1) but to turn it on set this parameter to a
                        bigger number, for example 10 or 20. (default: 1)
  --max-intron-size [INT]
                        Maximum allowed intron size when searching for splice
                        variants with STAR Splices alignments are disabled by
                        default (=1) but to turn it on set this parameter to a
                        big number, for example 10000 or 100000. (default: 1)
  --star-two-pass-mode  Activates the 2-pass mode in STAR to improve mapping
                        accuracy
  --star-genome-loading [STRING]
                        Similar to the STAR option --genomeLoad. It allows to
                        load the genome index into memory so it can easily be
                        shared by other jobs to save loading time. Read the
                        STAR manual for more info on this. (default:
                        NoSharedMemory)
  --star-sort-mem-limit STAR_SORT_MEM_LIMIT
                        The maximum available RAM for sorting BAM during
                        mapping with STAR. Default is 0 which means that it
                        will be set to the genome index size
  --demultiplexing-mismatches [INT]
                        Number of allowed mismatches when demultiplexing the
                        reads against the barcodes with TaggD (default: 2)
  --demultiplexing-kmer [INT]
                        KMer size to use when demultiplexing against the
                        barcodes with TaggD (default: 6)
  --demultiplexing-overhang [INT]
                        Extra flanking bases added on each side of the barcode
                        when demultiplexing against the barcodes with TaggD
                        (default: 0)
  --demultiplexing-start [INT]
                        Start position of the IDs (Barcodes) in R1 (counting
                        from 0) (default: 0)
  --demultiplexing-metric [STRING]
                        Distance metric to use for TaggD demultiplexing:
                        Options: Subglobal, Levenshtein or Hamming (default:
                        Subglobal)
  --demultiplexing-multiple-hits-keep-one
                        When multiple ambiguous hits with same score are found
                        in the demultiplexing step, keep only one (random).
  --demultiplexing-trim-sequences DEMULTIPLEXING_TRIM_SEQUENCES
                        Trim the barcodes in the input file when doing
                        demultiplexing. The input given is a list of tuples
                        START END START END where START is the integer
                        position of the first base (0 based) and END is the
                        integer position of the last base (1 based). The final
                        barcode will be obtained by combining all the
                        sequences given in the input. This is useful when
                        having a barcode composed of multiple sequences in the
                        reador when the barcode needs to be trimmed out.
                        Trimmng sequences can be given several times.
  --demultiplexing-chunk-size [INT]
                        Chunk size for parallel processing (number of reads assigned to each thread) (default: 10000)
  --htseq-mode [STRING]
                        Mode of annotation when using htseq-count. Modes =
                        {union, intersection-nonempty(default), intersection-
                        strict}
  --htseq-no-ambiguous  When using htseq-count discard reads annotating
                        ambiguous genes (default False)
  --htseq-features HTSEQ_FEATURES [HTSEQ_FEATURES ...]
                        Which feature types to use from the GTF/GFF file in the annotation.
                        Can be given more than one type (default exon)
  --strandness [STRING]
                        What strandness mode to use when annotating with
                        htseq-count [no, yes(default), reverse]
  --include-non-annotated
                        Do not discard un-annotated reads (they will be
                        labeled __no_feature)
  --umi-cluster-algorithm [STRING]
                        Type of clustering algorithm to use when performing
                        UMIs duplicates removal. Options = {naive,
                        hierarchical, Affinity, Adjacent and
                        AdjacentBi(default)} Note that for the affinity method
                        the umi allowed mismatches parameter will be ignored.
  --umi-allowed-mismatches [INT]
                        Number of allowed mismatches (hamming distance) that
                        UMIs of the same gene-spot must have in order to
                        cluster together (default: 1)
  --umi-start-position [INT]
                        Position in R1 (base wise) of the first base of the
                        UMI (starting by 0) (default: 18)
  --umi-end-position [INT]
                        Position in R1 (base wise) of the last base of the UMI
                        (starting by 1) (default: 27)
  --umi-filter          Enables the UMI quality filter based on the template
                        given in --umi-filter-template
  --umi-filter-template [STRING]
                        UMI template (IUPAC nucleotide code) for the UMI
                        filter, default = WSNNWSNNV
  --umi-quality-bases [INT]
                        Maximum number of low quality bases allowed in an UMI
                        (default: 6)
  --umi-counting-offset [INT]
                        UMI count for each gene-spot combination is computed
                        as the number of unique UMIs in each strand/start
                        position. However some reads might have slightly
                        different start positions due to amplification
                        artifacts. This parameters allows to define an offset
                        window from where to count unique UMIs. You can set it
                        to a very high value +9999 to count unique UMIs for
                        the whole gene (default: 250)
  --compute-saturation  Performs a saturation curve computation by sub-
                        sampling the annotated reads, computing unique UMIs
                        and adding the stats to the log file (this can be used
                        to plot saturation curves)
  --saturation-points SATURATION_POINTS [SATURATION_POINTS ...]
                        Saturation points for the saturation curve computation
                        can be provided instead of using default values.
                        Provide a list of values like this for example: 10000
                        20000 50000 100000
  --disable-trimming    Use this flag if you want to skip the trimming step
  --disable-mapping     Use this flag if you want to skip the mapping step
  --disable-annotation  Use this flag if you want to skip the annotation
  --disable-barcode     Use this flag if you want to skip the barcode
                        demultiplexing step
  --disable-umi         Use this flag if you want to skip the UMI filtering
                        step
  --transcriptome       Use this flag if you want to use transcriptome instead
                        of a genome, the gene tag will be obtained from the
                        transcriptome file
  --version             show program's version number and exit
```

## Example

An example run would be:

```console
st_pipeline_run --expName test --ids ids_file.txt \
  --ref-map path_to_index --htseq-no-ambiguous --log-file log_file.txt --output-folder /home/me/results \
  --ref-annotation annotation_file.gtf --contaminant-index path_to_cont_index file1.fastq file2.fastq
```

## Visium

To process Visium datasets it is recommended to use these options:

```console
--demultiplexing-mismatches 1
--demultiplexing-kmer 4
--umi-allowed-mismatches 2
--umi-start-position 16
--umi-end-position 28
```

## Extra tools

The ST Pipeline ships many scripts that will be installed automatically and that
can be very useful to pre/post process the data and general QC stats and plots.

### Emsembl ids

If you used an Ensembl annotation file and you would like change
the output file so it contains gene ids/names instead of Ensembl ids.
You can use the script `convertEnsemblToNames` that comes with the ST Pipeline

```console
convertEnsemblToNames --annotation path_to_annotation_file --output st_data_updated.tsv st_data.tsv
```

### Merge demultiplexed FASTQ files

If you used different indexes to sequence and need to merge the fastq files
you can use the script `merge_fastq` that comes with the ST Pipeline

```console
merge_fastq --run-path path_to_run_folder --out-path path_to_output --identifiers S1 S2 S3 S4
```

Where `--identifiers` will be strings that identify each demultiplexed sample.

## Filter out genes by gene type

If you want to remove from the dataset (matrix in TSV) genes corresponding
to certain gene types (For instance to keep only protein_coding). You can do
so with the script `filter_gene_type_matrix` that comes with the ST Pipeline

```console
filter_gene_type_matrix --gene-types-keep protein-coding --annotation path_to_annotation_file stdata.tsv
```

You may include the parameter `--ensembl-ids` if your genes are represented as Emsembl ids instead.

The value of `--gene-types-keep` must match the annotation file provided.

### Remove spots from dataset

If you want to remove spots from a dataset (matrix in TSV) for instance
to keep only spots inside the tissue. You can do so with the script `adjust_matrix_coordinates`
that comes with the ST Pipeline

```console
adjust_matrix_coordinates --outfile new_stdata.tsv --coordinates-file coordinates.txt stdata.tsv
```

Where `coordinates.txt` must be a tab delimited file with 6 columns:

```console
orig_x orig_y new_x new_y new_pixel_x new_pixel_y
```

Only spots whose coordinates in the file will be kept and then optionally you
can update the coordinates in the matrix choosing for the new array or pixel coordinates.

### Quality stats

The ST Pipeline generate useful stats/QC information in the LOG file but if you
want to obtain more detailed information about the quality of the data, you can run the following script:

```console
st_qa stdata.tsv
```

If you want to perform quality stats on multiple datasets you can run:

```console
multi_qa stdata1.tsv stadata2.tsv stdata3.tsv stdata4.tsv
```

`multi_qa` generates violing plots, correlation plots/tables and more useful information and
it allows to log the counts for the correlation analaysis.
