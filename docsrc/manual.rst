Manual
------

ST Pipeline is a tool to process the Spatial Transcriptomics and Visium data.
The data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers.
The output contains the counts matrix and a log file with useful
stats in information.

The ST Pipeline requires two FASTQ files, an IDs files (BARCODE, X, Y),
the path to a STAR genome index, the path to a annotation file in GTF
format and a dataset name.

The ST Pipeline has many parameters, you can see a description of them
by typing : st_pipeline_run.py --help

Note that the minimum read length is dependant on the type of kit used, and
should be adjusted accordingly, i.e. a 150bp kit should have a different
minimum read length than a 75bp kit.

Soft clipping is also not recommended when using the 75bp kit, due to the
shorter length.

The UMI filter can be used for array batches 1000L6 and earlier. It is
not recommended to use it for array batches 1000L7 and newer as the UMI in
these arrays is fully randomised.

Author Jose Fernandez Navarro <jc.fernandez.navaro@gmail.com>


``st_pipeline_run.py [options] fastq_files fastq_files``

**positional arguments**

.. code-block:: bash

  fastq_file_fw
    Read_1 containing the spatial barcodes and UMIs for each sequence.

  fastq_file_rv
    Read_2 containing the gene sequence corresponding to the sequence in
    Read_1.

**optional arguments**

.. code-block:: bash

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
  --htseq-mode [STRING]
                        Mode of annotation when using htseq-count. Modes =
                        {union, intersection-nonempty(default), intersection-
                        strict}
  --htseq-no-ambiguous  When using htseq-count discard reads annotating
                        ambiguous genes (default False)
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
