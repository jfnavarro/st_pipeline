Manual
------

ST Pipeline is a tool to process the Spatial Transcriptomics raw data
or single cell data.
The data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers.
The output contains the counts matrix, a stats file, a log file
and a BED file with all the transcripts.

The ST Pipeline requires two fastq files, an IDs files (BARCODE, X, Y),
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
unable to be used for array batches 1000L7 and newer as the UMI in
these arrays is fully randomised.

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

.. code-block:: bash

  usage: st_pipeline_run.py [-h]
                          --ids [FILE]
                          --ref-map [FOLDER]
                          --ref-annotation [FILE]
                          --expName [STRING]
                          [--allowed-missed [INT]]
                          [--allowed-kmer [INT]]
                          [--overhang [INT]]
                          [--min-length-qual-trimming [INT]]
                          [--mapping-rv-trimming [INT]]
                          [--length-id [INT]]
                          [--contaminant-index [FOLDER]]
                          [--qual-64]
                          [--htseq-mode [STRING]]
                          [--htseq-no-ambiguous]
                          [--start-id [INT]]
                          [--no-clean-up]
                          [--verbose]
                          [--mapping-threads [INT]]
                          [--min-quality-trimming [INT]]
                          [--bin-path [FOLDER]]
                          [--log-file [STR]]
                          [--output-folder [FOLDER]]
                          [--temp-folder [FOLDER]]
                          [--umi-allowed-mismatches [INT]]
                          [--umi-start-position [INT]]
                          [--umi-end-position [INT]]
                          [--keep-discarded-files]
                          [--remove-polyA [INT]]
                          [--remove-polyT [INT]]
                          [--remove-polyG [INT]]
                          [--remove-polyC [INT]]
                          [--remove-polyN [INT]]
                          [--filter-AT-content [INT%]]
                          [--filter-GC-content [INT%]]
                          [--disable-multimap]
                          [--disable-clipping]
                          [--umi-cluster-algorithm [STRING]]
                          [--min-intron-size [INT]]
                          [--max-intron-size [INT]]
                          [--umi-filter]
                          [--umi-filter-template [STRING]]
                          [--compute-saturation]
                          [--include-non-annotated]
                          [--inverse-mapping-rv-trimming [INT]]
                          [--low-memory]
                          [--two-pass-mode]
                          [--strandness [STRING]]
                          [--umi-quality-bases [INT]]
                          [--umi-counting-offset [INT]]
                          [--demultiplexing-metric [STRING]]
                          [--demultiplexing-multiple-hits-keep-one]
                          [--demultiplexing-trim-sequences [INT]]
                          [--version]
                          fastq_file_fw fastq_file_rv




**positional arguments**

.. code-block:: bash

  fastq_file_fw
    Read_1 containing the spatial barcodes and UMIs for each sequence.

  fastq_file_rv
    Read_2 containing the gene sequence corresponding to the sequence in
    Read_1.

**optional arguments**

.. code-block:: bash

  -h, --help                          Show this help message and exit.
  --ids [FILE]                        Path to the file containing the map of
                                      barcodes to the array coordinates.
  --ref-map [FOLDER]                  Path to the folder with the STAR index
                                      for the genome that you want to use to
                                      align the reads.
  --ref-annotation [FILE]             Path to the reference annotation file
                                      (GTF or GFF format is required) to be
                                      used to annotated the reads.
  --expName [STRING]                  Name of the experiment/dataset
                                      (The output files will prepend this
                                      name).
  --allowed-missed [INT]              Number of allowed mismatches when
                                      demultiplexing against the barcodes
                                      with TaggD (default: 2).
  --allowed-kmer [INT]                KMer length when demultiplexing against
                                      the barcodes with TaggD (default: 6).
  --overhang [INT]                    Extra flanking bases added when
                                      demultiplexing against the barcodes.
  --min-length-qual-trimming [INT]    Minimum length of the reads after
                                      trimming, shorter reads will be
                                      discarded (default: 25).
  --mapping-rv-trimming [INT]         Number of bases to trim in the reverse
                                      reads for the mapping step (5' end)
                                      (default: 0).
  --length-id [INT]                   Length of IDs
                                      (the length of the barcodes)
                                      (default: 18).
  --contaminant-index [FOLDER]        Path to the folder with a STAR index
                                      with a contaminant genome. Reads will
                                      be filtered against the specified
                                      genome and mapping reads will be
                                      discarded.
  --qual-64                           Use phred-64 quality instead of
                                      phred-33(default).
  --htseq-mode [STRING]               Mode of Annotation when using HTSeq.
                                      Modes = {union ,
                                      intersection-nonempty(default),
                                      intersection-strict}.
  --htseq-no-ambiguous                When using htseq discard reads
                                      annotating ambiguous genes
                                      (default False).
  --start-id [INT]                    Start position of the IDs (Barcodes)
                                      in the R1 (counting from 0)
                                      (default: 0).
  --no-clean-up                       Do not remove temporary/intermediary
                                      files (useful for debugging).
  --verbose                           Show extra information on the log file.
  --mapping-threads [INT]             Number of threads to use in the mapping
                                      step (default: 4).
  --min-quality-trimming [INT]        Minimum phred quality a base must have
                                      in the trimming step (default: 20).
  --bin-path [FOLDER]                 Path to folder where binary executables
                                      are present (system path by default).
  --log-file [STR]                    Name of the file that we want to use to
                                      store the logs
                                      (default output to screen).
  --output-folder [FOLDER]            Path of the output folder.
  --temp-folder [FOLDER]              Path of the location for temporary
                                      files.
  --umi-allowed-mismatches [INT]      Number of allowed mismatches
                                      (hamming distance) that UMIs of the
                                      same gene-spot must have in order to
                                      cluster together (default: 1).
  --umi-start-position [INT]          Position in R1 (base wise) of the first
                                      base of the UMI (starting by 0)
                                      (default: 18).
  --umi-end-position [INT]            Position in R1 (base wise) of the last
                                      base of the UMI (starting by 1)
                                      (default: 27).
  --keep-discarded-files              Writes down unaligned, un-annotated
                                      and un-demultiplexed reads to files.
  --remove-polyA [INT]                Remove PolyA stretches of the given
                                      length from R2 (default: 15).
  --remove-polyT [INT]                Remove PolyT stretches of the given
                                      length from R2 (default: 15).
  --remove-polyG [INT]                Remove PolyG stretches of the given
                                      length from R2 (default: 15).
  --remove-polyC [INT]                Remove PolyC stretches of the given
                                      length from R2 (default: 15).
  --remove-polyN [INT]                Remove PolyN stretches of the given
                                      length from R2 (default: 15).
  --filter-AT-content [INT%]          Discards reads whose number of A and T
                                      bases in total are more or equal than
                                      the number given in percentage
                                      (default: 90).
  --filter-GC-content [INT%]          Discards reads whose number of G and C
                                      bases in total are more or equal than
                                      the number given in percentage
                                      (default: 90).
  --disable-multimap                  If activated, multiple aligned reads
                                      obtained during mapping will be all
                                      discarded. Otherwise the highest scored
                                      one will be kept.
  --disable-clipping                  If activated, disable soft-clipping
                                      (local alignment) in the mapping step.
  --umi-cluster-algorithm [STRING]    Type of clustering algorithm to use
                                      when performing UMIs duplicates
                                      removal.
                                      Modes = {naive(default), hierarchical, Adjacent and AdjacentBi}.
  --min-intron-size [INT]             Minimum allowed intron size when
                                      searching for splice variants in the
                                      mapping step (default: 20).
  --max-intron-size [INT]             Maximum allowed intron size when
                                      searching for splice variants in the
                                      mapping step (default: 100000).
  --umi-filter                        Enables the UMI quality filter based on
                                      the template given in
                                      --umi-filter-template.
  --umi-filter-template [STRING]      UMI template (IUPAC nucleotide code)
                                      for the UMI filter, default = WSNNWSNNV
  --compute-saturation                Performs a saturation curve computation
                                      by sub-sampling the annotated reads,
                                      computing unique molecules and then a
                                      saturation curve
                                      (included in the log file).
  --include-non-annotated             Do not discard un-annotated reads
                                      (they will be labeled __no_feature)
  --inverse-mapping-rv-trimming [INT] Number of bases to trim in the reverse
                                      reads for the mapping step on the
                                      3' end.
  --low-memory                        Writes temporary records into disk in
                                      order to save memory but gaining a
                                      speed penalty.
  --two-pass-mode                     Activates the 2 pass mode in STAR to
                                      also map against splice variants.
  --strandness [STRING]               What strandness mode to use when
                                      annotating with htseq-count
                                      [no, yes(default), reverse].
  --umi-quality-bases [INT]           Maximum number of low quality bases
                                      allowed in an UMI (default: 8).
  --umi-counting-offset [INT]         Expression count for each gene-spot
                                      combination is expressed as the number
                                      of unique UMIs in each strand/start
                                      position. However some reads might have
                                      slightly different start positions due
                                      to amplification artifacts. This
                                      parameters allows one to define an
                                      offset from where to count unique UMIs
                                      (default: 150).
  --demultiplexing-metric             Distance metric for TaggD demultiplexing: 
                                      Subglobal, Levenshtein or Hamming 
                                      (default: Subglobal)
  --demultiplexing-multiple-hits-keep-one  When multiple ambiguous hits with same score are 
                                      found in the demultiplexing, keep one (random)
  --demultiplexing-trim-sequences     Trims from the barcodes in the input file when doing demultiplexing.
                            	      The bases given in the list of tuples as START END START END .. where
                                      START is the integer position of the first base (0 based) and END is the integer
                                      position of the last base (1 based).
                                      Trimmng sequences can be given several times.
  --version                           Show program's version number and exit
