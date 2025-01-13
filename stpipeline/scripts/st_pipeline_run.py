#! /usr/bin/env python
"""
ST Pipeline is a tool to process Spatial Transcriptomics raw datasets (FASTQ).
The raw data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers (UMI).
The output contains the counts matrix (TSV), a stats file, a log file
and a BED file with all the transcripts.

The ST Pipeline requires two FASTQ files, an IDs files (BARCODE, X, Y),
the path to a STAR genome index, the path to a annotation file in GTF format
an a dataset name.

The ST Pipeline has many parameters and options, you can see a description of them
by typing : st_pipeline_run.py --help

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys

from stpipeline.core.pipeline import Pipeline


def main() -> int:
    # Create pipeline object
    pipeline = Pipeline()

    # Create a parser
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Parse parameters, sanity check and run the pipeline
    try:
        parser = pipeline.createParameters(parser)

        # Parse arguments
        options = parser.parse_args()
        pipeline.load_parameters(options)
        print("ST Pipeline, parameters loaded")

        # Create logger
        pipeline.createLogger()
        print("ST Pipeline, logger created")

        # Sanity check
        pipeline.sanityCheck()
        print("ST Pipeline, sanity check passed. Starting the run...")

        # Run the pipeline
        pipeline.run()
        print("ST Pipeline, run completed!")
    except Exception as e:
        print("Error running the pipeline")
        print(str(e))
        return 1
    finally:
        pipeline.clean_filenames()

    return 0


if __name__ == "__main__":
    sys.exit(main())
