Examples
--------

The following is an example of an SBATCH file to run the ST pipeline. An SBATCH
file is used with the SLURM job scheduler.
This example is for version 1.3.0 of the ST pipeline.

.. code-block:: bash

	#!/bin/bash

	# FASTQ reads
	FW=/fastdisk/INBOX/YOUR_RUN/R1.fastq.gz
	RV=/fastdisk/INBOX/YOUR_RUN/R2.fastq.gz

	# References for mapping, annotation and nonRNA-filtering
	MAP=/fastdisk/mouse/GRCm38_86v2/StarIndex
	ANN=/fastdisk/mouse/GRCm38_86v2/annotation/annotation.gtf
	CONT=/fastdisk/mouse/GRCm38_86v2/ncRNA/StarIndex

	# Barcodes settings
	ID=/fastdisk/ids/YOUR_IDs.txt

	# Output folder and experiment name
	# Do not use / or \ in the experiment name
	OUTPUT=/your_experiment_folder
	EXP=YOUR_EXP_NAME

	# Running the pipeline
	# Add this if you want to keep the intermediate files  --no-clean-up
	# Add this if you want to keep discarded files --keep-discarded-files
	st_pipeline_run.py \
	  --output-folder $OUTPUT \
	  --ids $ID \
	  --ref-map $MAP \
	  --ref-annotation $ANN \
	  --expName $EXP \
	  --remove-polyA 10 \
	  --remove-polyT 10 \
	  --remove-polyG 10 \
	  --remove-polyC 10 \
	  --htseq-no-ambiguous \
	  --verbose \
	  --mapping-threads 8 \
	  --log-file $OUTPUT/${EXP}_log.txt \
	  --two-pass-mode \
	  --filter-AT-content 90 \
	  --filter-GC-content 90 \
	  --contaminant-index $CONT \
	  --min-length-qual-trimming 30 \
	  --disable-clipping \
	  $FW $RV
