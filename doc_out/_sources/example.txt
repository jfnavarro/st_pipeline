Examples
--------

The following is an example of an BASH file to run the ST pipeline. 
This example is for version 1.3.2 of the ST pipeline.

.. code-block:: bash

	#!/bin/bash

	# FASTQ reads
	FW=/YOUR_RUN/R1.fastq.gz
	RV=/YOUR_RUN/R2.fastq.gz

	# References for mapping, annotation and nonRNA-filtering
	MAP=/mouse/GRCm38_86v2/StarIndex
	ANN=/mouse/GRCm38_86v2/annotation/annotation.gtf
	CONT=/mouse/GRCm38_86v2/ncRNA/StarIndex

	# Barcodes settings
	ID=/stpipeline/ids/YOUR_IDs.txt

	# Output folder and experiment name
	# Do not use / or \ in the experiment name
	OUTPUT=/your_experiment_folder
	EXP=YOUR_EXP_NAME

	# Running the pipeline
	st_pipeline_run.py \
	  --output-folder $OUTPUT \
	  --ids $ID \
	  --ref-map $MAP \
	  --ref-annotation $ANN \
	  --expName $EXP \
	  --htseq-no-ambiguous \
	  --verbose \
	  --log-file $OUTPUT/${EXP}_log.txt \
	  --contaminant-index $CONT \
	  $FW $RV
