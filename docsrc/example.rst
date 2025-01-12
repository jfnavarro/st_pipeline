Examples
--------

The following is an example of an BASH file to run the ST pipeline.

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

If you want to process Visium datasets it is recommended to use these settings


.. code-block:: bash

    --allowed-missed 1 \
    --allowed-kmer 4 \
    --umi-allowed-mismatches 2 \
    --umi-start-position 16 \
    --umi-end-position 28 \
