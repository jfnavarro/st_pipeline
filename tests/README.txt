*** Copyright Spatial Transcriptomics, 2014. ***

This folder contains test data for the pipeline.
Reads were extracted from a MOB tissue dataset and the genome version is GRCm38 with
ENSEMBLE annotation file and ArrayJet 1002 chip.

- A repeat-unmasked version of the shortest mouse chromosome, 19, was obtained from Ensembl.
- A STAR, 'chromosome19', was built for the chromosome.
- All reads that were mapped against chromosome 19 *for the original dataset* were retained as new input.
- The HTSeq original annotation file was stripped of all other content than chromosome 19
- The chip file (150204_arrayjet_1000L2_probes.txt) was unaltered, as was so the contaminant genome (R45S5_R5S1).



