The scripts in the analysis folder.

barcode_distance_stat.py
Gets a the list of barcode ids as input, measures the hamming, levenshtein, or smith-waterman distance between then and 
outputs the statistics of distance between them.

------------------------------
synthetic_data_generator.py
generates synthetic fastq files for measuring the performance of the demultiplxer.

Inputs: 
	Source fastq, it can be any fastq file. It sequences accesions and other parts are used for generating the test file.
	Source designed barcodes
	Number of desired test reads
	barcode length
	distance between the generated barcodes.
	
Outputs:
	Solution file, which is a csv file that contains the real barcode that is attached to each sequence.
	Generated fastq file
-------------------------------
barcode_demultiplxer_functions.py
	Uses seed matching blast like method for demultiplxeing barcodes

Inputs:
	ids file
	reads file in fastq
	kmer length
	barcode length
	barcode offset
	distance measure, can be either hamming or levenshtein
	
Output:
	outputfile 
	
-------------------------------
demultiplxer_evaluator.py
	Evaluates the performance of the demultiplxer function.
	
Input: 
	output from the demultiplxe
	The txt file from the synthetic data generator that contains the true barcode for each read.
	
Output:
	The percentage of the correctly demultiplxed barcodes.
	
-------------------------------
analyser2.py
	Gets the json file and makes matrices for each separate gene in the feature and accumulative expression of all genes. 
	
	
	
	
	
