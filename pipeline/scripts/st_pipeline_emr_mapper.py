#!/usr/bin/env python
""" 
This is the Amazon EMR wrapper for the mapper function to run the Map Reduce version of the ST pipeline
"""

import sys
sys.path.append("./")
from stpipeline.core.pipeline import *
from stpipeline.common.json_utils import json_iterator
import argparse
 
def createFastqFilesFromChunk(chunk):
    #create two unique temp names
    temp_name = tempfile.mktemp(prefix='stpipeline_temp_', suffix=str(random.random()), dir='')
    fastq_forward = temp_name + "_forward.fastq"
    fastq_reverse = temp_name + "_reverse.fastq"
    
    #create two fastq writers
    outF = safeOpenFile(fastq_forward, 'w')
    outF_writer = writefq(outF)
    outF2 = safeOpenFile(fastq_reverse, 'w')
    outF_writer2 = writefq(outF2)
    
    #interate the lines in the chunk
    #each line has a fastq record (forward and reverse)
    for line in chunk.split("\n"):
        cols = line.replace("\t", " ").split(" ")
        assert len(cols) == 5
        header = cols[0]
        seq1 = cols[1]
        qual1 = cols[2]
        seq2 = cols[3]
        qual2 = cols[4]
        outF_writer.send((header, seq1, qual1))
        outF_writer2.send((header, seq2, qual2))
    outF.close()
    outF2.close()
    outF_writer.close()
    outF_writer2.close()
    
    return fastq_forward, fastq_reverse, temp_name
            
            
def run_pipeline_on_chunks(pipeline, chunks):
    """ 
    This function is called for Map Reduce jobs, when we want to run the pipeline,
    once all the streaming has been done
    in the input, it will iterate trough the chunks, create temp fastq files
    and call the pipeline on them, it will then parse the output and send
    the json formated features to the reducer
    """
    
    #TODO do mapping using gene as KEY
        
    for chunk in chunks:
        #create fastq files from tab delimited chunk
        reverse, forward, temp_name = createFastqFilesFromChunk(chunk)
        #tell the pipeline the location of the files
        pipeline.Fastq_fw = os.path.abspath(reverse)
        pipeline.Fastq_rv = os.path.abspath(forward)
        pipeline.expName = temp_name
        #run pipeline with new files and name 
        try:
            pipeline.sanityCheck()
            pipeline.run()
        except Exception as e:
            print e
            # remove temp files
            safeRemove(reverse)
            safeRemove(forward)
            sys.stderr.write("Error: Running the pipeline \n")
            sys.exit(-1)
                        
        # now we parse the output of pipeline (json) 
        outputPipeline = temp_name + "_barcodes.json"  
        outputPipelineReads = temp_name + "_reads.json"  
        it = json_iterator(outputPipeline)
        # send features to the reducer jobs (reads are ignored for now)
        for doc in it:
            feature_gene = (doc['y'], doc['x'], doc['gene'], doc['barcode'])
            hits = doc['hits']
            doc_json_formated = {}
            doc_json_formated['y'], doc_json_formated['x'], doc_json_formated['gene'], doc_json_formated['barcode'] = feature_gene
            #TODO would be nice to be able to send json objects trough
            yield doc_json_formated,hits
                
        # remove temp files
        safeRemove(reverse)
        safeRemove(forward)
        safeRemove(outputPipeline)
        safeRemove(outputPipelineReads)   
                             
def main(argv):
    
    pipeline = Pipeline()
    
    #create a parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser = pipeline.createParameters(parser)
    
    #parse arguments
    options = parser.parse_args()
    
    #load parameters
    pipeline.load_parameters(options)
    pipeline.createLogger()

    #local variables
    batch = []
    chunks = [] 
        
    #parse input streaming to create chunks
    line = sys.stdin.readline()
    try:
        while line:
            batch.append(line + '\n')
            if len(batch) == int(options.chunks): 
                chunks.append(''.join(batch))
                batch = []
            line = sys.stdin.readline()
    except "end of file":
        return None
    
    #some lines left behind?
    if len(batch) > 0:
        chunks.append(''.join(batch))
        
    #call pipeline and output streams
    for key,value in run_pipeline_on_chunks(pipeline, chunks):
        print '%s%s%d' % (key, "\t", value)
    
if __name__ == "__main__":
    main(sys.argv)
