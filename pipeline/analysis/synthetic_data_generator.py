import sys
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import random 

def barcodefile_reader(infile_path):
    """ reads the barcodes file and returns a dictionary of the barcodes and key and corrdinate tupple as value """
    res_dict={}
    infile=open(infile_path)
    for line in infile:
        tmp=line.strip().split('\t')
        res_dict[tmp[0]]=(tmp[1],tmp[2])
    return res_dict



def read_fastq(fname):
    """Provide read info from fastq file, potentially not existing."""
    if fname:
        with open(fname) as in_handle:
            for info in FastqGeneralIterator(in_handle):
                yield info
    else:
        for info in itertools.repeat(("", None, None)):
            yield info
         
         
def mutate_barcode(bc,distance,dist_range='up_to'):
    """ gets a barcode and mutates it to a desired distance.
        This function can be improved to allow for levenshtein distance as well (padding). 
        The parameter dist_range takes two value up_to means distance is in range [0,distance], only means all bcs are at the same distance 
        """
    
    nucls=['A','C','T','G']
    #print bc
    l= [char for char in bc]
    #print l
    for i in range(distance):
        mut_pos=random.randint(0,len(bc)-1) 
        m = random.sample(nucls,1)[0]
        # making sure that the nucleotide always mutates to something different
        if dist_range=='only':
            while m == l[mut_pos]:
                m = random.sample(nucls,1)[0]
        if dist_range=='up_to':
            m = random.sample(nucls,1)[0]
        l[mut_pos]=m
            
    # converting the list to a string again
    res=''
    for nu in l:
#         print l
#         print nu
        res+=nu
    return res

       
def test_fastq_generator(source_fastq,test_reads_no,output_fastq_file,ids_file_path,distance,bc_length,solution_file):
    """ for a desired number of tests, attaches a mutated barcode to the existing fastq sequence. """
    ids_dict=barcodefile_reader(ids_file_path)
    outfile=open(output_fastq_file,'w')
    sol_file=open(solution_file,'w')
    
    i=0
    for (title, seq, quality) in read_fastq(source_fastq):
        i+=1
        if i > int(test_reads_no):
            break
        
        seq_no_bc=seq[bc_length:]
        # choosing a random bc_id, mutate it, and attach it to the read
        id=random.sample(ids_dict.keys(),1)[0]
        #mutate the chosen id
        mutated_id=mutate_barcode(id, distance)
        
        #writing the resulting fastq file
        outfile.write('@'+title+'\n')
        outfile.write(mutated_id+seq_no_bc+'\n')
        outfile.write('+\n')
        outfile.write(quality+'\n')
        
        #writing the solution file
        sol_file.write('@'+title+'\t'+id+'\n')
        
    outfile.close()
    sol_file.close()
    
    
if __name__=='__main__':
    source_fastq='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/raw_reads/miseqf1/F1-i1_S1_L001_R1_001_formated_withTranscript.fastq'
    ids_file_path='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/ids/130307_Design3_27mer.txt'
    output_fastq_file='/Users/hosseinshahrabifarahani/Documents/ST_experiments/demultiplexing_experiments/mutated_fastq.fastq'
    solution_file='/Users/hosseinshahrabifarahani/Documents/ST_experiments/demultiplexing_experiments/solfile.txt'
    test_reads_no=100
    bc_length=27
    distance=3
    
    test_fastq_generator(source_fastq,test_reads_no,output_fastq_file,ids_file_path,distance,bc_length,solution_file)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    