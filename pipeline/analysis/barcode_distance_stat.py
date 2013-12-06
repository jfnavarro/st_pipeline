import random
import sys
import itertools
from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def id_file_reader(filepath):
    infile=open(filepath)
    ids_list=[]
    for line in infile:
        id=line.strip().split('\t')[0]
        ids_list.append(id)
    return ids_list


def min_dist_finder(ids_list, sample_no,distance_type):
    #gets a list of ids and sample no and finds the distance of the nearest barcode to that sample
    sample_barcodes=random.sample(ids_list,sample_no)
    barcode_len=len(sample_barcodes[0])
    d_list=[]
        
    for sample_bc in sample_barcodes:
        min_distance=1000
        for bc in ids_list:
            if bc!=sample_bc:
                if distance_type=='SW':
                    alignments = pairwise2.align.globalmx(bc,sample_bc,1,0)
                    matches_no=alignments[0][2]
                    d=barcode_len - matches_no   
                if distance_type=='hamming':
                    d= hamming_distance(bc,sample_bc)
                if distance_type=='levenshtein':
                    d= levenshtein_distance(bc,sample_bc)
                if d < min_distance:
                    min_distance=d
        d_list.append(min_distance)
                    #print alignments
    return d_list

             
def outfile_writer(outfile_path,d_list):
    outfile=open(outfile_path,'w')
    for d in d_list:
        outfile.write(str(d)+',')
    outfile.write('\n')
    outfile.write('minimum distance is '+str(min(d_list)))
    outfile.close()
    

def hamming1(str1, str2):
    return sum(itertools.imap(str.__ne__, str1, str2))
  
  
def hamming_distance(s1, s2):
    "Return the Hamming distance between equal-length sequences."
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def levenshtein_distance(seq1, seq2):
    # measures the Levenstein distance between two sequences
    oneago = None
    thisrow = range(1, len(seq2) + 1) + [0]
    for x in xrange(len(seq1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in xrange(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(seq2) - 1]
    

if __name__=='__main__':
    filepath=sys.argv[1]
    sample_no=int(sys.argv[2])
    output_path=sys.argv[3]
    
    ids_list=id_file_reader(filepath)
    d_list=min_dist_finder(ids_list, sample_no)
    outfile_writer(output_path, d_list)
    
    
    
    #/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/ids/130208_Design4_24mer.txt
#     filepath='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/ids/130307_Design3_27mer.txt'
#     sample_no=4
#     ids_list=id_file_reader(filepath)
#     d_list=min_dist_finder(ids_list, sample_no)
#     
    print d_list
    print 'minimum distance is '+str(min(d_list))

#     alignments = pairwise2.align.globalxx("ACTG", "ACTA")
#     print alignments
#     print alignments[0][2]
