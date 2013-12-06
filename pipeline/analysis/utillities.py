
import sys
import os
import itertools

# this file contains some generic necessary functions
class Barcode:
    
    def __init__(self,seq,x,y):
        
        self.x = x
        self.y = y
        self.seq = seq
        
        
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

def kmer_producer(seq,kmer_len):
    """ gets string seq and length of the kmer and produces a 
        set of kmers with length kmer_len"""
        
    res_set=set([])
    for i in range(len(seq)-kmer_len+1):
        res_set.add(seq[i:i+kmer_len])
    return res_set
          
# for each id a set of reads and their accumulative distance, picking the nearest one and doing full
# Smith-Waterman alignment. 
def kmer_comparer_2(kmerSet_1,kmerSet_2,max_mismatch,distance_type):
    """ gets two sets of kmers and comparers them
        returns the sum of distances and max distance"""
    
    dist_sum=0
    max_dist=0
    for k1 in kmerSet_1:
        for k2 in kmerSet_2:
            if distance_type=='hamming':
                distance=hamming_distance(k1, k2)
            if distance_type=='levenshtein':
                distance=levenshtein_distance(k1, k2)
            dist_sum+=distance
            if distance > max_dist:
                max_dist=distance
    return dist_sum,max_dist


def kmer_comparer(kmerSet_1,kmerSet_2,max_mismatch,distance_type):
    dist_sum=0
    max_dist=0
    for k1 in kmerSet_1:
        for k2 in kmerSet_2:
            if distance_type=='hamming':
                distance=hamming_distance(k1,k2)
            if distance_type=='levenshtein':
                distance=levenshtein_distance(k1, k2)
            dist_sum+=distance
            if distance > max_dist:
                max_dist=distance
    return dist_sum, max_dist
    

def best_match(reads_set, barcodes,max_mismatch,barcode_len,kmer_len=1):
    
    
    for read in reads_set:
        bc= extract_barcode_from_read(read, barcode_len)
        kmer_producer(bc, kmer_len)
        
    

def extract_barcode_from_read(seq1, bc_offset = 0):
    
    """Function which pulls a barcode of a provided size from paired seqs.

    This respects the provided details about location of the barcode, returning
    items of the specified size to check against the read.
    
    """
    def _get_end(size):
        assert size > 0
        return seq1[bc_offset:size+bc_offset]
    return _get_end


    
              
def read_barcodes(fname,kmer = 1):
    
    """ exctracts barcodes and cordinate from id file, 
    creates two maps one with full barcodes and another with sliced barcodes into kmer """
    
    barcodes = {}
    barcodes_kmer = {}
    with open(fname) as in_handle:
        for line in (l for l in in_handle if not l.startswith("#")):
            seq, x, y = line.rstrip("\t").split()
            barcodes[seq] = Barcode(seq,x,y)
            for sliced in [seq[i:i+kmer] for i in range(0, len(seq),kmer)]:
                if(len(sliced) >= kmer):
                    barcodes_kmer[sliced] = Barcode(seq,x,y)
            
    return barcodes,barcodes_kmer
           
           
                
def _barcode_has_ambiguous(barcodes):
    
    """ consider replace N for T """
    
    for seq in barcodes.keys():
        if "N" in seq:
            return True
    return False
    


if __name__=='__main__':
    seq='0123456789'
    kmer_len=4
    mers=kmer_producer(seq, kmer_len)
