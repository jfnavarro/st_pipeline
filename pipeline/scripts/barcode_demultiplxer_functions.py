import time
import itertools
from optparse import OptionParser
from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from operator import itemgetter, attrgetter
from Bio.Blast.Applications import NcbiblastxCommandline


class Barcode:
    
    def __init__(self,seq,x,y):
        
        self.x = x
        self.y = y
        self.seq = seq
        
           
            
def OLD_demultiplex(id_filepath, reads, output_path ,  bc_offset=0, bc_len=27,alignment_method='h_distance_based',local_dist_measure='hamming',kmer_length=9,max_local_distance=3,max_global_distance=9):
    
    barcodes_id_dict=barcodefile_reader(id_filepath)
    res_list=[]
    # comparing each read with the rest of the barcodes
    for (title, seq, quality) in read_fastq(reads):
        #print title
        tmp=[]
        bc=seq[bc_offset:bc_offset+bc_len]
        
        for id in barcodes_id_dict.keys():
            
            if alignment_method=='h_distance_based':
                d=hamming_distance(bc, id)
                tmp.append((d,id))
                
            if alignment_method=='l_distance_based':
                d=levenshtein_distance(bc, id) 
                tmp.append((d,id))
                
            if alignment_method=='kmer_method':
                ids_kmer_dict=sliced_seq_maker(id_filepath,9)
                kmer_method(reads, barcodes_id_dict,ids_kmer_dict, kmer_length , bc_len, max_local_distance, max_global_distance, local_dist_measure,bc_offset=0)
                
        tmp.sort()
        nearest_dist=tmp[0][0]
        if nearest_dist <= max_distance:
            res_list.append(title+'\t'+id+'\t'+barcodes_id_dict[id][0]+'\t'+barcodes_id_dict[id][1]+'\t'+quality+'\t'+seq+'\n')

    return res_list
            
            
def demultiplex(id_filepath, reads, output_path ,  bc_offset=0, bc_len=27,alignment_method='h_distance_based',local_dist_measure='hamming',kmer_length=9,max_local_distance=3,max_global_distance=9):
    
    barcodes_id_dict=barcodefile_reader(id_filepath)
    res_list=[]
    # comparing each read with the rest of the barcodes
     
    if alignment_method=='kmer_method':
        ids_kmer_dict=sliced_seq_maker(id_filepath,9)
        res_list=kmer_method(reads, barcodes_id_dict,ids_kmer_dict, kmer_length , bc_len, max_local_distance, max_global_distance, local_dist_measure,bc_offset=0)
                
    return res_list

    
def extract_barcode_from_read(seq1, bc_len, bc_offset = 0):
    bc=seq1[bc_offset:bc_offset+bc_len]
    return bc
    
    
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


def give_distance(s1,s2,distance_type):
    if distance_type=='levenshtein':
        d=levenshtein_distance(s1, s2)
    if distance_type=='hamming':
        d=hamming_distance(s1, s2)
    return d

    
    
def read_fastq(fname):
    """Provide read info from fastq file, potentially not existing."""
    if fname:
        with open(fname) as in_handle:
            for info in FastqGeneralIterator(in_handle):
                yield info
    else:
        for info in itertools.repeat(("", None, None)):
            yield info
     
            
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


def kmer_producer(seq,kmer_len):
    """ gets string seq and length of the kmer and produces a 
        set of kmers with length kmer_len"""
        
    res_set=set([])
    for i in range(len(seq)-kmer_len+1):
        res_set.add(seq[i:i+kmer_len])
    return res_set


def barcode_extractor(read,barcode_length,barcode_start=0):
    barcode=read[barcode_start:barcode_start+barcode_length]
    return barcode


def barcodefile_reader(infile_path):
    """ reads the barcodes file and returns a dictionary of the barcodes and key and corrdinate tupple as value """
    res_dict={}
    infile=open(infile_path)
    for line in infile:
        tmp=line.strip().split('\t')
        res_dict[tmp[0]]=(tmp[1],tmp[2])
    return res_dict

def find_barcode_kmersearch(read_barcode, barcodes_dict):
    pass


def find_barcode_SW(barcode_read,barcodes_dict,gap_open_penalty,gap_extend_penalty,verbose=False):
    """uses Smith-Waterman dynamic prograkmming alignment"""
    # First searching for the best matches
    if barcode_read in barcodes_dict.keys():
        if(verbose): 
            print "Perfect match!" 
            return barcode_read,barcodes_dict[barcode_read].x,barcodes_dict[barcode_read].y
        
    # if no perfect match found pairwise alignment is done
    for bc in barcodes_dict.keys():
        alignments=pairwise2.align.localxs(bc,barcode_read,float(gap_open_penalty),float(gap_extend_penalty),one_alignment_only=True)
        matches_no=alignments[0][2]
        d=len(barcode_read) - matches_no  
        
    return d


def sliced_seq_maker(infile_path,kmer_length):
    """ slices the sequences in the id file"""
    res_dict={}
    infile=open(infile_path)
    for line in infile:
        tmp=line.strip().split('\t')
        seq=tmp[0]
        res_dict[seq]=[]
        for sliced in [seq[i:i+kmer_length] for i in range(0, len(seq),kmer_length)]:
            res_dict[seq].append(sliced)
    return res_dict


def back_kmer_comparer(kmerlist_1, kmerlist_2, min_local_dist, min_gloab_dist,distance_measure='levenshtein'):
    """ get 2 sets of kmers and measures the distance between them.
        If any of the local distances or the global distance is above a threshold it returns no_match """
    
    local_d_list=[]
    for i in range(len(kmerlist_1)):
        local_d=give_distance(kmerlist_1[i], kmerlist_2[i],distance_measure)
        #print kmerlist_1[i], kmerlist_2[i]
        if local_d > min_local_dist:
            r='no_match'
            break
        else:
            local_d_list.append(local_d)          
    glob_dist=sum(local_d_list)
    if glob_dist > min_gloab_dist:
        r='no_match' 
    
    return r


def kmer_comparer(kmerlist_1, kmerlist_2, min_local_dist, min_gloab_dist,distance_measure='levenshtein'):
    """ get 2 sets of kmers and measures the distance between them.
        If any of the local distances or the global distance is above a threshold it returns no_match """
        
    local_d_list=[]
    r=''
    for i in range(len(kmerlist_1)):
        local_d=give_distance(kmerlist_1[i], kmerlist_2[i],distance_measure)
        if local_d > min_local_dist:
            r='no_match'
            break
        else:
            local_d_list.append(local_d) 
                       
    if r!='no_match':
        glob_dist=sum(local_d_list)
        if glob_dist < min_gloab_dist:
            return glob_dist
        else:
            return 'no_match'
    else:
        return 'no_match'


def kmer_method(reads,ids_dict,ids_kmer_dict,kmer_length, bc_len,  min_local_dist, min_gloab_dist, dist_measure='levenshtein',bc_offset = 0):
    """ gets two sets of kmers from reads and ids and compares them.
        For each read it makes a list of distances from various ids and picks the nearest one. """
    res_list=[]
    for (title, seq, quality) in read_fastq(reads):
        #making a set of kmers from the read       
        score_list=[]
        read_kmer_list=[]
        bc=extract_barcode_from_read(seq,bc_len,bc_offset)
        for sliced in [bc[i:i+kmer_length] for i in range(0, len(bc),kmer_length)]:
            read_kmer_list.append(sliced)       
        # comparing two sets
        
        "each id"
        for k,v in ids_kmer_dict.iteritems():
            d= kmer_comparer(v, read_kmer_list, min_local_dist, min_gloab_dist, dist_measure)
            if d!='no_match':
                print d
                print title+'\t'+k+'\t'+ids_dict[k][0]+'\t'+ids_dict[k][1]+'\t'+quality+'\t'+seq+'\n'
                res_list.append(title+'\t'+k+'\t'+ids_dict[k][0]+'\t'+ids_dict[k][1]+'\t'+quality+'\t'+seq+'\n')
        return res_list
                
                
        
if __name__=='__main__':
    
    parser = OptionParser()
    parser.add_option("-i", "--ids_file", dest="path_to_ids_file")
    parser.add_option("-r", "--reads_file", dest="path_to_reads_file")
    parser.add_option("-k", "--kmer_length", dest="kmer_length", default=1)
    parser.add_option("-l", "--bc_length", dest="bc_length", default=27)
    parser.add_option("-e", "--max_local_dist", dest="max_local_dist", default=3)
    parser.add_option("-g", "--max_global_dist", dest="max_global_dist", default=9)
    parser.add_option("-b", "--bc_offset", dest="bc_offset", default=0)
    parser.add_option("-a", "--alignment_method", dest="alignment_method", default="kmer_based")
    parser.add_option("-m", "--dist_measure", dest="dist_measure", default="hamming")
    parser.add_option("-o", "--output_file", dest="output_file_path")
     
    options, args = parser.parse_args()
    demultiplex(options.path_to_ids_file, options.path_to_reads_file , options.output_file_path, int(options.bc_offset) , int(options.bc_length) , \
                options.alignment_method, options.dist_measure, int(options.kmer_length) , int(options.max_local_dist), int(options.max_global_dist))
     
    demultiplex(options.path_to_ids_file, options.path_to_reads_file , options.output_file_path, 0, 27, 'kmer_method', 'hamming', 9, 7, 12)
    
    
    
    
    

    
    
#     """ lines for testing pre pars"""
#     reads='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/raw_reads/miseqf1/F1-i1_S1_L001_R1_001_formated_withTranscript.fastq'
#     reads='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/raw_reads/miseqf1/1_sample.fastq'
#     id_filepath='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/ids/130307_Design3_27mer.txt'
#     output_path=''
#     demultiplex(id_filepath, reads, output_path, 0, 27, 'kmer_method', 'hamming', 9, 7, 12)
    


         

