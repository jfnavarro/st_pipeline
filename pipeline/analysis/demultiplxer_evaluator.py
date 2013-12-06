import sys

def file_reader(infilepath):
    """Both the solution file and the fastq file have the same  tab separated structure.
        Returns a dictionary  with read title as key and the real or assigned bc as the value."""
    res_dict={}
    infile=open(infilepath)
    for line in infile:
        #print line
        tmp=line.strip().split('\t')
        res_dict[tmp[0]]=tmp[1]
    infile.close()
    return res_dict

def performance_calculator(learned_bcs_dict,real_answer_dict):
    """ Returns the percentage of the correctly demultiplexed reads. """
    
    total_no_reads=len(learned_bcs_dict.keys())
    correctly_demultiplxed=0
    for read_title, bc in learned_bcs_dict.iteritems():
        if real_answer_dict[read_title]==bc:
            correctly_demultiplxed+=1
            
    performance=float(correctly_demultiplxed)/total_no_reads
    return performance

def main(answer_file_path, result_file_path):
    learned_bcs_dict=file_reader(result_file_path)
    real_answer_dict=file_reader(answer_file_path)
    res=performance_calculator(learned_bcs_dict, real_answer_dict)
    return res



if __name__=='__main__':
    #result_file_path='/Users/hosseinshahrabifarahani/Documents/ST_experiments/demultiplexing_experiments/mutated_fastq.fastq'
    result_file_path='/Users/hosseinshahrabifarahani/Documents/ST_experiments/input_files/raw_reads/miseqf1/F1-i1_S1_L001_R1_001_formated_withTranscript_nameMap.txt'
    answer_file_path='/Users/hosseinshahrabifarahani/Documents/ST_experiments/demultiplexing_experiments/solfile.txt'
    performance=main(answer_file_path,result_file_path)