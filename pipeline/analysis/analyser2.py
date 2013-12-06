import json
import numpy as np
import matplotlib.pyplot as plt
import os


def jason_reader(data_path):
    # reads the json file and returns a list of reads numbers
    out_list=[]
    json_data=open(data_path).read()
    data = json.loads(json_data)
#     for cell in data:
#         print 'x', cell['x']
#         print 'y', cell['y']
#         print 'hits',cell['hits']
#         print 'gene',cell['gene']
#         print '---------------'
    return data
        
        
def dimention_finder(data):
    #gets data from json and finds the maximum x and y dimentions
    max_x=0
    max_y=0  
    genes_set=set([])
    for cell in data:
        if int(cell['x'])>max_x:
            max_x=int(cell['x'])
        if int(cell['y'])>max_y:
            max_y=int(cell['y'])    
        if cell['gene'] not in genes_set:
            genes_set.add(cell['gene'])       
    return max_x,max_y,genes_set

def np_array_maker(data,max_x,max_y,genes_set):
    genes_exp_dict={}
    
    #initializing the dictionary of exp. matrices
    for gene in genes_set:
        genes_exp_dict[gene]=np.array(np.zeros((max_x+1,max_y+1)))
        
    #populating the expression matrices
    for cell in data:
        genes_exp_dict[cell['gene']][int(cell['x']),int(cell['y'])]=int(cell['hits'])
        #print cell['x'],cell['y']
        
    return genes_exp_dict

def file_writer(genes_exp_dict,output_folder):
    #writes the expression matrix of each gene in a text file
    for k,v in genes_exp_dict.iteritems():
        print k
        np.savetxt(output_folder+'/'+k+'.csv',v, delimiter=' ,', newline='\n',fmt='%10.0f')
        
        
def accumulative_array_maker(data,max_x,max_y,output_path):
    
    #initializing the expression array
    exp_array=np.array(np.zeros((max_x+1,max_y+1)))
    
    for cell in data:
        exp_array[(int(cell['x']),int(cell['y']))]+=int(cell['hits'])
    
    np.savetxt(output_path,exp_array, delimiter=' ,', newline='\n',fmt='%10.0f')
    
def max_hit_finder(data):
    
    max_hit=0
    max_barcode=''
    for cell in data:
        if int(cell['hits']) > max_hit:
            max_hit=int(cell['hits'])
            max_barcode=cell['barcode']
            max_gene=cell['gene']
    return max_hit,max_barcode,max_gene
        
if __name__=='__main__':
    results_folder='/Users/hosseinshahrabifarahani/Documents/ST_experiments/findindexes_km_experiments/miseqf1/m_k_experiments'
     
    for file in os.listdir(results_folder):
        #print file
        if file.split('.')[-1]=='json' and file.find('barcodes')!=-1:
            data_path=results_folder+'/'+file
            output_path=results_folder+'/'+file.split('exp')[0]+'csv'
            print output_path
            print file
            data=jason_reader(data_path)
            x,y,genes=dimention_finder(data)
            accumulative_array_maker(data,x,y,output_path)
            
            
            
    
    
    
    
#     data_path='/Users/hosseinshahrabifarahani/Documents/ST_experiments/findindexes_km_experiments/miseqf1/m_k_experiments/m_5_k7.exp_barcodes.json'
#     output_path='/Users/hosseinshahrabifarahani/Documents/ST_experiments/findindexes_km_experiments/miseqf1/m_k_experiments'
#     data = jason_reader(data_path)
#  
#     x,y,genes=dimention_finder(data)
#     accumulative_array_maker(data,x,y,output_path)
# 
# 
#     mat_dict=np_array_maker(data,x,y,genes)
#     file_writer(mat_dict,output_path)
    