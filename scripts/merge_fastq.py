#! /usr/bin/env python
""" 
Script that merges the FASTQ files present in an Illumina run folder path.
The script merges the the FASTQ files based on the indexes/identifiers given
as input and puts the merged files in the given output folder.

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys
import os
import glob
import shutil
import subprocess

def run_command(command, out=subprocess.PIPE):
    try:
        print "Running command: {}".format(" ".join(x for x in command).rstrip())
        proc = subprocess.Popen(command,
                                stdout=out, stderr=subprocess.PIPE,
                                close_fds=True, shell=False)
        (stdout, errmsg) = proc.communicate()
        print stdout
        print errmsg
    except Exception as e:
        raise e
               
def main(run_path, indexes, out_path):

    if not os.path.isdir(run_path) or not os.path.isdir(out_path):
        sys.stderr.write("Error, run path or output path folders do not exist\n")
        sys.exit(1)
                
    # First gunzip all the FASTQ files
    os.chdir(run_path)
    try:
        for file in glob.glob("*.gz"):
            run_command(["gunzip", "-f", file])
    except Exception as e:
        sys.stderr.write("Error, gunziping FASTQ files\n" + str(e))
        sys.exit(1)
       
    # Second merge the FASTQ files
    for index in indexes:
        r1_files = sorted(glob.glob("*{}*R1*.fastq".format(index)))
        r2_files = sorted(glob.glob("*{}*R2*.fastq".format(index)))
        try:
            with open("{}_R1.fastq".format(index), "w") as file1:
                run_command(["cat"] + r1_files, out=file1)
            with open("{}_R2.fastq".format(index), "w") as file2:
                run_command(["cat"] + r2_files, out=file2)
        except Exception as e:
            sys.stderr.write("Error, merging FASTQ files\n" + str(e))
            sys.exit(1)
    
    # Third gzip everything again
    try:
        for file in glob.glob("*.fastq"): 
            run_command(["gzip", "-f", file])
    except Exception as e:
        sys.stderr.write("Error, gziping FASTQ files\n" + str(e))
        sys.exit(1)
        
    # Move merged FASTQ files to output path
    if run_path != out_path:
        for index in indexes:
            for file in glob.glob("{}_R*".format(index)):
                shutil.move(file, out_path)
                   
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-path", required=True, help="Path to the run folder")
    parser.add_argument("--out-path", required=True, help="Path to the output folder")
    parser.add_argument("--identifiers", required=True, nargs='+', type=str,
                        help="List of identifiers for each sample (E.x. S1 S2 S3 S4)")
    args = parser.parse_args()
    main(args.run_path, args.identifiers, args.out_path)
