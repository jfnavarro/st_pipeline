#! /usr/bin/env python
""" 
Script that merges the BCL files in an Illumina run folder path.
The script merges the the BCL files based on the indexes given
as input and puts the merged files in the given output folder.
@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import argparse
import sys
import os
import glob
import shutil
import subprocess

def run_command(command, out=subprocess.PIPE):
    try:
        print "running command: {}".format(" ".join(x for x in command).rstrip())
        proc = subprocess.Popen(command,
                                stdout=out, stderr=subprocess.PIPE,
                                close_fds=True, shell=False)
        (stdout, errmsg) = proc.communicate()
        print errmsg
    except Exception as e:
        print str(e)
        raise e
               
def main(run_path, indexes, out_path):

    if not os.path.isdir(run_path) or not os.path.isdir(out_path):
        sys.stderr.write("Error, run path or output path folders do not exist\n")
        sys.exit(1)
    
    run_path = os.path.join(run_path, "Demultiplexing")
    if not os.path.isdir(run_path):
        sys.stderr.write("Error, run path does not contain a Demultiplexing folder\n")
        sys.exit(1)
                
    # First gunzip all the bcl files
    os.chdir(run_path)
    try:
        for file in glob.glob("*.gz"):
            run_command(["gunzip", "-f", file])
    except Exception as e:
        sys.stderr.write("Error, gunziping BCL files\n")
        sys.exit(1)
       
    # Second merge the BCL files
    for index in indexes:
        r1_files = sorted(glob.glob("*{}*R1*.fastq".format(index)))
        r2_files = sorted(glob.glob("*{}*R2*.fastq".format(index)))
        try:
            with open("{}_R1.fastq".format(index), "w") as file1:
                run_command(["cat"] + r1_files, out=file1)
            with open("{}_R2.fastq".format(index), "w") as file2:
                run_command(["cat"] + r2_files, out=file2)
        except Exception as e:
            sys.stderr.write("Error, merging BCL files\n")
            sys.exit(1)
    
    # Third gzip everything again
    try:
        for file in glob.glob("*.fastq"): 
            run_command(["gzip", "-f", file])
    except Exception as e:
        sys.stderr.write("Error, gziping BCL files\n")
        sys.exit(1)
        
    # Move merged BCL files to output path
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
