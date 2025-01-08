#! /usr/bin/env python
""" 
Script that merges the FASTQ files present in an Illumina run folder path.
The script merges the FASTQ files based on the identifiers (typically indexes) given
as input (one identifier to each sample) and puts the merged files in the given output folder.

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys
import os
import glob
import shutil
import subprocess
from typing import List, Union, IO

def run_command(command: List[str], out: Union[int, IO[bytes]] = subprocess.PIPE) -> None:
    """
    Executes a shell command and prints its stdout and stderr.

    Args:
        command: The command to execute, represented as a list of strings.
        out: The output stream for the command's stdout. Defaults to subprocess.PIPE.

    Raises:
        Exception: If an error occurs during command execution.
    """
    try:
        print(f"Running command: {" ".join(x for x in command).rstrip()}")
        proc = subprocess.Popen(
            command,
            stdout=out,
            stderr=subprocess.PIPE,
            close_fds=True,
            shell=False
        )
        stdout, errmsg = proc.communicate()
        print(stdout)
        print(errmsg)
    except Exception as e:
        raise e


def main(run_path: str, indexes: List[str], out_path: str) -> int:
    if not os.path.isdir(run_path) or not os.path.isdir(out_path):
        print("Error, either run_path or out_path folders do not exist")
        return 1

    # First gunzip all the FASTQ files
    os.chdir(run_path)
    for file in glob.glob("*.gz"):
        try:
            run_command(["gunzip", "-f", file])
        except Exception as e:
            print(f"Error, gunziping FASTQ file {file}, {str(e)}")
            return 1

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
            print(f"Error, merging FASTQ files, {str(e)}")
            return 1

    # Third gzip everything again
    for file in glob.glob("*.fastq"):
        try:
            run_command(["gzip", "-f", file])
        except Exception as e:
            print(f"Error, gziping FASTQ file {file}, {str(e)}")
            return 1

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
    sys.exit(main(args.run_path, args.identifiers, args.out_path))
