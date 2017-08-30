""" 
This file contains some general utilities
"""
import threading
from datetime import datetime
import os
import subprocess
import stat

def which_program(program):
    """ 
    Checks that a program exists and is executable
    :param program: the program name
    :type program: str
    :return: The program name if the program is in the system and is executable
    """
    def is_exe(fpath):
        return fpath is not None and os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate
    return None
             
class TimeStamper(object):
    """
    Thread safe time stamper 
    """
    def __init__(self):
        self.lock = threading.Lock()
        self.prev = None
        self.count = 0

    def getTimestamp(self):
        with self.lock:
            ts = datetime.now()
            if ts == self.prev:
                ts += '.%04d' % self.count
                self.count += 1
            else:
                self.prev = ts
                self.count = 1
        return ts

def safeRemove(filename):
    """
    Safely remove a file
    :param filename: the path of the file
    :type filename: str
    """
    try:
        if filename is not None and os.path.isfile(filename):
            os.remove(filename)
    except UnboundLocalError:
        pass
        
def safeOpenFile(filename, atrib):
    """
    Safely opens a file
    For writing mode it removes the previous file if it exits
    For reading mode it check that the file exists
    :param filename: the path of the file
    :param atrib: the file open/write attribute
    :type filename: str
    :type atrib: str
    :return: the file descriptor
    :raises: IOError
    """
    if filename is None:
        raise IOError("Error, no valid file name given\n")
    if atrib.find("w") != -1:
        safeRemove(filename)
    elif atrib.find("r") != -1:
        if not (os.path.isfile(filename) or is_fifo(filename)):
            raise IOError("Error, the file does not exist {}\n".format(filename))
    else:
        raise IOError("Error, incorrect attribute {}\n".format(atrib))

    return open(filename, atrib)

def fileOk(_file):
    """
    Checks file exists and is not zero size
    :param file: a file name
    :return: True if the file is correct
    """
    return _file is not None and os.path.isfile(_file) and not os.path.getsize(_file) == 0
        
def getSTARVersion():
    """
    Tries to find the STAR binary
    and makes a system call to get its
    version and return it
    """
    version = ""
    try:
        proc = subprocess.Popen(["STAR", "--version"], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE,
                                shell=False, close_fds=True)
        (stdout, errmsg) = proc.communicate()
        version = stdout
    except Exception:
        version = "Not available"
    return version.rstrip()

def getTaggdCountVersion():
    """
    Tries to find the Taggd binary
    and makes a system call to get its
    version and return it
    """
    version = ""
    try:
        proc = subprocess.Popen(["pip", "show", "taggd"], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE,
                                shell=False, close_fds=True)
        (stdout, errmsg) = proc.communicate()
        for line in stdout.split("\n"):
            if line.find("Version:") != -1:
                version = str(line.split()[-1])
    except Exception:
        version = "Not available"
    return version.rstrip()

def getHTSeqCountVersion():
    """
    Tries to find the HTSeqCount binary
    and makes a system call to get its
    version and return it
    """
    version = ""
    try:
        proc = subprocess.Popen(["pip", "show", "htseq"], 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE,
                                shell=False, close_fds=True)
        (stdout, errmsg) = proc.communicate()
        for line in stdout.split("\n"):
            if line.find("Version:") != -1:
                version = str(line.split()[-1])
    except Exception:
        version = "Not available"
    return version.rstrip()

def is_fifo(file_name):
    """
    Checks if the file name is a FIFO
    :param file_name: a file name
    :return: True if the file is a FIFO
    """
    return (os.path.exists(file_name) and stat.S_ISFIFO(os.stat(file_name).st_mode))