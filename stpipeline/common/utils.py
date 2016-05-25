""" 
This file contains some general utilities
"""
import resource
import threading
from datetime import datetime
import os
import subprocess
import gc
from collections import namedtuple
_ntuple_diskusage = namedtuple('usage', 'total used free')

def which_program(program):
    """ 
    Checks that a program exists and is executable
    :param program: the program name
    :type program: str
    :returns: The program name if the program is in the system and is executable
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

def memory_use(point):
    """ 
    Returns memory usage at a certain time point
    :param point: a time point
    :type point: str
    :returns: a tuple (time, user memory, sys memory, total memory)
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           '''%(point,usage[0],usage[1],
                (usage[2]*resource.getpagesize()) / 1000000.0) 
           
           
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
            
def disk_usage(path):
    """
    Return disk usage statistics on the given path 
    :param path: the path to a folder
    :type path: str
    :returns: a tuple (total, used, free)
    """
    if not os.path.isdir(path):
        return
    st = os.statvfs(path)
    free = st.f_bavail * st.f_frsize
    total = st.f_blocks * st.f_frsize
    used = (st.f_blocks - st.f_bfree) * st.f_frsize
    return _ntuple_diskusage(total, used, free)

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
    In writing the file if it exists already and check for free space
    In reading mode check that the file exists
    :param filename: the path of the file
    :param atrib: the file open attribute
    :type filename: str
    :type atrib: str
    :returns: the file descriptor
    :raises: RuntimeError
    """
    if atrib.find("w") != -1:
        safeRemove(filename)
        usage = disk_usage(os.path.dirname(filename))
        if usage.free <= 4073741824: # at least 4GB
            raise IOError("Error, no free space available to open file %s\n" % (filename))
    elif atrib.find("r") != -1:
        if filename is None or not os.path.isfile(filename): # is it present?
            raise IOError("Error, the file does not exist %s\n" % (filename))
    else:
        raise IOError("Error, incorrect attribute %s\n" % (atrib))

    return open(filename, atrib)

def fileOk(_file):
    """
    Checks file exists and is not zero size
    """
    return _file is not None and os.path.isfile(_file) and not os.path.getsize(_file) == 0
    
def replaceExtension(filename,extension):
    """
    Replace the extension of filename 
    for the extension given, returns the new filename
    including the path 
    extension must be like .ext 
    """
    base = os.path.splitext(filename)[0]
    return base + extension
 
def stripExtension(filename):
    """
    Remove the extension from a file name
    :param filename: the file name
    :type filename: str
    :returns: the file without the extension
    """
    f = filename.rsplit('.', 1)
    if f[0].find("/") != -1:
        return f[0].rsplit('/', 1)[1]
    else:
        return f[0]

def getExtension(filename):
    """
    Gets the extension of a filename
    :param filename: the file name
    :type filename: str
    :returns: the extension of the filename
    """
    f = filename.rsplit('.', 1)
    return f[1]
        
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