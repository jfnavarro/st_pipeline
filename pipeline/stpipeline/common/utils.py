#!/usr/bin/env python
""" 
This file contains general utils and some file utils
"""

import resource
import threading
from datetime import datetime
import os
from collections import namedtuple
_ntuple_diskusage = namedtuple('usage', 'total used free')

def which(program):
    """ 
    check that a program exists and is executable 
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

def Using(point):
    """ 
    returns memory usage at a certain point 
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           '''%(point,usage[0],usage[1],
                (usage[2]*resource.getpagesize()) / 1000000.0) 
           
           
class TimeStamper(object):
    """
    thread safe time stamper 
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
    Return disk usage statistics about the given path 
    """   
    st = os.statvfs(path)
    free = st.f_bavail * st.f_frsize
    total = st.f_blocks * st.f_frsize
    used = (st.f_blocks - st.f_bfree) * st.f_frsize
    return _ntuple_diskusage(total, used, free)

def safeRemove(filename):
    """
    safely remove a file 
    """
    try:
        if filename is not None and os.path.isfile(filename):
            os.remove(filename)
    except UnboundLocalError:
        pass
        
def safeOpenFile(filename, atrib):
    """
    safely opens a file 
    """
    if atrib.find("w") != -1:
        safeRemove(filename)
        usage = disk_usage('/')
        if(usage.free <= 1073741824): ## at least 1GB
            raise RuntimeError("Error : no free space available\n")
    elif(atrib.find("r") != -1):
        if filename is None or not os.path.isfile(filename):  # is present?
            raise RuntimeError("Error : wrong filename\n")
    else:
        raise RuntimeError("Error : wrong attribute " + atrib + " opening file\n")

    handler = open(filename, atrib)
    return handler


def fileOk(_file):
    """
    checks file exists and is not zero size
    """
    return _file is not None and os.path.isfile(_file) and not os.path.getsize(_file) == 0
    
def replaceExtension(filename,extension):
    """
    replace the extesion of filename 
    for the extension given, returns the new filename
    including the path 
    extension must be like .ext 
    """
    base = os.path.splitext(filename)[0]
    return base + extension
 
def stripExtension(string):
    """
    remove the extension from string
    and returns it
    """
    f = string.rsplit('.', 1)
    if(f[0].find("/") != -1):
        return f[0].rsplit('/', 1)[1]
    else:
        return f[0]

def getExtension(string):
    """
    gets the filename extension of a filename
    """
    f = string.rsplit('.', 1)
    return f[1]

def getCleanFileName(path):
    """
    extracts and returns the filename from a complete path 
    """
    head, tail = os.path.split(path)
    return tail
