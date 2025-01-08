""" 
This file contains some general utilities
"""
import threading
from datetime import datetime
import os
import subprocess
from typing import Optional


def which_program(program: str) -> Optional[str]:
    """
    Checks if a program exists and is executable.

    Args:
        program: The program name.

    Returns:
        The full path to the program if found, otherwise None.
    """

    def is_exe(fpath: str) -> bool:
        return fpath and os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath: str):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, _ = os.path.split(program)
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

class TimeStamper:
    """
    Thread-safe time stamper to generate unique timestamps.
    """

    def __init__(self):
        self.lock = threading.Lock()
        self.prev: Optional[str] = None
        self.count: int = 0

    def get_timestamp(self) -> str:
        """
        Generates a unique timestamp.

        Returns:
            A unique timestamp.
        """
        with self.lock:
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
            if ts == self.prev:
                ts += f".{self.count:04d}"
                self.count += 1
            else:
                self.prev = ts
                self.count = 1
        return ts


def safe_remove(filename: Optional[str]) -> None:
    """
    Safely removes a file if it exists.

    Args:
        filename: Path to the file.
    """
    if filename and os.path.isfile(filename):
        try:
            os.remove(filename)
        except Exception as e:
            print(f"Error removing file {filename}: {e}")


def safe_open_file(filename: str, mode: str):
    """
    Safely opens a file.

    For write mode, removes the previous file if it exists.
    For read mode, checks that the file exists.

    Args:
        filename: Path to the file.
        mode: File open mode.

    Returns:
        The opened file descriptor.

    Raises:
        IOError: If the file does not exist for read mode or invalid mode is provided.
    """
    if "w" in mode:
        safe_remove(filename)
    elif "r" in mode and not os.path.isfile(filename):
        raise IOError(f"Error: File does not exist: {filename}")
    else:
        raise IOError(f"Error: Invalid mode: {mode}")

    return open(filename, mode)


def file_ok(file: Optional[str]) -> bool:
    """
    Checks if a file exists and is not empty.

    Args:
        file: Path to the file.

    Returns:
        bool: True if the file exists and is not empty, otherwise False.
    """
    return bool(file) and os.path.isfile(file) and os.path.getsize(file) > 0


def get_star_version() -> str:
    """
    Gets the version of the STAR binary.

    Returns:
        The version of STAR or "Not available" if not found.
    """
    try:
        proc = subprocess.Popen(
            ["STAR", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False,
            close_fds=True
        )
        stdout, _ = proc.communicate()
        return stdout.decode().strip()
    except Exception:
        return "Not available"


def get_taggd_count_version() -> str:
    """
    Gets the version of the Taggd binary.

    Returns:
        The version of Taggd or "Not available" if not found.
    """
    try:
        proc = subprocess.Popen(
            ["pip", "show", "taggd"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False,
            close_fds=True
        )
        stdout, _ = proc.communicate()
        for line in stdout.decode().splitlines():
            if "Version:" in line:
                return line.split()[-1]
    except Exception:
        pass
    return "Not available"


def get_htseq_count_version() -> str:
    """
    Gets the version of the HTSeqCount binary.

    Returns:
        The version of HTSeqCount or "Not available" if not found.
    """
    try:
        proc = subprocess.Popen(
            ["pip", "show", "htseq"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False,
            close_fds=True
        )
        stdout, _ = proc.communicate()
        for line in stdout.decode().splitlines():
            if "Version:" in line:
                return line.split()[-1]
    except Exception:
        pass
    return "Not available"