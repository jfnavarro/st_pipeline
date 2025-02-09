"""
This file contains some general utilities
"""

import os
import shutil
import subprocess
import threading
from datetime import datetime
from typing import IO, Any, Optional


def which_program(program: str) -> bool:
    """
    Check if a program is installed and available in the system's PATH.

    Args:
        program: The name of the program to check.

    Returns:
        True if the program is found and executable, False otherwise.
    """
    return shutil.which(program) is not None


class TimeStamper:
    """
    Thread-safe time stamper to generate unique numeric timestamps.
    """

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.prev: Optional[datetime] = None
        self.count: int = 0

    def get_timestamp(self) -> datetime:
        """
        Generates a unique numeric timestamp.

        Returns:
            A unique timestamp as a datetime object.
        """
        with self.lock:
            ts = datetime.now()
            if ts == self.prev:
                ts += ".%04d" % self.count  # type: ignore
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


def safe_open_file(filename: str, mode: str) -> IO[Any]:
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
    if mode not in ["w", "r"]:
        raise IOError(f"Error: Invalid mode: {mode}")
    if "w" in mode:
        safe_remove(filename)
    elif "r" in mode and not os.path.isfile(filename):
        raise IOError(f"Error: File does not exist: {filename}")

    return open(filename, mode)


def file_ok(file: Optional[str]) -> bool:
    """
    Checks if a file exists and is not empty.

    Args:
        file: Path to the file.

    Returns:
        True if the file exists and is not empty, otherwise False.
    """
    return file is not None and os.path.isfile(file) and os.path.getsize(file) > 0


def get_star_version() -> str:
    """
    Gets the version of the STAR binary.

    Returns:
        The version of STAR or "Not available" if not found.
    """
    try:
        proc = subprocess.Popen(
            ["STAR", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, close_fds=True
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
            ["pip", "show", "taggd"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, close_fds=True
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
            ["pip", "show", "htseq"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, close_fds=True
        )
        stdout, _ = proc.communicate()
        for line in stdout.decode().splitlines():
            if "Version:" in line:
                return line.split()[-1]
    except Exception:
        pass
    return "Not available"
