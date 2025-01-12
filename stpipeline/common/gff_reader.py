"""
A simple module to parse a GFF/GTF file and query it
"""
import gzip
import re
from typing import Generator, Optional, Union, Dict, Any, List

# Code snipped from:
# https://gist.github.com/slowkow/8101481

GTF_HEADER = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame"]
R_SEMICOLON = re.compile(r"\s*;\s*")
R_COMMA = re.compile(r"\s*,\s*")
R_KEYVALUE = re.compile(r"(\s+|\s*=\s*)")


def gff_lines(filename: str) -> Generator[Dict[str, Any], None, None]:
    """
    Opens an optionally gzipped GTF/GFF file and generates a dictionary for each line.

    Args:
        filename: Path to the GTF/GFF file. The file can be gzipped.

    Yields:
       Parsed fields from each GTF/GFF line.
    """
    fn_open = gzip.open if filename.endswith(".gz") else open
    # 'rt' ensures reading text from gzipped file
    with fn_open(filename, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            else:
                yield gff_parse(line)


def _get_value(value: Optional[str]) -> Optional[Union[str, List[Optional[str]]]]:
    """
    Processes a value from the GTF/GFF file, stripping quotes and handling lists.

    Args:
        value: The value to process.

    Returns:
        Processed value, or None if the value is equivalent to null.
    """
    if not value:
        return None
    # Strip double and single quotes
    value = value.strip("\"'")
    if "," in value:
        # Return a list if the value contains commas
        value = re.split(R_COMMA, value)  # type: ignore
    # Handle equivalent-to-null values
    elif value in ["", ".", "NA"]:
        return None
    return value


def gff_parse(line: str) -> Dict[str, Any]:
    """
    Parses a single GTF/GFF line and returns a dictionary of its fields.

    Args:
        line: A single line from a GTF/GFF file.

    Returns:
        Parsed fields from the line as key-value pairs.
    """
    result: Dict[str, Any] = {}
    fields = line.rstrip().split("\t")

    # Parse standard fields
    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # Parse INFO field
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]
    for i, info in enumerate(infos, 1):
        try:
            # Parse "key=value"
            key, _, value = re.split(R_KEYVALUE, info, maxsplit=1)
        except ValueError:
            # Use INFO1, INFO2, etc. for unnamed values
            key = f"INFO{i}"
            value = info
        # Ignore fields with no value
        if value:
            result[key] = _get_value(value)

    return result
