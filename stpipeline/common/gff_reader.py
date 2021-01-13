"""
A simple module to parse a GFF/GTF file and query it
"""
import gzip
import re

# Code snipped from:
# https://gist.github.com/slowkow/8101481

GTF_HEADER = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA = re.compile(r'\s*,\s*')
R_KEYVALUE = re.compile(r'(\s+|\s*=\s*)')


def gff_lines(filename):
    """
    Opens an optionally gzipped GTF/GFF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield gff_parse(line)


def gff_parse(line):
    """
    Parses a single GTF/GFF line and returns a dict.
    """
    result = {}
    fields = line.rstrip().split('\t')
    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])
    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]
    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)
    return result


def _get_value(value):
    if not value:
        return None
    # Strip double and single quotes.
    value = value.strip('"\'')
    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None
    return value
