""" 
This module contains some specific functions for
to parse and modify FASTQ files
"""
import re
import regex
from typing import Tuple, Union

def remove_adaptor(
    sequence: str,
    quality: str,
    adaptor: str,
    missmatches: int = 2
) -> Tuple[str, str]:
    """
    Trims a given adaptor sequence from a FASTQ read if it is found.

    Args:
        sequence: The sequence of the read.
        quality: The quality string corresponding to the read.
        adaptor: The adaptor sequence to search for.
        missmatches: The allowed number of mismatches when searching for the adaptor. Defaults to 2.

    Returns:
       A tuple containing the trimmed sequence and quality strings.

    Raises:
        ValueError: If the input sequence and quality lengths do not match.
    """
    if len(sequence) < len(adaptor) or len(sequence) != len(quality):
        return sequence, quality

    try:
        if missmatches == 0:
            pos = sequence.find(adaptor)
        else:
            candidates = regex.findall(
                rf'(?:{adaptor}){{s<={missmatches}}}',
                sequence,
                overlapped=False
            )
            if candidates:
                local_seq = candidates[0]
                local_pos = 0
                if adaptor[0] != local_seq[0]:
                    local_pos = local_seq.find(adaptor[0])
                pos = sequence.find(local_seq[local_pos:])
            else:
                pos = -1

        if pos != -1:
            return sequence[:pos], quality[:pos]
        else:
            return sequence, quality
    except Exception as e:
        raise RuntimeError(f"Failed to trim adaptor: {e}")

def quality_trim_index(bases: str, qualities: str, cutoff: int, base: int = 33) -> int:
    """
    Find the position at which to trim a low-quality end from a nucleotide sequence.

    Qualities are assumed to be ASCII-encoded as chr(qual + base).

    This algorithm is derived from BWA's 'bwa_trim_read':
    - Subtract the cutoff value from all qualities.
    - Compute partial sums from all indices to the end of the sequence.
    - Trim sequence at the index at which the sum is minimal.

    Args:
        bases: Nucleotide sequence.
        qualities: ASCII-encoded quality scores.
        cutoff: Quality cutoff value.
        base: Base value for ASCII encoding. Defaults to 33.

    Returns:
        Index position to trim the sequence.

    Note:
        This function handles Illumina NextSeq data specifically by treating high-quality 'G' bases
        at the end of reads as having a quality of (cutoff - 1).

    References:
        CutAdapt (https://github.com/marcelm/cutadapt/)
    """
    s = 0
    max_qual = 0
    max_i = len(qualities)
    for i in reversed(range(max_i)):
        q = ord(qualities[i]) - base
        if bases[i] == 'G':
            q = cutoff - 1
        s += cutoff - q
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            max_i = i
    return max_i


def trim_quality(
    sequence: str,
    quality: str,
    min_qual: int = 20,
    min_length: int = 30,
    phred: int = 33
) -> Tuple[Union[str, None], Union[str, None]]:
    """
    Quality trims a FASTQ read using a BWA-like approach.

    The function trims a nucleotide sequence and its quality scores based on a minimum quality threshold. 
    If the trimmed sequence is shorter than a minimum length, it returns None.

    Args:
        sequence: Nucleotide sequence of the read.
        quality: Quality scores of the read, ASCII-encoded.
        min_qual: Quality threshold to trim. Defaults to 20.
        min_length: Minimum valid length for a read after trimming. Defaults to 30.
        phred: Phred encoding format for quality scores (33 or 64). Defaults to 33.

    Returns:
        A tuple containing the trimmed sequence and quality scores, 
        or (None, None) if trimming results in a sequence shorter than `min_length`.
    """
    if len(sequence) < min_length:
        return None, None

    # Get the position at which to trim (number of bases to trim)
    cut_index = quality_trim_index(sequence, quality, min_qual, phred)

    # Check if the trimmed sequence would have at least the minimum length
    if (cut_index + 1) >= min_length:
        new_seq = sequence[:cut_index]
        new_qual = quality[:cut_index]
        return new_seq, new_qual
    else:
        return None, None


def check_umi_template(umi, template):
    """
    Checks that the UMI (molecular barcode) given as input complies
    with the pattern given in template.
    Returns True if the UMI complies
    :param umi: a molecular barcode
    :param template: a reg-based template with the same
                    distance of the UMI that should tell how the UMI should be formed
    :type umi: str
    :type template: str
    :return: True if the given molecular barcode fits the pattern given
    """
    p = re.compile(template)
    return p.match(umi) is not None
