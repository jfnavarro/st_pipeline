from dataclasses import dataclass


@dataclass
class Transcript:
    """
    Represents a transcript with associated genomic and metadata information.

    Attributes:
        chrom (str): Chromosome name.
        start (int): Start position of the transcript.
        end (int): End position of the transcript.
        clear_name (str): Clear name of the transcript.
        mapping_quality (int): Mapping quality score.
        strand (str): Strand information ('+' or '-').
        umi (str): Unique Molecular Identifier (UMI).
    """

    chrom: str
    start: int
    end: int
    clear_name: str
    mapping_quality: int
    strand: str
    umi: str
