"""
This module defines the GeneBuffer class and the parse_unique_events function.
"""

import logging
import pysam
import operator
from .gff_reader import gff_lines
from .dataset import Transcript
from typing import Dict, Generator, Optional, Tuple, List

logger = logging.getLogger("STPipeline")

class GeneBuffer:
    """
    This object defines a buffer by holding a dictionary 
    of genes, spot coordinates, and transcripts.
    It assumes the transcripts are added in a coordinate-ordered fashion.

    Attributes:
        buffer: Dictionary storing gene data.
        last_position: Last genomic position processed.
        last_chromosome: Last chromosome processed.
        gene_end_coordinates: Dictionary mapping gene IDs to their end coordinates and chromosomes.
    """

    def __init__(self, gff_filename: Optional[str]):
        """
        Initializes the GeneBuffer object and computes gene end coordinates from a GFF file.

        Args:
            gff_filename: Path to the GFF file containing gene annotations.
        """
        self.buffer = {}
        self.last_position = 0
        self.last_chromosome = 'chrom'
        self.gene_end_coordinates = {}
        if gff_filename:
            self.__compute_gene_end_coordinates(gff_filename)

    def __compute_gene_end_coordinates(self, gff_filename: str) -> None:
        """
        Reads the end coordinates and chromosomes of all genes present in the GFF file
        and saves them in a dictionary with the gene ID as the key.

        Args:
            gff_filename: Path to the GFF file.
        """
        logger.debug(f"Parsing GFF file {gff_filename} to compute gene end coordinates.")

        gene_end_coordinates = {}
        for line in gff_lines(gff_filename):
            seqname = line["seqname"]
            end = int(line["end"])
            gene_id = line.get("gene_id", None)
            if not gene_id:
                msg = f"The gene_id attribute is missing in the annotation file ({gff_filename})."
                loggger.error(msg)
                raise ValueError(msg)

            if gene_id[0] == '"' and gene_id[-1] == '"': 
                gene_id = gene_id[1:-1]

            if gene_id in gene_end_coordinates:
                if end > gene_end_coordinates[gene_id][1]:
                    gene_end_coordinates[gene_id] = (seqname, end)
            else:
                gene_end_coordinates[gene_id] = (seqname, end)

        gene_end_coordinates['__no_feature'] = (None, -1)
        self.gene_end_coordinates = gene_end_coordinates

    def get_gene_end_position(self, gene: str) -> Tuple[Optional[str], int]:
        """
        Returns the genomic end coordinate and chromosome of the given gene.

        Args:
            gene (str): Gene ID.

        Returns:
            Tuple[Optional[str], int]: Chromosome and end coordinate of the gene.

        Raises:
            ValueError: If the gene is not found in the annotation file or is ambiguous.
        """
        try:
            return self.gene_end_coordinates[gene]
        except KeyError:
            if '__ambiguous[' in gene:
                ambiguous_genes = gene[gene.index('[') + 1:gene.index(']')].split('+')
                try:
                    return max(
                        [self.gene_end_coordinates[amb_gene] for amb_gene in ambiguous_genes],
                        key=operator.itemgetter(1)
                    )
                except KeyError:
                    raise ValueError(f"Ambiguous gene {gene} not found in annotation file.")
            raise ValueError(f"Gene {gene} not found in annotation file.")

    def add_transcript(self, gene: str, spot_coordinates: Tuple[int, int], transcript: Transcript, position: int) -> None:
        """
        Adds a transcript to the gene buffer.

        Args:
            gene (str): Gene name.
            spot_coordinates (Tuple[int, int]): Spot coordinates (x, y).
            transcript (Transcript): Transcript information.
            position (int): Transcript's left-most genomic coordinate.
        """
        self.last_position = position
        self.last_chromosome = transcript.chrom

        self.buffer.setdefault(gene, {}).setdefault(spot_coordinates, []).append(transcript)

    def check_and_clear_buffer(self, empty: bool = False) -> Generator[Tuple[str, Dict], None, None]:
        """
        Checks and clears the buffer, yielding genes that are outside the current chromosome or position.

        Args:
            empty (bool): If True, forces clearing the buffer.

        Yields:
            Tuple[str, Dict]: Gene name and its buffer content.
        """
        for gene in list(self.buffer.keys()):
            if gene == '__no_feature' and not empty:
                continue

            chrom, end_position = self.get_gene_end_position(gene)

            if empty or self.last_position > end_position or self.last_chromosome != chrom:
                yield gene, self.buffer.pop(gene)

def parse_unique_events(input_file: str, gff_filename: Optional[str] = None) -> Generator[Tuple[str, Dict], None, None]:
    """
    Parses transcripts from a coordinate-sorted BAM file and organizes them by gene and spot coordinates.

    Args:
        input_file (str): Path to the input BAM file containing annotated records.
        gff_filename (Optional[str]): Path to the GFF file containing gene coordinates.

    Yields:
        Tuple[str, Dict]: Gene name and a dictionary mapping spot coordinates to transcripts.
    """
    genes_buffer = GeneBuffer(gff_filename) if gff_filename else None
    genes_dict: Dict[str, Dict[Tuple[int, int], List[Transcript]]] = {}

    sam_file = pysam.AlignmentFile(input_file, "rb")

    for rec in sam_file.fetch(until_eof=True):
        clear_name = rec.query_name
        mapping_quality = rec.mapping_quality
        start = rec.reference_start - rec.query_alignment_start
        end = rec.reference_end + (rec.query_length - rec.query_alignment_end)
        chrom = sam_file.get_reference_name(rec.reference_id)
        strand = "-" if rec.is_reverse else "+"

        if strand == "-":
            start, end = end, start

        x, y, gene, umi = -1, -1, 'None', 'None'
        for k, v in rec.tags:
            if k == "B1":
                x = int(v)
            elif k == "B2":
                y = int(v)
            elif k == "XF":
                gene = str(v)
            elif k == "B3":
                umi = str(v)

        transcript = Transcript(chrom, start, end, clear_name, mapping_quality, strand, umi)

        if genes_buffer:
            genes_buffer.add_transcript(gene, (x, y), transcript, rec.reference_start)
            for g, t in genes_buffer.check_and_clear_buffer():
                yield g, t
        else:
            genes_dict.setdefault(gene, {}).setdefault((x, y), []).append(transcript)

    sam_file.close()

    if genes_buffer:
        for g, t in genes_buffer.check_and_clear_buffer(True):
            yield g, t
    else:
        yield from genes_dict.items()