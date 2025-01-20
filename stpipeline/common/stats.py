"""
This shared object is used to collect
different statistics and QA parameters for
the pipeline run.
"""

import json
from dataclasses import asdict, dataclass, field
from typing import List


@dataclass
class Stats:
    """
    Stats collects information and statistics for use in the ST pipeline.
    """

    input_reads_forward: int = 0
    input_reads_reverse: int = 0
    reads_after_trimming_forward: int = 0
    reads_after_trimming_reverse: int = 0
    reads_after_rRNA_trimming: int = 0
    reads_after_mapping: int = 0
    reads_after_annotation: int = 0
    reads_after_demultiplexing: int = 0
    reads_after_duplicates_removal: int = 0
    genes_found: int = 0
    duplicates_found: int = 0
    pipeline_version: str = "-"
    mapper_tool: str = "-"
    annotation_tool: str = "-"
    demultiplex_tool: str = "-"
    input_parameters: List[str] = field(default_factory=list)
    max_genes_feature: int = 0
    min_genes_feature: int = 0
    max_reads_feature: int = 0
    min_reads_feature: int = 0
    average_gene_feature: float = 0.0
    average_reads_feature: float = 0.0

    def __str__(self) -> str:
        """
        Returns a string representation of the Stats object.

        Returns:
            A formatted string of all stats attributes.
        """
        return "\n".join(f"{field_name}: {getattr(self, field_name)}" for field_name in self.__dataclass_fields__)

    def write_json(self, filename: str) -> None:
        """
        Writes the stats to a JSON file.

        Args:
            filename: The path to the JSON file to write.
        """
        with open(filename, "w") as file:
            json.dump(asdict(self), file, indent=2, separators=(",", ": "))

    @classmethod
    def from_json(cls, filename: str) -> "Stats":
        """
        Creates a Stats object from a JSON file.

        Args:
            filename: The path to the JSON file to read.

        Returns:
            A Stats object populated with data from the JSON file.
        """
        with open(filename, "r") as file:
            data = json.load(file)
        return cls(**data)
