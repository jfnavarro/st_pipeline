"""
Module for distance metrics.
"""
import distance  # type: ignore


def hamming_distance(a: str, b: str) -> int:
    """
    Calculates the Hamming distance between two strings using the `distance` library.

    Args:
        a: First string.
        b: Second string.

    Returns:
        The Hamming distance between the two strings.

    Raises:
        ValueError: If the strings are of unequal length.
    """
    if len(a) != len(b):
        raise ValueError("Strings must be of equal length to calculate Hamming distance")
    return distance.hamming(a, b)  # type: ignore
