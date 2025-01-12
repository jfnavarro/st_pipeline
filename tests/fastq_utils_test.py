#! /usr/bin/env python
""" 
Unit-test the package fastq_utils
"""
import pytest
from stpipeline.common.fastq_utils import (
    remove_adaptor,
    quality_trim_index,
    trim_quality,
    check_umi_template,
    has_sufficient_content
)

# Test for remove_adaptor
def test_remove_adaptor():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFFFFFFFFFFF"
    adaptor = "TAGCTT"
    trimmed_seq, trimmed_qual = remove_adaptor(sequence, quality, adaptor, missmatches=0)
    assert trimmed_seq == "AGCT"
    assert trimmed_qual == "FFFF"

    # Test no adaptor found
    trimmed_seq, trimmed_qual = remove_adaptor(sequence, quality, "GATTACA")
    assert trimmed_seq == sequence
    assert trimmed_qual == quality

    # Test with mismatches
    trimmed_seq, trimmed_qual = remove_adaptor(sequence, quality, "TAGCTT", missmatches=1)
    assert trimmed_seq == "AGCT"
    assert trimmed_qual == "FFFF"

def test_quality_trim_index_basic():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFFFFFFFFFFF"  # ASCII 'F' -> Phred score 40
    cutoff = 20
    result = quality_trim_index(sequence, quality, cutoff)
    assert result == len(sequence)  # No trimming, all bases are high quality


def test_quality_trim_index_trimming():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFF!!!!!!!!!"  # Phred scores: 'F' (40), '!' (0)
    cutoff = 20
    result = quality_trim_index(sequence, quality, cutoff)
    assert result == 6  # Trims after the first 6 high-quality bases


def test_quality_trim_index_low_quality_g():
    sequence = "AGCTTAGCTTGGA"
    quality = "FFFFFF!!!!!!"  # Phred scores: 'F' (40), '!' (0)
    cutoff = 20
    result = quality_trim_index(sequence, quality, cutoff)
    assert result == 6  # Trims after the first 6 high-quality bases


def test_trim_quality_basic():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFFFFFFFFFFF"  # All high quality
    min_qual = 20
    min_length = 10
    trimmed_seq, trimmed_qual = trim_quality(sequence, quality, min_qual, min_length)
    assert trimmed_seq == "AGCTTAGCTTAGCTA"
    assert trimmed_qual == "FFFFFFFFFFFFFFF"


def test_trim_quality_trimming():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFF!!!!!!!!!"  # Low-quality bases at the end
    min_qual = 20
    min_length = 5
    trimmed_seq, trimmed_qual = trim_quality(sequence, quality, min_qual, min_length)
    assert trimmed_seq == "AGCTTA"
    assert trimmed_qual == "FFFFFF"


def test_trim_quality_below_min_length():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFF!!!!!!!!!"  # Low-quality bases at the end
    min_qual = 20
    min_length = 10
    trimmed_seq, trimmed_qual = trim_quality(sequence, quality, min_qual, min_length)
    assert trimmed_seq is None
    assert trimmed_qual is None


def test_trim_quality_low_quality_g():
    sequence = "AGCTTAGCTTGGA"
    quality = "FFFFFF!!!!!!"  # Phred scores: 'F' (40), '!' (0)
    min_qual = 20
    min_length = 5
    trimmed_seq, trimmed_qual = trim_quality(sequence, quality, min_qual, min_length)
    assert trimmed_seq == "AGCTTA"
    assert trimmed_qual == "FFFFFF"

def test_trim_quality_short():
    min_qual = 20
    min_length = 10

    # Test with sequence shorter than min_length
    trimmed_seq, trimmed_qual = trim_quality("AGCTT", "FFFFF", min_qual, min_length)
    assert trimmed_seq is None
    assert trimmed_qual is None

# Test for check_umi_template
def test_check_umi_template():
    umi = "ACGT1234"
    template = r"[ACGT]{4}\d{4}"
    assert check_umi_template(umi, template) is True

    umi = "ACGT12"
    assert check_umi_template(umi, template) is False

# Test for has_sufficient_content
def test_has_sufficient_content():
    sequence = "ATATGGCCATAT"
    chars_to_count = "AT"
    threshold = 50.0
    assert has_sufficient_content(sequence, chars_to_count, threshold) is True

    threshold = 90.0
    assert has_sufficient_content(sequence, chars_to_count, threshold) is False

    # Test empty sequence
    with pytest.raises(ValueError):
        has_sufficient_content("", chars_to_count, threshold)

    # Test empty chars_to_count
    with pytest.raises(ValueError):
        has_sufficient_content(sequence, "", threshold)
