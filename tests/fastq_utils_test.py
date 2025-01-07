#! /usr/bin/env python
""" 
Unit-test the package fastq_utils
"""
import pytest
from your_module_name import (
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
    adaptor = "TAGCTA"
    trimmed_seq, trimmed_qual = remove_adaptor(sequence, quality, adaptor)
    assert trimmed_seq == "AGCTT"
    assert trimmed_qual == "FFFFF"

    # Test no adaptor found
    trimmed_seq, trimmed_qual = remove_adaptor(sequence, quality, "GATTACA")
    assert trimmed_seq == sequence
    assert trimmed_qual == quality

    # Test with mismatches
    trimmed_seq, trimmed_qual = remove_adaptor(sequence, quality, "TAGCTT", missmatches=1)
    assert trimmed_seq == "AGCT"
    assert trimmed_qual == "FFFF"

# Test for quality_trim_index
def test_quality_trim_index():
    bases = "AGCTTAGCTTAGCTA"
    qualities = "FFFFFFFFFFFFFFF"
    cutoff = 20
    index = quality_trim_index(bases, qualities, cutoff)
    assert index == len(bases)  # No trimming needed for high-quality bases

    qualities = "FFFFF####FFFFFF"
    index = quality_trim_index(bases, qualities, cutoff)
    assert index == 5  # Trimming occurs at the first low-quality base

# Test for trim_quality
def test_trim_quality():
    sequence = "AGCTTAGCTTAGCTA"
    quality = "FFFFFFFFFFFFFFF"
    min_qual = 20
    min_length = 10

    trimmed_seq, trimmed_qual = trim_quality(sequence, quality, min_qual, min_length)
    assert trimmed_seq == "AGCTTAGCTT"
    assert trimmed_qual == "FFFFFFFFFF"

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
