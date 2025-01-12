#! /usr/bin/env python
"""
Unit-test the package utils
"""
import pytest
from unittest.mock import patch
import os
from stpipeline.common.utils import (
    which_program,
    TimeStamper,
    safe_remove,
    safe_open_file,
    file_ok,
    get_star_version,
    get_taggd_count_version,
    get_htseq_count_version,
)


@pytest.fixture
def temp_file(tmp_path):
    temp = tmp_path / "temp_file.txt"
    temp.write_text("Temporary file content.")
    return str(temp)


def test_which_program():
    program = "python"
    result = which_program(program)
    assert result is True


def test_which_program_not_found():
    program = "nonexistent_program"
    result = which_program(program)
    assert result is False


def test_timestamper():
    stamper = TimeStamper()
    ts1 = stamper.get_timestamp()
    ts2 = stamper.get_timestamp()
    assert ts1 != ts2


def test_safe_remove(temp_file):
    assert os.path.exists(temp_file)
    safe_remove(temp_file)
    assert not os.path.exists(temp_file)


def test_safe_remove_nonexistent():
    non_existent_file = "non_existent_file.txt"
    safe_remove(non_existent_file)
    assert not os.path.exists(non_existent_file)


def test_safe_open_file_read(temp_file):
    with safe_open_file(temp_file, "r") as f:
        content = f.read()
    assert content == "Temporary file content."


def test_safe_open_file_write(tmp_path):
    file_path = tmp_path / "test_file.txt"
    with safe_open_file(str(file_path), "w") as f:
        f.write("Test content.")
    assert file_path.exists()
    with open(file_path, "r") as f:
        assert f.read() == "Test content."


def test_safe_open_file_invalid_mode():
    with pytest.raises(IOError):
        safe_open_file("invalid_file.txt", "x")


def test_file_ok(temp_file):
    assert file_ok(temp_file)


def test_file_ok_empty_file(tmp_path):
    empty_file = tmp_path / "empty.txt"
    empty_file.touch()
    assert not file_ok(str(empty_file))


def test_get_star_version():
    with patch("subprocess.Popen") as mock_popen:
        mock_popen.return_value.communicate.return_value = (b"STAR_2.7.9a", b"")
        mock_popen.return_value.returncode = 0
        version = get_star_version()
        assert version == "STAR_2.7.9a"


def test_get_star_version_not_found():
    with patch("subprocess.Popen", side_effect=FileNotFoundError):
        version = get_star_version()
        assert version == "Not available"


def test_get_taggd_count_version():
    with patch("subprocess.Popen") as mock_popen:
        mock_popen.return_value.communicate.return_value = (b"Name: taggd\nVersion: 1.0.0", b"")
        mock_popen.return_value.returncode = 0
        version = get_taggd_count_version()
        assert version == "1.0.0"


def test_get_taggd_count_version_not_found():
    with patch("subprocess.Popen", side_effect=FileNotFoundError):
        version = get_taggd_count_version()
        assert version == "Not available"


def test_get_htseq_count_version():
    with patch("subprocess.Popen") as mock_popen:
        mock_popen.return_value.communicate.return_value = (b"Name: htseq\nVersion: 0.11.3", b"")
        mock_popen.return_value.returncode = 0
        version = get_htseq_count_version()
        assert version == "0.11.3"


def test_get_htseq_count_version_not_found():
    with patch("subprocess.Popen", side_effect=FileNotFoundError):
        version = get_htseq_count_version()
        assert version == "Not available"
