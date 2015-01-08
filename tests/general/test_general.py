import os

import pytest

import pipelines as ps

class TestPicardIndex:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        index = 'test.bam.bai'
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        temp_dir = 'path/to/temp/dir'
        lines = ps.general._picard_index(in_bam, index, picard_memory,
                                         picard_path, temp_dir)
