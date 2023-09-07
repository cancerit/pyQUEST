# pyQUEST
#
# Copyright (C) 2023 Genome Research Ltd.
#
# Author: Luca Barbon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import pytest
import pysam

from pyquest.readers.fq import get_fastq_read_info
from pyquest.readers.hts import get_hts_read_info
from pyquest.errors import InputReadError

@pytest.mark.parametrize('seq,exp_discard', [
    ('GACGTAGATCTGAT', False),
    ('', True)
])
def test_get_fq_read_info(seq:str, exp_discard: bool):
    read_info = get_fastq_read_info(min_length=0, seq=seq, qc_fail=False)

    # Test read parsed and marked for discard correctly
    assert read_info.to_discard(discard_qc=False) == exp_discard


def test_get_fq_read_info_error():
    seq = "AGG$$$44"

    # Test read parsing throws error correctly
    pytest.raises(InputReadError, get_fastq_read_info, 0, seq, False)

def _make_test_read_from_seq(seq:str) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_sequence = seq
    return read

@pytest.mark.parametrize('seq,exp_discard', [
    ('GACGTAGATCTGAT', False),
    ('', True)
])
def test_get_hts_read_info(seq:str, exp_discard: bool):
    read = _make_test_read_from_seq(seq)

    read_info = get_hts_read_info(min_length=0, skip_read_flag=False, read=read)

    # Test read parsed and marked for discard correctly
    assert read_info.to_discard(discard_qc=False) == exp_discard


def test_get_hts_read_info_error():
    seq = "AGG$$$44"
    read = _make_test_read_from_seq(seq)

    # Test read parsing throws error correctly
    pytest.raises(InputReadError, get_hts_read_info, 0, False, read)