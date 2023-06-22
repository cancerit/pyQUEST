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

from pyquest.readers.read_info import ReadInfo
from pyquest.readers.read_qc import is_dna, is_masked, is_ambiguous
from pyquest.stats import LibraryIndependentStats


@pytest.mark.parametrize('seq,exp', [
    ('GACGTAGATCTGAT', True),
    ('AT', True),
    ('AGGTXAT', False),
    ('', False)
])
def test_is_dna(seq: str, exp: bool):
    assert is_dna(seq) is exp


@pytest.mark.parametrize('seq,exp', [
    ('GACGTagaATCTGAT', True),
    ('at', True),
    ('AGGTXAT', False),
    ('xT', False),
    ('', False)
])
def test_is_masked(seq: str, exp: bool):
    assert is_masked(seq) is exp


@pytest.mark.parametrize('seq,exp', [
    ('GACGTagaATCTGAT', False),
    ('ATNGT', True),
    ('AGGTXAT', False),
    ('GnnTA', True),
    ('', False)
])
def test_is_ambiguous(seq: str, exp: bool):
    assert is_ambiguous(seq) is exp


@pytest.mark.parametrize('seq,short,masked,ambiguous', [
    ('ACGgtTAC', False, True, False),
    ('ACGgtNNN', False, True, True),
    ('AC', True, False, False),
    ('ac', True, True, False),
    ('an', True, True, True),
    ('ACNgtTAC', False, True, True),
    ('ACGgnTAC', False, True, True)
])
def test_eval_read_update_stats(seq: str, short: bool, masked: bool, ambiguous: bool):
    stats = LibraryIndependentStats.empty('A')
    min_length = 3
    qc_fail_discard = False

    # Fill read info
    read_info = ReadInfo(
        sequence=seq,
        is_masked=masked,
        is_ambiguous=ambiguous,
        is_short=len(seq) < min_length,
        is_qc_fail=False,
        is_flagged=False)

    # Check read info
    assert read_info.to_discard(qc_fail_discard) == short or masked or ambiguous

    # Update stats
    assert stats.eval_read_info(qc_fail_discard, read_info) != short or masked or ambiguous

    # Check stats
    assert stats.length_excluded_reads == (1 if short else 0)
    assert stats.masked_reads == (1 if masked else 0)
    assert stats.ambiguous_nt_reads == (1 if ambiguous else 0)
    assert stats.total_invalid_reads == (1 if masked or ambiguous else 0)
    assert stats.total_excluded_reads == (1 if masked or ambiguous or short else 0)
