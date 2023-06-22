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

from collections import Counter

from .readers.main import parse_reads
from .readers.read_file_info import ReadFileInfo
from .stats import LibraryIndependentStats


def run_library_independent_counting(
    read_file_info: ReadFileInfo,
    sample: str | None,
    cpus: int,
    reference_fp: str | None,
    min_length: int = 0
) -> tuple[int, LibraryIndependentStats, Counter[str]]:
    unique_count: int
    counts: Counter[str]

    unique_count, stats, counts = parse_reads(  # type: ignore
        read_file_info,
        sample=sample,
        cpus=cpus,
        reference=reference_fp,
        exclude_by_len=min_length,
        exclude_qcfail=False)

    return unique_count, stats, counts
