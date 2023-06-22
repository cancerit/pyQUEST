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

# Adapted from pyCROQUET (readparser module)

from collections import Counter
import logging
import re
from typing import TextIO

from pyquest.readers.read_info import ReadInfo
from pyquest.readers.read_qc import is_ambiguous, is_dna, is_masked

from ..errors import InputReadError, InvalidFASTQError
from ..stats import LibraryIndependentStats
from .constants import LOAD_INFO_THRESHOLD


ILLUMINA_SINGLE_FASTQ_HEADER_PATTERN = re.compile(r"^@([^\s/]+)$")
ILLUMINA_FASTQ_HEADER_PATTERN = re.compile(r"^@(\S+)/([12])$")
CASAVA_FASTQ_HEADER_PATTERN = re.compile(r"^@(\S+)\s([012]):([YN])+:[\d+]+:\S+$")


OFFSET_ILLUMINA = 64
OFFSET_CASAVA = 33


def _parse_fq_header(header: str):
    qc_fail = False
    phred_offset = OFFSET_ILLUMINA
    if (match := ILLUMINA_SINGLE_FASTQ_HEADER_PATTERN.match(header)) is not None:
        # illumina format, unparied
        groups = match.groups()
        name = groups[0]
        pair_member = None
    elif (match := ILLUMINA_FASTQ_HEADER_PATTERN.match(header)) is not None:
        # illumina format, paired or single end wit read identifier
        groups = match.groups()
        name = groups[0]
        pair_member = int(groups[1])
    elif (match := CASAVA_FASTQ_HEADER_PATTERN.match(header)) is not None:
        # casava1.8+ format
        phred_offset = OFFSET_CASAVA
        groups = match.groups()
        name = groups[0]
        pair_member = int(groups[1])
        if groups[2] == "Y":
            qc_fail = True
    else:
        raise InvalidFASTQError(f"Unsupported FastQ header format: {header}")
    return name, pair_member, qc_fail, phred_offset


def get_fastq_read_info(min_length: int, seq: str, qc_fail: bool) -> ReadInfo:
    if not is_dna(seq):
        raise InputReadError(f"Invalid sequence: '{seq}'!")

    return ReadInfo(
        sequence=seq,
        is_flagged=False,
        is_qc_fail=qc_fail,
        is_ambiguous=is_ambiguous(seq),
        is_masked=is_masked(seq),
        is_short=(len(seq) < min_length))


def parse_fastq(
    ifh: TextIO,
    sample: str,
    exclude_qcfail: bool = False,
    exclude_by_len: int | None = None
) -> tuple[LibraryIndependentStats, Counter[str]]:
    """
    Closes received file handle
    Only used for the reads seq minimization process
    """

    min_length: int = exclude_by_len if exclude_by_len is not None else 0
    stats = LibraryIndependentStats.empty(sample)
    reads: Counter[str] = Counter()

    def eval_read(seq: str, qc_fail: bool) -> None:
        read_info: ReadInfo = get_fastq_read_info(min_length, seq, qc_fail)
        if stats.eval_read_info(exclude_qcfail, read_info):
            reads[read_info.sequence] += 1

    seq: str
    qc_fail: bool

    header = ifh.readline()
    while header:

        seq = ifh.readline().strip()
        _ = ifh.readline()  # throw away separator
        _ = ifh.readline()  # throw away qual at this point
        _, _, qc_fail, _ = _parse_fq_header(header.strip())  # for validation only here
        # looks odd, but best way to handle end of file
        header = ifh.readline()

        eval_read(seq, qc_fail)

        if stats.total_reads % LOAD_INFO_THRESHOLD == 0:  # pragma: no cover
            logging.debug(f"Parsed {stats.total_reads} reads, {len(reads)} were unique...")

    return stats, reads
