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
from contextlib import contextmanager
import logging
from typing import Final, Generator

import pysam
from pysam.libcalignedsegment import SAM_FLAGS

from ..stats import LibraryIndependentStats
from .constants import LOAD_INFO_THRESHOLD
from ..errors import InputReadError, InvalidHTSError, MissingMetadataError
from .read_file_info import ReadFileInfo
from .read_info import ReadInfo
from .read_qc import is_ambiguous, is_dna


MAX_HTS_CPUS = 4
SKIP_READ_FLAG: Final[int] = (
    SAM_FLAGS.FSECONDARY |
    SAM_FLAGS.FSUPPLEMENTARY
)


def _cap_hts_cpus(cpus: int) -> int:
    return cpus if cpus < MAX_HTS_CPUS else MAX_HTS_CPUS


def _get_hts_sample_name(sam: pysam.AlignmentFile, sample: str | None = None) -> str | None:
    if sample is not None:
        return sample

    header = sam.header.as_dict()  # type: ignore

    if "RG" not in header:
        return sample

    for rg in header["RG"]:
        if "SM" not in rg:
            continue
        if sample is None:
            sample = rg["SM"]
        elif sample != rg["SM"]:
            raise InvalidHTSError("Multiple different sample names found in header")

    return sample


@contextmanager
def hts_reader(seq_file: str, mode, cpus: int, reference: str | None = None) -> Generator[pysam.AlignmentFile, None, None]:
    verbosity = pysam.set_verbosity(0)
    try:
        sam = pysam.AlignmentFile(
            seq_file,
            mode=mode,
            reference_filename=reference,
            require_index=False,
            threads=_cap_hts_cpus(cpus),
            check_sq=False)
    except ValueError as ex:
        raise InvalidHTSError("HTS file error: " + ex.args[0])

    pysam.set_verbosity(verbosity)
    try:
        yield sam
    finally:
        sam.close()


def get_hts_read_info(min_length: int, skip_read_flag: int, read: pysam.AlignedSegment) -> ReadInfo:

    # Fetch sequence
    seq: str | None = read.get_forward_sequence()
    if not seq:
        raise InputReadError("Error while parsing HTS file: missing sequence!")

    if not is_dna(seq):
        raise InputReadError(f"Invalid sequence: '{seq}'!")

    return ReadInfo(
        sequence=seq,
        is_flagged=(read.flag & skip_read_flag != 0),
        is_qc_fail=read.is_qcfail,
        is_ambiguous=is_ambiguous(seq),
        is_masked=(
            read.cigarstring is not None and
            'S' in read.cigarstring
        ),
        is_short=(len(seq) < min_length),
        is_empty=False) # HTS reads can't be empty? not sure here


def parse_htsfile(
    read_file_info: ReadFileInfo,
    sample: str | None,
    cpus: int,
    reference: str | None = None,
    exclude_qcfail: bool = False,
    exclude_by_len: int | None = None
) -> tuple[LibraryIndependentStats, Counter[str]]:
    min_length: int = exclude_by_len if exclude_by_len is not None else 0
    stats: LibraryIndependentStats
    reads: Counter[str] = Counter()
    skip_read_flag: int = (
        SKIP_READ_FLAG if not exclude_qcfail else
        (SKIP_READ_FLAG | SAM_FLAGS.FQCFAIL)
    )

    def get_read_info(read: pysam.AlignedSegment) -> ReadInfo:
        return get_hts_read_info(min_length, skip_read_flag, read)

    def eval_read(read: pysam.AlignedSegment) -> None:
        read_info: ReadInfo = get_read_info(read)
        if stats.eval_read_info(exclude_qcfail, read_info):
            reads[read_info.sequence] += 1

    with hts_reader(read_file_info.fp, read_file_info.read_mode, cpus, reference) as sam:

        sample_name = _get_hts_sample_name(sam, sample=sample)
        if sample_name is None:
            raise MissingMetadataError(
                "No sample name found in input file header, please provide via '--sample'!")

        stats = LibraryIndependentStats.empty(sample_name)

        try:
            for read in sam.fetch(until_eof=True):
                eval_read(read)

                if stats.total_reads % LOAD_INFO_THRESHOLD == 0:  # pragma: no cover
                    logging.debug(f"Parsed {stats.total_reads} reads, {len(reads)} were unique...")

        except OSError as ex:
            raise InvalidHTSError("HTS file error: " + ex.args[0])

    return stats, reads
