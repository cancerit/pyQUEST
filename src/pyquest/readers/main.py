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
import gzip
import logging
from time import time

import magic

from ..errors import MissingMetadataError
from ..stats import LibraryIndependentStats
from .read_file_info import ReadFileFormat, ReadFileInfo
from .hts import parse_htsfile
from .fq import parse_fastq


def is_gzip(seq_file: str) -> bool:
    magic_types = magic.from_file(seq_file)
    return (
        "gzip compressed data" in magic_types or
        "gzip compatible" in magic_types
    )


def parse_reads(
    read_file_info: ReadFileInfo,
    sample: str | None,
    cpus: int,
    reference: str | None = None,
    exclude_qcfail: bool = False,
    exclude_by_len: int | None = None
) -> tuple[int, LibraryIndependentStats, Counter[str]]:
    """
    This function is for the initial collation of unique read sequences in the original orientation only (hts will do revcomp).
    There is no chunking of data, so this relies on large memory lookups at present.

    Selecting correct underlying parser is via file extension:
    - cram/bam/sam -> htslib processing
    - gz assume gzip FASTQ
    - anything else assume uncompressed FASTQ
    """

    start = time()
    match read_file_info.fmt:

        case ReadFileFormat.HTS:
            logging.info(f"Sequence input detected as *{read_file_info.ext}")
            stats, counts = parse_htsfile(
                read_file_info,
                sample,
                cpus,
                reference=reference,
                exclude_qcfail=exclude_qcfail,
                exclude_by_len=exclude_by_len)

        case ReadFileFormat.FASTQ:
            is_compressed: bool = is_gzip(read_file_info.fp)

            logging.info(
                "Sequence input detected as gzip (assume FASTQ)" if is_compressed else
                "Uncompressed data (assume FASTQ)")

            if not sample:
                raise MissingMetadataError("Sample name required for FASTQ inputs!")

            with (gzip.open if is_compressed else open)(read_file_info.fp, 'rt') as fq_fh:
                stats, counts = parse_fastq(
                    fq_fh,
                    sample,
                    exclude_qcfail=exclude_qcfail,
                    exclude_by_len=exclude_by_len)

        case _:
            raise RuntimeError("Invalid reads format!")

    unique_count: int = len(counts)
    logging.info(f"Parsed {stats.total_reads} reads, {unique_count} were unique...")
    logging.info(f"Read parsing took: {int(time() - start)}s")

    return unique_count, stats, counts
