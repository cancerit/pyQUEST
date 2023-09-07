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
from contextlib import contextmanager
import gzip
import json
import logging

from .app_info import AppInfo
from .stats import LibraryIndependentStats


def _get_header_line(k: str, v: str) -> str:
    return f"##{k}: {v}\n"


def _get_fields_line(fields: list[str]) -> str:
    return '#' + '\t'.join(fields) + '\n'


def write_full_header(fh, app_info: AppInfo, fields: list[str]) -> None:
    fh.write(_get_header_line('Command', app_info.command))
    fh.write(_get_header_line('Version', app_info.version))
    fh.write(_get_fields_line(fields))


@contextmanager
def open_output(fp: str, compress: bool = True):
    with (gzip.open if compress else open)(fp, 'wt') as fh:
        yield fh


def write_unique_counts(counts: Counter[str], app_info: AppInfo, output_dir: str, compress: bool = False) -> None:
    """Generate a file with the number of incidents of the same query sequence"""

    fp: str = f"{output_dir}.query_counts.tsv"
    if compress:
        fp += '.gz'

    logging.info(f"Writing query counts file: {fp}")

    with open_output(fp, compress=compress) as fh:

        # Write header
        write_full_header(fh, app_info, [
            'SEQUENCE',
            'LENGTH',
            'COUNT'
        ])

        # Write counts
        for k in sorted(counts.keys()):
            fh.write(f"{k}\t{len(k)}\t{counts[k]}\n")


def write_stats(app_info: AppInfo, stats: LibraryIndependentStats, fp: str) -> None:
    logging.info(f"Writing statistics file: {fp}")
    with open(fp, 'w') as fh:
        json.dump(stats.to_dict(app_info), fh)


def write_most_common_reads(n: int, sample: str | None, counts: Counter[str], output_dir: str, compress: bool = False) -> None:
    fp: str = f"{output_dir}.{sample if sample else 'pyQUEST'}.top{n}.fasta"
    if compress:
        fp += '.gz'

    logging.info(f"Writing most common reads file: {fp}")

    with open_output(fp, compress=compress) as fh:
        for i, (read, count) in enumerate(counts.most_common(n)):
            # if n is greater than the number of reads it will silently write all reads
            fh.write(f"> pyQUEST|{i+1}|{count}\n")
            fh.write(f"{read}\n")
