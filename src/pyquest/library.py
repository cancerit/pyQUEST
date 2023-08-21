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
from dataclasses import dataclass
import logging
from typing import Final

import numpy as np

from .app_info import AppInfo
from .errors import InvalidLibraryError
from .stats import LibraryDependentStats, LibraryIndependentStats
from .utils import get_stats
from .writer import open_output, write_full_header


LIB_FIELD_ID: int = 0
LIB_FIELD_NAME: int = 1
LIB_FIELD_SEQ: int = 2


@dataclass(slots=True)
class Target:
    id: str
    name: str
    seq: str


@dataclass(slots=True)
class TargetLibrary:
    min_length: int
    total_target_count: int
    short_target_count: int
    target_counts: Counter[str]
    targets: list[Target]

    @property
    def unique_target_count(self) -> int:
        return len(self.target_counts)

    @classmethod
    def load(cls, fp: str, min_length: int = 0):
        if min_length < 0:
            raise InvalidLibraryError("Invalid minimum length!")

        total_target_count: int = 0
        target_counts: Counter[str] = Counter()
        targets: list[Target] = []
        target: Target
        short_sequences: set[str] = set()

        with open(fp) as fh:

            # Skip header
            lib_line = fh.readline().strip()
            while lib_line[0] == '#':
                lib_line = fh.readline().strip()
                if not lib_line:
                    # TODO: handle correctly
                    raise RuntimeError("Empty library!")

            while lib_line:
                total_target_count += 1
                t = lib_line.strip().split('\t', maxsplit=3)
                target = Target(t[LIB_FIELD_ID], t[LIB_FIELD_NAME], t[LIB_FIELD_SEQ])
                target_counts[target.seq] += 1
                if len(target.seq) < min_length:
                    short_sequences.add(target.seq)
                targets.append(target)
                lib_line = fh.readline().strip()

        short_target_count: int = len(short_sequences)
        return cls(min_length, total_target_count, short_target_count, target_counts, targets)

    def map_and_write(
        self,
        app_info: AppInfo,
        stats: LibraryIndependentStats,
        queries: Counter[str],
        fp: str,
        custom_count_threshold: int | None = None,
        compress: bool = False
    ) -> LibraryDependentStats:
        """
        Assign the query sequence counts to the matching library sequences
        while collecting statistics and list the resulting mapping to file
        """
        count: int

        if self.short_target_count > 0:
            logging.warning(f"{self.short_target_count} unique library sequences are below the minimum length ({self.min_length})!")

        with open_output(fp, compress=compress) as fh:

            # Write header to library-dependent count file
            write_full_header(fh, app_info, [
                'ID',
                'NAME',
                'SEQUENCE',
                'LENGTH',
                'COUNT',
                'UNIQUE',
                'SAMPLE'
            ])

            sample_name: Final[str] = stats.sample_name
            low_counts: Counter[int] = Counter()

            for target in self.targets:
                count = queries.get(target.seq, 0)

                # Evaluate the number of synonymous targets
                is_unique: int = 1 if self.target_counts[target.seq] == 1 else 0

                # Write to library-dependent count file
                fh.write(f"{target.id}\t{target.name}\t{target.seq}\t{len(target.seq)}\t{count}\t{is_unique}\t{sample_name}\n")

        # Initialise stats
        multimap_reads: int = 0
        mapped_to_template_reads: int = 0
        mean_count_per_template: float = 0.0
        median_count_per_template: float = 0.0
        gini_coefficient: float = 0.0

        # Verify at least one target is passing the length filter
        long_target_count: int = self.unique_target_count - self.short_target_count
        if long_target_count > 0:

            # Preallocate the read counts
            template_counts = np.zeros(long_target_count, dtype=np.uint64)

            i: int = 0
            for seq, n in self.target_counts.items():
                if len(seq) >= self.min_length:
                    count = queries.get(seq, 0)
                    template_counts[i] = count
                    if n > 1:
                        multimap_reads += count
                    i += 1

            # Sort the counts (required by `get_stats`)
            template_counts.sort()

            # Count templates with low read counts
            # TODO: take better advantage of the sorting...?
            low_counts[0] = np.count_nonzero(template_counts == 0)
            count_thresholds: list[int] = [15, 30]
            if custom_count_threshold is not None:
                count_thresholds.append(custom_count_threshold)
            for t in count_thresholds:
                low_counts[t] = np.count_nonzero(template_counts < t)

            # Generate all stats
            mapped_to_template_reads, mean_count_per_template, median_count_per_template, gini_coefficient = get_stats(
                template_counts, gini_corr=False)

        # Populate library-dependent stats
        unmapped_reads: int = stats.total_reads - mapped_to_template_reads
        return LibraryDependentStats(

            # Library-independent stats
            sample_name=stats.sample_name,
            input_reads=stats.input_reads,
            total_reads=stats.total_reads,
            vendor_failed_reads=stats.vendor_failed_reads,
            length_excluded_reads=stats.length_excluded_reads,
            ambiguous_nt_reads=stats.ambiguous_nt_reads,
            masked_reads=stats.masked_reads,
            total_invalid_reads=stats.total_invalid_reads,
            total_excluded_reads=stats.total_excluded_reads,
            total_zero_reads=stats.total_zero_reads,

            mapped_to_template_reads=mapped_to_template_reads,
            multimap_reads=multimap_reads,
            unmapped_reads=unmapped_reads,
            total_templates=self.total_target_count,
            total_unique_templates=self.unique_target_count,
            length_excluded_templates=self.short_target_count,

            # Derived stats
            mean_count_per_template=mean_count_per_template,
            median_count_per_template=median_count_per_template,
            gini_coefficient=gini_coefficient,

            # Count thresholds
            zero_count_templates=low_counts.get(0, 0),
            low_count_templates_lt_15=low_counts.get(15, 0),
            low_count_templates_lt_30=low_counts.get(30, 0),
            low_count_templates_user=(
                (
                    custom_count_threshold,
                    low_counts.get(custom_count_threshold, 0)
                ) if custom_count_threshold is not None else
                None
            ))
