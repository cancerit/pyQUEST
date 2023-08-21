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

from dataclasses import dataclass
from typing import Any

from .app_info import AppInfo
from .readers.read_info import ReadInfo


def _round_float(x: float) -> float:
    return round(x, ndigits=2)


@dataclass
class LibraryIndependentStats:

    # Optional metadata
    sample_name: str

    # Reads
    input_reads: int
    total_reads: int
    vendor_failed_reads: int
    length_excluded_reads: int
    ambiguous_nt_reads: int
    masked_reads: int
    total_invalid_reads: int
    total_excluded_reads: int
    total_zero_reads: int

    @classmethod
    def empty(cls, sample_name: str):
        return cls(sample_name, 0, 0, 0, 0, 0, 0, 0, 0)

    def eval_read_info(self, discard_qc: bool, read_info: ReadInfo) -> bool:
        if read_info.is_qc_fail:
            self.vendor_failed_reads += 1
        if read_info.is_ambiguous:
            self.ambiguous_nt_reads += 1
        if read_info.is_masked:
            self.masked_reads += 1
        if read_info.is_sequence_invalid:
            self.total_invalid_reads += 1
        if read_info.is_short:
            self.length_excluded_reads += 1
        if read_info.is_empty:
            self.total_zero_reads += 1

        self.input_reads += 1
        keep: bool = not read_info.to_discard(discard_qc)
        if keep:
            self.total_reads += 1
        else:
            self.total_excluded_reads += 1
        return keep

    def to_dict(self, app_info: AppInfo) -> dict[str, Any]:
        return {
            **app_info.to_dict(),
            'sample_name': self.sample_name,
            'input_reads': self.input_reads,
            'total_reads': self.total_reads,
            'discarded_reads': self.total_excluded_reads,
            'vendor_failed_reads': self.vendor_failed_reads,
            'length_excluded_reads': self.length_excluded_reads,
            'ambiguous_nt_reads': self.ambiguous_nt_reads,
            'masked_reads': self.masked_reads,
            'zero_length_reads': self.total_zero_reads,
        }


@dataclass
class LibraryDependentStats(LibraryIndependentStats):

    # Reads
    mapped_to_template_reads: int
    mean_count_per_template: float
    median_count_per_template: float
    multimap_reads: int
    unmapped_reads: int

    # Templates
    total_templates: int
    total_unique_templates: int
    length_excluded_templates: int
    zero_count_templates: int
    low_count_templates_lt_15: int
    low_count_templates_lt_30: int
    low_count_templates_user: tuple[int, int] | None

    gini_coefficient: float

    def to_dict(self, app_info: AppInfo) -> dict[str, Any]:
        return {
            **super().to_dict(app_info),
            'mapped_to_template_reads': self.mapped_to_template_reads,
            'mean_count_per_template': _round_float(self.mean_count_per_template),
            'median_count_per_template': _round_float(self.median_count_per_template),
            'multimap_reads': self.multimap_reads,
            'unmapped_reads': self.unmapped_reads,
            'total_templates': self.total_templates,
            'total_unique_templates': self.total_unique_templates,
            'length_excluded_templates': self.length_excluded_templates,
            'zero_count_templates': self.zero_count_templates,
            'low_count_templates_lt_15': self.low_count_templates_lt_15,
            'low_count_templates_lt_30': self.low_count_templates_lt_30,
            'low_count_templates_user': {
                'lt': self.low_count_templates_user[0],
                'count': self.low_count_templates_user[1]
            } if self.low_count_templates_user is not None else None,
            'gini_coefficient': _round_float(self.gini_coefficient)
        }
