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
from tempfile import NamedTemporaryFile

from pyquest.app_info import AppInfo
from pyquest.library import Target, TargetLibrary
from pyquest.stats import LibraryIndependentStats


def parse_lib_dep_counts(fp):
    def parse_row(s):
        id, name, seq, seq_len, count, unique, sample = s.rstrip().split('\t')
        return id, name, seq, int(seq_len), int(count), bool(int(unique)), sample

    with open(fp) as fh:
        return [
            parse_row(line)
            for line in fh.readlines()
            if line[0] != '#'
        ]


def library_from_dict(d):
    seq_count = Counter({k: len(v) for k, v in d.items()})
    return TargetLibrary(0, sum(seq_count.values()), 0, seq_count, [
        Target(str(j), name, seq)
        for seq, names in d.items()
        for j, name in enumerate(names)
    ])


def test_library_map():

    # Initialise library
    d = {
        'AAA': {'p1', 'p2'},
        'CCC': {'q1'},
        'GGG': {'r1'},
        'TTT': {'s1', 's2', 's3'}
    }
    library = library_from_dict(d)
    print(library)
    assert library.unique_target_count == len(d)

    # Set query read counts
    queries = Counter({
        'AAA': 5,
        'GGG': 2,
        'TTT': 10,
        'GTG': 50,
        'TGT': 50
    })

    # Initialise the library-independent stats
    lib_indep_stats = LibraryIndependentStats.empty('X')
    lib_indep_stats.total_reads = sum(queries.values())

    app_info = AppInfo('1.0.0', 'pyquest')

    with NamedTemporaryFile() as temp_file:

        # Map queries to targets
        stats = library.map_and_write(
            app_info,
            lib_indep_stats,
            queries,
            temp_file.name,
            compress=False)

        # Load library-dependent counts
        counts = parse_lib_dep_counts(temp_file.name)

    # Evaluate the library-dependent stats
    matches = 3
    assert stats.zero_count_templates == 1
    assert stats.low_count_templates_lt_15 == matches + 1
    assert stats.low_count_templates_lt_30 == matches + 1
    assert stats.unmapped_reads == 100
    assert stats.total_unique_templates == len(d)

    # Evaluate the matching read counts
    assert len(counts) == library.total_target_count
    for id, name, seq, seq_len, count, unique, sample in counts:
        assert id.isdigit()
        assert count == queries[seq]
        assert unique == (library.target_counts[seq] == 1)
        assert sample == lib_indep_stats.sample_name
        assert len(seq) == seq_len
