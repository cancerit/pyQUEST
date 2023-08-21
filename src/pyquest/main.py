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

import logging
import os
import platform
import sys
from typing import NoReturn

import click
from click_option_group import OptionGroup

from . import __version__
from .app_info import AppInfo
from .counting import run_library_independent_counting
from .errors import CustomException, InvalidInputFormatError, InvalidLibraryError, UnsupportedData
from .library import TargetLibrary
from .readers.read_file_info import ReadFileFormat, ReadFileInfo
from .writer import write_stats, write_unique_counts, write_most_common_reads


LOG_LEVELS = [
    'WARNING',
    'INFO',
    'DEBUG'
]

HELP_LIBRARY = "Expanded library definition TSV file with optional headers (common format for single/dual/other)"
HELP_MIN_LENGTH = "Minimum read length"
HELP_MOST_COMMON = "Output top X most common unique read sequences in FASTA format"
HELP_SAMPLE = "Sample name to apply to count column, required for fastq, interrogate header for others when not defined."
HELP_OUTPUT = "Final output to this filename prefix"
HELP_LOW_COUNT = "*.stats.json includes low_count_guides_lt_{15,30}, this option allow specification of an additional cut-off."
HELP_REFERENCE = "Required for CRAM"
HELP_CPUS = "CPUs to use (0 to detect)"


existing_file_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
    resolve_path=True)


class PlatformError(Exception):
    pass


def abort(ex: CustomException) -> NoReturn:
    logging.error(ex.message)
    sys.exit(1)


def get_cpus(cpus: int) -> int:
    if cpus > 0:
        return cpus
    if platform.system() == 'Darwin':
        logging.error("Process affinity detection not available on macOS: please set the CPU count!")
        raise PlatformError
    return len(os.sched_getaffinity(0))  # type: ignore[attr]


def setup_fs(output_dir: str) -> None:
    # TODO: verify file vs. directory as output
    outfolder = os.path.dirname(os.path.abspath(output_dir))
    if not os.path.exists(outfolder):
        os.makedirs(outfolder, exist_ok=True)


lib_dep_opts = OptionGroup("\nLibrary-dependent", help="Options specific to library-dependent counting")
sample_opts = OptionGroup("\nInput sample metadata", help="Options adding information to the input")
perf_opts = OptionGroup("\nPerformance", help="Options to tune the performance")
debug_opts = OptionGroup("\nDebug", help="Options specific to troubleshooting, testing and debugging")


@click.command()
@click.argument('queries', required=True, type=existing_file_path)
@click.option(
    '-o',
    '--output',
    required=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help=HELP_OUTPUT
)
@click.option('--min-length', type=click.IntRange(min=1), default=1, help=HELP_MIN_LENGTH, show_default=True)
@click.option('--most-common', type=click.IntRange(min=1, max=50), default=None, help=HELP_MOST_COMMON, show_default=True)
@sample_opts.option('-s', '--sample', type=str, help=HELP_SAMPLE)
@sample_opts.option('-r', '--reference', type=existing_file_path, help=HELP_REFERENCE)
@lib_dep_opts.option('-l', '--library', default=None, type=existing_file_path, help=HELP_LIBRARY)
@lib_dep_opts.option('--low-count', type=click.IntRange(min=0), default=None, help=HELP_LOW_COUNT, show_default=True)
@perf_opts.option('-c', '--cpus', type=click.IntRange(min=0), default=1, show_default=True, help=HELP_CPUS)
@debug_opts.option(
    '--loglevel',
    default='INFO',
    show_default=True,
    type=click.Choice(LOG_LEVELS, case_sensitive=False),
    help="Set logging verbosity"
)
@debug_opts.option(
    '--no-compression',
    is_flag=True,
    default=False,
    help="Disable output compression"
)
@click.version_option(__version__)
def main(
    library: str | None,
    sample: str | None,
    output: str,
    low_count: int | None,
    most_common: int | None,
    min_length: int,
    reference: str | None,
    queries: str,
    cpus: int,
    loglevel: str,
    no_compression: bool
):
    """
    Count reads and optionally map them to a library.

    QUERIES: Query sequence file (fastq[.gz], sam, bam, cram)
    """

    # Setup logger
    logging.basicConfig(
        level=logging._nameToLevel[loglevel.upper()],
        format="%(levelname)s: %(message)s")

    # Warn on ignored options
    if not library and low_count is not None:
        logging.warning("Low count option ignored in library-independent mode.")

    # Get usable CPU's
    try:
        usable_cpus: int = get_cpus(cpus)
    except PlatformError:
        sys.exit(1)

    # Validate sample name
    read_file_info = ReadFileInfo.from_path(queries)
    if read_file_info.fmt == ReadFileFormat.FASTQ and not sample:
        logging.error("When a FASTQ is provided, a sample name is required!")
        sys.exit(1)

    # Setup output directory
    setup_fs(output)

    app_info = AppInfo(
        __version__,
        ' '.join([os.path.basename(sys.argv[0]), *sys.argv[1:]])
    )

    # Library-independent counting
    logging.info("Loading reads...")
    try:
        unique_count, stats, query_counts = run_library_independent_counting(
            read_file_info, sample, usable_cpus, reference, min_length=min_length)
    except InvalidInputFormatError as ex:
        abort(ex)

    # Write unique counts to file
    logging.info("Writing library-independent counts...")
    write_unique_counts(query_counts, app_info, output, compress=not no_compression)

    # Write X most common unique reads to file if specified
    if most_common is not None:
        write_most_common_reads(most_common, sample, query_counts, output, compress=not no_compression)

    # Warn user with number of 0-length reads
    if stats.total_zero_reads > 0:
        logging.warning(f"{stats.total_zero_reads} zero-length reads")

    if library:

        # Load library
        # TODO: verify stat's are not affected by pre-filtering the targets by length
        logging.info("Loading library...")
        try:
            targets = TargetLibrary.load(library, min_length=min_length)
        except InvalidLibraryError as ex:
            abort(ex)

        logging.info("Finding exact matches and writing library-dependent counts...")
        fp: str = f"{output}.lib_counts.tsv"
        if not no_compression:
            fp += '.gz'
        try:
            stats = targets.map_and_write(
                app_info,
                stats,
                query_counts,
                fp,
                compress=not no_compression,
                custom_count_threshold=low_count)
        except UnsupportedData as ex:
            abort(ex)

    # Write stats to file
    logging.info("Writing stats report...")
    stats_fp: str = f"{output}.stats.json"
    write_stats(app_info, stats, stats_fp)
