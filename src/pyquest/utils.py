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
import numpy as np

from .errors import UnsupportedData


def get_median_from_sorted(sa: np.ndarray) -> float:
    """
    Calculate the median, assuming the input is sorted
    """

    n: int = sa.shape[0]

    if n == 0:
        raise UnsupportedData("Median: empty array!")

    if n == 1:
        return float(sa[0])

    # Even length
    if n % 2 == 0:
        i: int = n // 2
        return float(sa[[i - 1, i]].mean())

    # Odd length
    return float(sa[(n - 1) // 2])


def get_stats(sa: np.ndarray, gini_corr: bool = False) -> tuple[int, float, float, float]:
    """
    Calculate the following statistics, assuming the array is sorted:
    - total
    - mean
    - median
    - Gini coefficient
    """

    n: int = sa.shape[0]

    if n == 0:
        raise UnsupportedData("Stats calculation: empty array!")

    m: np.uint64 = sa.sum()

    if m == 0:
        logging.warning("No library matches!")
        return 0, 0.0, 0.0, 0.0

    # Mean
    mean: float = int(m) / n

    # Median
    median: float = get_median_from_sorted(sa)

    # Gini coefficient
    gini_coeff: float = float(
        (2 * np.sum(sa * np.arange(1, n + 1)) / m - (n + 1)) /
        ((n - 1) if gini_corr else n)
    )

    return int(m), mean, median, gini_coeff
