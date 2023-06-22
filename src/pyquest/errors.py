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

import abc


class CustomException(Exception, abc.ABC):
    def __init__(self, message: str, *args: object) -> None:
        self.message = message
        super().__init__(*args)


class InvalidInputFormatError(CustomException, abc.ABC):
    pass


class InvalidFASTQError(InvalidInputFormatError):
    pass


class InvalidHTSError(InvalidInputFormatError):
    pass


class InputReadError(InvalidInputFormatError):
    pass


class MissingMetadataError(CustomException):
    pass


class InvalidLibraryError(InvalidInputFormatError):
    pass


class UnsupportedData(CustomException):
    pass
