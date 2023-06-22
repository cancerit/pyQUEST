#!/bin/sh

pycodestyle --ignore=E501,W504 src/pyquest
pycodestyle --ignore=E501,W504 tests
