#-------------------------------------------------------------------------------
# Name:        ppgfunctions
# Purpose:     general purpose functions
#
# Author:      Pavel Polishchuk
#
# Created:     11.01.2013
# Copyright:   (c) Pavel Polishchuk 2013-2015
# Licence:     BSD 3-clause
#-------------------------------------------------------------------------------

import os


def SortTwoLists(primary, secondary):
    # sort two lists by order of elements of the primary list
    paired_sorted = sorted(zip(primary, secondary), key=lambda x: x[0])
    return (map(list, zip(*paired_sorted)))  # two lists


def GetFileNameNoExt(path):
    return (os.path.splitext(os.path.basename(path)))[0]


def GetFileExt(path):
    return (os.path.splitext(os.path.basename(path)))[1]


def GetWorkDir(filename):
    return (os.path.abspath(os.path.dirname(filename)))
