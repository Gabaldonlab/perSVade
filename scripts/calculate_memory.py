#!/usr/bin/env python

# This function writes the available RAM memory in Gb into a file. It is done in chunks of 3Gb

# module imports
import argparse, os, sys
import psutil

# parse args
parser = argparse.ArgumentParser(description="Calculate memory")
parser.add_argument("--outfile", dest="outfile", required=True, type=str, help="The outfile")
opt = parser.parse_args()

# init object
int_object = int(1e8)*5
unit_str = " "

# keep writing the memory to a file
while True:

    unit_str += " "*int_object
    unit_str_Gb = sys.getsizeof(unit_str)/1e9

    open(opt.outfile, "w").write(str(unit_str_Gb))
