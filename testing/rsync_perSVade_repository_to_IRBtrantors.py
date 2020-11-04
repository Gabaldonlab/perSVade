#!/usr/bin/env python

# module imports
import os
import sys
import pandas as pd

# define the parent dir. It has to be the bsc machine
ParentDir = "%s/samba"%(os.getenv("HOME")); # local
if os.path.exists(ParentDir): threads = 4
else: raise ValueError("This has to be run from the BSC machine")

# This will rsync all the files from the perSVade repository into the trantors of the IRB.