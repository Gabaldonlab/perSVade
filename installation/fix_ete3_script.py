#!/usr/bin/env python

# This script fixes an ete3 script that gives problems in version 3.0.0
import sys
import os

# get the script
script = sys.argv[1]

# get the wrong line
wrong_str = 'db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))'
good_str = 'db.execute("INSERT OR REPLACE INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))'

wrong_lines = [l for l in open(script, "r").readlines() if wrong_str in l]

if len(wrong_lines)>1: raise ValueError("there should be only one wrong line")

# get the good lines
good_lines = [l.replace(wrong_str, good_str) for l in open(script, "r").readlines()]

# write script again
script_tmp = "%s.tmp"%script
open(script_tmp, "w").write("".join(good_lines))

os.rename(script_tmp, script)