#!/usr/bin/env python

import sys

from kepler_utils.database.database import DumpFileEntry, CNVFileEntry

if len (sys.argv) > 1:
    path = sys.argv [1]
else:
    path = '.' 
    
if len (sys.argv) > 2:
    globstring = sys.argv [2]
else:
    globstring = "*"

if len (sys.argv) > 3:
    tags = sys.argv [3:]
else:
    tags = []
    
DumpFileEntry.scan_for_updates (path, globstring + '#*', tags = tags, log_info = True)

CNVFileEntry.scan_for_updates (path, globstring + '.cnv', tags = tags)
