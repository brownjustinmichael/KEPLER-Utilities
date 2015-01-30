import sys

from kepler_utils.database.database import DumpFileEntry, CNVFileEntry

if len (sys.argv) > 1:
    path = sys.argv [1]
else:
    path = '.' 
    
if len (sys.argv) > 2:
    tags = sys.argv [2:]
else:
    tags = []
    
DumpFileEntry.scan_for_updates (path, '*#*', tags = tags)

CNVFileEntry.scan_for_updates (path, '*.cnv', tags = tags)
