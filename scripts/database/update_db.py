import sys

import database.database

if len (sys.argv) > 1:
    path = sys.argv [1]
else:
    path = '.' 
    
if len (sys.argv) > 2:
    tags = sys.argv [2:]
else:
    tags = []
    
database.database.DumpFileEntry.scan_for_updates (path, '*#*', tags = tags)

database.database.CNVFileEntry.scan_for_updates (path, '*.cnv', tags = tags)
