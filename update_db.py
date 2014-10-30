import sys

import database.database

if len (sys.argv) > 1:
    path = sys.argv [1]
else:
    path = '.' 

database.database.DumpFileEntry.scan_for_updates (path, '*#*')

database.database.CNVFileEntry.scan_for_updates (path, '*.cnv')
