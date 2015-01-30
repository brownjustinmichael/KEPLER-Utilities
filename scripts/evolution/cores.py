import sys

import astropy.units as u

import records.dump as dump
import database.cache as cache

if len (sys.argv) < 2:
    print ("Correct usage: python cores.py file-name (core-type-1 core-type-2 ...)")

file = sys.argv [1]

cores = []
if len (sys.argv) > 2:
    for coreType in sys.argv [2:]:
        cores.append (coreType)
# if len (cores) == 0:
#     cores += ["he4", "c12", "o16", "ne20", "si28", "fe56"]

record = dump.DataDump (file)

# for core in cores:
print ("He", cache.calculate_he_core (record).to (u.solMass))
print ("C/O", cache.calculate_co_core (record).to (u.solMass))
print ("Ne", cache.calculate_ne_core (record).to (u.solMass))
print ("Si", cache.calculate_si_core (record).to (u.solMass))
print ("Fe", cache.calculate_fe_core (record).to (u.solMass))