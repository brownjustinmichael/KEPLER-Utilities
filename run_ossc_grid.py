import os
import random
import math
import numpy as np

import jobs.generate as generate
import jobs.kepler_jobs as kepler_jobs

sets = []

semirange = np.arange (-3,4,1)
oshtrange = np.arange (0.1,1.1,0.1)

for semi in semirange:
    for osht in oshtrange:
        sets.append ((osht, semi))

random.shuffle (sets)
            
generator = os.path.join(os.path.dirname(os.path.realpath(__file__)), "jobs/generator/s15hg")

command = "/Users/justinbrown/Codes/kepler/run/kepler"
run_location = "/Users/justinbrown/Codes/kepler/run/s15"

for osht, semi in sets:
    name = "s15n" + str (osht).replace ('0', '') + str (semi)
    kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': False, 'scpower': 2.0 ** semi, 'osfactor': osht, 'tags': ['OS/SC Grid', 'Stabilized']}, queue = 'default')

# for sc in np.arange (-8, 1, 2):
#     name = "s20hs" + ("%d" % sc)
#     kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': False, 'scpower': 0.0, 'scfactor': 10.0 ** sc, 'osherwig': 0.01, 'tags': ['For Kevin'], 'query': True}, queue = 'priority')
