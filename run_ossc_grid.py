import os
import random
import math
import numpy as np

import jobs.generate as generate
import jobs.kepler_jobs as kepler_jobs

sc_range = np.arange (-3, 4)
os_range = np.arange (0.1, 1.1, 0.1)

sets = []

for sc in 2.0 ** sc_range:
    for osh in os_range:
        sets.append ((osh, sc))
            
random.shuffle (sets)
 
generator = os.path.join(os.path.dirname(os.path.realpath(__file__)), "jobs/generator/s15g")

command = "/Users/justinbrown/Codes/kepler/run/kepler"
run_location = "/Users/justinbrown/Codes/kepler/run/s15"

for os, sc in sets:
    name = "s15o" + ("%.1f" % os).lstrip ("0").rstrip ("0") + "%d" % math.log (sc, 2)
    kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': False, 'scpower': sc, 'osfactor': os, 'tags': ['OS/SC Grid'], 'query': True}, queue = 'default')

# for sc in np.arange (-8, 1, 2):
#     name = "s20hs" + ("%d" % sc)
#     kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': False, 'scpower': 0.0, 'scfactor': 10.0 ** sc, 'osherwig': 0.01, 'tags': ['For Kevin'], 'query': True}, queue = 'priority')
