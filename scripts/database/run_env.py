#!/usr/bin/env python

import os
import random
import math
import numpy as np
import astropy.units as u

import kepler_utils
import kepler_utils.jobs.generate as generate
import kepler_utils.jobs.kepler_jobs as kepler_jobs

generator = os.path.join(os.path.dirname(os.path.realpath(kepler_utils.__file__)), "jobs/generator/s15envg")

command = "/Users/justinbrown/Codes/kepler/gfortran/keplery"
run_location = "/Users/justinbrown/Dropbox/Research/Stan/kepler/env"

ls = 10. ** np.arange (37.7, 38.3, 0.1)
mcs = 10. ** np.arange (33.75, 34.05, 0.05) * u.g
ms = u.Quantity ([15.], u.solMass)
rs = 10.0 ** np.arange (10.3, 10.5, 0.1)
habuns = np.arange (0.611, 0.311, -0.1)

sets = []
for habun in habuns:
    for m in ms:
        for r in rs:
            for mc in mcs:
                for l in ls:
                    sets.append ((l, mc, m, r, habun))

for habun in habuns:
    file = open ("env/" + str (habun) + "split", "w")
    
    file.write ("copycomp 1\n")
    file.write ("setcomp 2 %f\n" % (habun))
    file.write ("setcomp 5 %f\n" % (0.985 - habun))
    file.write ("chngcomp %d %d\n" % (1, 700))
    
    file.write ("\n")
    file.close ()

for i, (l, mc, m, r, habun) in enumerate (sets):
    name = "s15e%04d" % (i + 2700)
    kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': False, 'query': False, 'p izonezms': 0, 'p q1faczms': 1000., 'p nstop': 10000, 'p radius0': r, 'p xlum0': l, 'p summ0': mc.to (u.g).value, 'p timezms': 1.0e-5, 'p tstop': 3.15e13, 'p ipup': 0, 'rescalem': str ((m - mc).to (u.solMass).value) + " mult", 'alias zams': "\"link \'%s\'\"" % (str (habun) + "split"), 'tags': ["Envelope", "Habun = " + str (habun), "Split"], "query": True, "goal": "final"}, queue = 'default')
