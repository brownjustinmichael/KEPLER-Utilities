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

ls = 10. ** np.arange (37.7, 38.5, 0.1)
mcs = 10. ** np.arange (33.75, 34.05, 0.05) * u.g
ms = 10. ** np.arange (34.45, 34.55, 0.05) * u.g
rs = 10.0 ** np.arange (10.3, 10.5, 0.1)
ns = [i * 200 for i in range (4)]

sets = []
for mc in mcs:
    for m in ms:
        for r in rs:
            for n in ns:
                for l in ls:
                    sets.append ((l, mc, m, r, n))
    #             break
    #         break
    #     break
    # break

for n in ns:
    file = open ("env/" + str (n) + "comp", "w")
    
    # file.write ("cutsurf -3\n")
    file.write ("\n")
    
    for i in range (n, 0, -1):
        file.write ("copycomp 1\n")
        file.write ("setcomp 2 %f\n" % (0.711 * i / n))
        file.write ("setcomp 5 %f\n" % (0.985 - 0.711 * i / n))
        file.write ("chngcomp %d %d\n" % (i, i + 1))
        
    # file.write ("newe\n")
    file.close ()

for i, (l, mc, m, r, n) in enumerate (sets):
    name = "s15e%04d" % i
    kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': True, 'query': True, 'p izonezms': 0, 'p nstop', 3000, 'p radius0': r, 'p xlum0': l, 'p summ0': mc.to (u.g).value, 'p timezms': 1.0e6, 'p tstop': 3.15e13, 'p ipup': 0, 'rescalem': str ((m - mc).to (u.solMass).value) + " mult", 'alias zams': "\"link \'%s\'\"" % (str (n) + "comp"), 'tags': ["Envelope", "NGrad = " + str (n)], "query": True, "goal": "final"}, queue = 'default')
