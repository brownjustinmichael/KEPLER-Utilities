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

l = 10. ** 38.3
mc = 10. ** 33.90 * u.g - 1.0 * u.solMass
m = 15 * u.solMass
r = 10. ** 10.1
n = 0
name = "s15ea%03d" % 23

file = open ("env/" + str (n) + "comp", "w")

ninner = 400

for i in range (n + ninner, ninner, -1):
    file.write ("copycomp 1\n")
    file.write ("setcomp 2 %f\n" % (0.711 * (i - ninner) / n))
    file.write ("setcomp 4 %f\n" % (0.0))
    file.write ("setcomp 5 %f\n" % (0.985 - 0.711 * (i - ninner) / n))
    file.write ("chngcomp %d %d\n" % (i, i + 1))
    
for i in range (ninner, 0, -1):
    file.write ("copycomp 1\n")
    file.write ("setcomp 2 %f\n" % (0.0))
    file.write ("setcomp 4 %f\n" % (0.0))
    file.write ("setcomp 5 %f\n" % (0.985))
    file.write ("chngcomp %d %d\n" % (i, i + 1))
    
# file.write ("cutsurf -3\n")
# file.write ("newe\n")
file.close ()

kepler_jobs.run.apply_async ([name, generator, run_location, command], kwargs = {'force': False, 'query': True, 'p izonezms': 0, 'p q1faczms': 1000.0, 'p nstop': 10000, 'p radius0': r, 'p xlum0': l, 'p summ0': mc.to (u.g).value, 'p timezms': 1.0e-5, 'p tstop': 3.15e13, 'p ipup': 0, 'rescalem': str ((m - mc).to (u.solMass).value) + " mult", 'alias zams': "\"link \'%s\'\"" % (str (n) + "comp"), 'tags': ["Envelope", "NGrad = " + str (n), "Test"], "query": True, "goal": "final"}, queue = 'priority')

