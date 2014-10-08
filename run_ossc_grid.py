import os
import random
import math
import numpy as np

import database.database as database
import jobs.generate as generate
import jobs.kepler_jobs as kepler_jobs

session = database.Session ()
mass_query = session.query (database.SimulationEntry).filter (database.SimulationEntry.binm10 > 14.9).filter (database.SimulationEntry.binm10 < 15.1).filter (database.SimulationEntry.brumoson > 0.0).filter (database.SimulationEntry.woodscon > 0.0)

sc_range = np.arange (-3, 4)
os_range = np.arange (0.1, 1.1, 0.1)

generator = os.path.join(os.path.dirname(os.path.realpath(__file__)), "jobs/generator/s15g")

command = "/Users/justinbrown/Codes/kepler/run/kepler"
run_location = "/Users/justinbrown/Codes/kepler/run/s15"

sets = []

for sc in 2.0 ** sc_range:
    for os in os_range:
        query = mass_query.filter (database.SimulationEntry.osfactor > os * 0.99).filter (database.SimulationEntry.osfactor < os * 1.01)
        query = query.filter (database.SimulationEntry.scpower > sc * 0.99).filter (database.SimulationEntry.scpower < sc * 1.01)
        if query.count () == 0:
            print ("Will run SC = " + str (sc) + ", OS = " + str (os))
            sets.append ((os, sc))
            
random.shuffle (sets)

sims = []

for os, sc in sets:
    name = "s15o" + ("%.1f" % os).lstrip ("0").rstrip ("0") + "%d" % math.log (sc, 2)
    kepler_jobs.run.delay (name, generator, run_location, command, force = True, scpower = sc, osfactor = os)