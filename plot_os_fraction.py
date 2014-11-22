import database.database as db
import database.cache as cache

import matplotlib.pyplot as plt
import numpy as np

session = db.Session ()

query = session.query (db.SimulationEntry)

query = query.filter (db.SimulationEntry.binm10 < 21.0).filter (db.SimulationEntry.binm10 > 19.0).filter (db.SimulationEntry.osfactor > 0.99).filter (db.SimulationEntry.complete == True)

dumpsMS = []
for sim in query.all ():
    nHDep = sim.getStateDump ("hdep").ncyc
    for dump in sim.dumpfiles:
        if dump.ncyc <= nHDep:
            dumpsMS.append (dump)
            
dumpsHe = []
for sim in query.all ():
    nHDep = sim.getStateDump ("hdep").ncyc
    nHeDep = sim.getStateDump ("hedep").ncyc
    for dump in sim.dumpfiles:
        if dump.ncyc <= nHeDep and dump.ncyc > nHDep:
            dumpsHe.append (dump)
            
coreFracs = [dump.cache (session, "core_os", cache.core_overshoot) for dump in dumpsMS]
coreFracsHe = [dump.cache (session, "core_os", cache.core_overshoot) for dump in dumpsHe]
envFracsHe = [dump.cache (session, "env_os", cache.env_overshoot) for dump in dumpsHe]
weights = [dump.dt for dump in dumpsMS]
weightsHe = [dump.dt for dump in dumpsHe]

plt.figure ()
ax = plt.subplot (111)

plt.hist (coreFracs, weights = weights, normed = True, stacked = True, bins = 20)
ax.set_xlabel ("Main Sequence Core Overshoot Length ($H_{p}$)")

plt.figure ()
ax = plt.subplot (111)

plt.hist (coreFracsHe, weights = weightsHe, normed = True, stacked = True, bins = 20)
ax.set_xlabel ("Helium Burning Core Overshoot Length ($H_{p}$)")

plt.figure ()
ax = plt.subplot (111)

plt.hist (envFracsHe, weights = weightsHe, normed = True, stacked = True, bins = 20)
ax.set_xlabel ("Helium Burning Envelope Overshoot Length ($H_{p}$)")

plt.show ()