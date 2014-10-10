import matplotlib.pyplot as plt
import plots.abundances
import records.dump

import sys

sys.argv.append ("/Users/justinbrown/Codes/kepler/run/s15/s15o1s1#heign")

if (len (sys.argv) < 2):
    print ("Usage: plot_abundances dump_file (output_file)")
    exit ()

record = records.dump.DataDump (sys.argv [1], False)

fig = plt.figure ()
ax = fig.add_axes ([0.1, 0.1, 0.6, 0.75])

abun = plots.abundances.AbundancePlot (ax, record)
plots = abun.plotAll (10.**-4)

ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)

fig = plt.figure ()
ax = plt.subplot (111)

record2 = records.dump.DataDump ("/Users/justinbrown/Codes/kepler/run/s15/s15o0s1#heign", False)

print (record.parameters ['teff'])
print (record2.parameters ['teff'])

print (record.parameters ['radius'])
print (record2.parameters ['radius'])

ax.plot (record ['dn'], record ['tn'], label = "s15o.12")
ax.plot (record2 ['dn'], record2 ['tn'], label = "s15o.9-3")
ax.set_xscale ('log')
ax.set_yscale ('log')

ax.legend ()

# Alternately, this could be done equivalently with 
# abundances.jTDPlot (sys.argv [1], 10.**-4)
    
if (len (sys.argv) > 2):
    plt.savefig (sys.argv [2])
else:
    # Plot the result
    plt.show ()