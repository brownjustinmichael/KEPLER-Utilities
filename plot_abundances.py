import matplotlib.pyplot as plt
import plots.abundances
import records.dump

import sys

# sys.argv.append ("/Users/justinbrown/Codes/kepler/run/s15/s15o1s1#heign")

if (len (sys.argv) < 2):
    print ("Usage: plot_abundances dump_file (output_file)")
    exit ()

record = records.dump.DataDump (sys.argv [1], False)

fig = plt.figure ()
ax = fig.add_axes ([0.1, 0.1, 0.6, 0.75])

abun = plots.abundances.AbundancePlot (ax, record)
plots = abun.plotAll (10.**-4)

ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)

# Alternately, this could be done equivalently with 
# abundances.jTDPlot (sys.argv [1], 10.**-4)
    
if (len (sys.argv) > 2):
    plt.savefig (sys.argv [2])
else:
    # Plot the result
    plt.show ()