import records.dump
import matplotlib.pyplot as plt
import astropy.units as u
import numpy
from itertools import cycle

class AbundancePlot (object):
    """docstring for AbundancePlot """
    def __init__ (self, axis, record):
        super(AbundancePlot , self).__init__()
        self.axis = axis
        self.record = record
        self.lines = ["-","--","-.",":"]
        self.linecycler = cycle (self.lines)
        
    def plotAbundance (self, isotope, **kwargs):
        if not isinstance (isotope, records.dump.Isotope):
            isotope = self.record.getIsotope (isotope)
        return self.axis.plot (self.record ['mass coordinate'].to (u.solMass), self.record [str (isotope)], next (self.linecycler), label = isotope.getLabel (), **kwargs) [0]
    
    def plotAll (self, ymin = 10.**-2.5, **kwargs):
        self.axis.set_ylim (ymin,10.**0)
        self.axis.set_yscale ('log')

        self.axis.set_xlabel ("Mass (Solar Masses)")
        self.axis.set_ylabel ("Abundance")
        
        maxes = {}
        
        for isotope in self.record.getIsotopes ():
            newmax = numpy.max (self.record [str (isotope)])
            if newmax > ymin:
                maxes [str (isotope)] = numpy.max (self.record [str (isotope)])
        maxkey = []
        while maxes != {}:
            current = 0.0
            currentkey = None
            for key in maxes:
                if maxes [key] > current:
                    currentkey = key
                    current = maxes [key]
            maxkey.append (currentkey)
            maxes.pop (currentkey)
        
        plots = {}
        for isotope in self.record.getIsotopes ():
            if str (isotope) in maxkey:
                plots [str (isotope)] = self.plotAbundance (isotope, **kwargs)

        return plots

def jTDPlot (record, ymin = 10.**-2.5, **kwargs):
    if not isinstance (record, records.dump.DataDump):
        record = records.dump.DataDump (record, False)
    
    fig = plt.figure ()
    ax = fig.add_axes ([0.1, 0.1, 0.6, 0.75])

    abun = AbundancePlot (ax, record)
    plots = abun.plotAll (ymin)

    ax.legend (bbox_to_anchor=(1.05, 1.2), loc=2, borderaxespad=1)
    
    return fig