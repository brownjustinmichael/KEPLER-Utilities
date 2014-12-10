import astropy.units as u
import numpy as np

import yields.abundances as ab
from records.dump import Isotope
from yields.integrator import IMFIntegrator

class YieldPlot (object):
    """Plots numerically integrated yields"""
    def __init__(self, yieldReader, imfIntegrator = None):
        super(YieldPlot , self).__init__()
        self.yieldReader = yieldReader
        if imfIntegrator is not None:
            self.imfIntegrator = imfIntegrator
        else:
            self.imfIntegrator = IMFIntegrator (self.yieldReader.get_masses ())
            
    def plot (self, ax, removeIsotopes = None, imfLowerLimit = None, imfUpperLimit = None):
        abundances = self.imfIntegrator.getAbundances (self.yieldReader, imfUpperLimit = imfUpperLimit, imfLowerLimit = imfLowerLimit)

        results = {}
        for isotope in self.yieldReader.isotopes:
            if isotope in ab.solar:
                results [isotope.string] = abundances.productionFactor (isotope)
        
        if removeIsotopes is not None:
            for iso in removeIsotopes:
                if isinstance (iso, Isotope):
                    iso = iso.string
                if iso in results:
                    results.pop (iso)
        
        masses = {}
        pFactors = {}
        isos = {}
        for isotope in self.yieldReader.isotopes:
            if isotope in ab.solar and isotope.string in results:
                if isotope.z not in masses:
                    masses [isotope.z] = []
                    pFactors [isotope.z] = []
                    isos [isotope.z] = isotope
                masses [isotope.z].append (isotope.a)
                pFactors [isotope.z].append (results [isotope.string])

        keys = list (masses.keys ())
        keys.sort ()
        x = [masses [i] for i in keys]
        y = [pFactors [i] for i in keys]
        label = [isos [i].getElementLabel () for i in keys]

        lines = []
        for z in zip (x, y, label):
            lines.append (ax.plot (z [0], z [1], marker = "o") [0])
            ax.annotate (z [2], xy = (z [0] [0], z [1] [0]))

        ax.set_yscale ("log")

        ax.set_xlabel ("Atomic Mass")
        ax.set_ylabel ("Production Factor")
        
        return lines