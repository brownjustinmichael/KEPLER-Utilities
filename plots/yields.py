import astropy.units as u
import numpy as np

import yields.abundances as ab
from records.dump import Isotope
from yields.integrator import IMFIntegrator

class YieldPlot (object):
    """Plots numerically integrated yields"""
    def __init__(self, yieldReader, imfIntegrator = None, windReader = None, windIntegrator = None):
        super(YieldPlot , self).__init__()
        self.yieldReader = yieldReader
        if imfIntegrator is not None:
            self.imfIntegrator = imfIntegrator
        else:
            self.imfIntegrator = IMFIntegrator (self.yieldReader.get_masses ())
        
        if windReader is not None:
            self.windReader = windReader
        else:
            self.windReader = self.yieldReader
            
        if windIntegrator is not None:
            self.windIntegrator = windIntegrator
        else:
            self.windIntegrator = IMFIntegrator (self.windReader.get_masses ())
            
    def plot (self, ax, windMultiplier = 1.0, removeIsotopes = None, imfLowerLimit = None, imfUpperLimit = None):
        yd = {}
        mask = None
        if imfUpperLimit is not None:
            mask = self.windReader.get_masses () < imfUpperLimit * 1.01
        if imfLowerLimit is not None:
            if mask is None:
                mask = self.windReader.get_masses () > imfLowerLimit * 0.99
            else:
                mask = np.logical_and (mask, self.windReader.get_masses () > imfLowerLimit * 0.99)
            
        for isotope in self.yieldReader.isotopes:
            yd [isotope.string] = self.imfIntegrator (self.yieldReader.get_yield (isotope), mask = mask)
            
            if isotope in self.windReader.isotopes:
                yd [isotope.string] += windMultiplier * self.windIntegrator (self.windReader.get_wind (isotope))

        abundances = ab.Abundances (yd)

        results = {}
        for isotope in self.yieldReader.isotopes:
            if isotope in ab.solar:
                results [isotope.string] = abundances.productionFactor (isotope)
        
        if removeIsotopes is not None:
            for iso in removeIsotopes:
                if isinstance (iso, Isotope):
                    iso = iso.string
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

        for z in zip (x, y, label):
            ax.plot (z [0], z [1], marker = "o")
            ax.annotate (z [2], xy = (z [0] [0], z [1] [0]))

        ax.set_yscale ("log")

        ax.set_ylim ((0.1,100))
        ax.set_xlim ((0.0, 120))

        ax.set_xlabel ("Atomic Mass")
        ax.set_ylabel ("Production Factor")