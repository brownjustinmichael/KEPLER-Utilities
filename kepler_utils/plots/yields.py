import astropy.units as u
import numpy as np
import re

from kepler_utils.yields.abundances import solar
from kepler_utils.records.dump import Isotope
from kepler_utils.yields.integrator import IMFIntegrator

class YieldPlot (object):
    """Plots numerically integrated yields"""
    def __init__(self, yieldReader, imfIntegrator = None):
        super(YieldPlot , self).__init__()
        self.yieldReader = yieldReader
        if imfIntegrator is not None:
            self.imfIntegrator = imfIntegrator
        else:
            self.imfIntegrator = IMFIntegrator (self.yieldReader.get_masses ())
            
    def getData (self, removeIsotopes = None, imfLowerLimit = None, imfUpperLimit = None, record = None, elements = True, relativeIso = None, **kwargs):
        abundances = self.imfIntegrator.getAbundances (self.yieldReader, imfUpperLimit = imfUpperLimit, imfLowerLimit = imfLowerLimit)

        results = {}
        tocheck = {}
        
        if elements:
            for isotope in self.yieldReader.isotopes:
                if isotope in solar:
                    if isotope.z not in tocheck or solar [tocheck [isotope.z]] < solar [isotope]:
                        tocheck [isotope.z] = isotope
                    else:
                        continue
            tocheck = [tocheck [key] for key in tocheck]
        else:
            tocheck = self.yieldReader.isotopes
        
        for isotope in tocheck:
            results [isotope.string] = abundances.productionFactor (isotope, element = elements, relativeIso = relativeIso)
        
        if removeIsotopes is not None:
            for iso in removeIsotopes:
                if isinstance (iso, Isotope):
                    iso = iso.string
                if iso in results:
                    results.pop (iso)
        
        masses = {}
        pFactors = {}
        isos = {}
        if record is not None:
            output = open (record, "w")
            output.write ("# iso z a proFac\n")
        for isotope in self.yieldReader.isotopes:
            if isotope in solar and isotope.string in results:
                if isotope.z not in masses:
                    masses [isotope.z] = []
                    pFactors [isotope.z] = []
                    isos [isotope.z] = isotope
                masses [isotope.z].append (isotope.a)
                pFactors [isotope.z].append (results [isotope.string])
                if record is not None:
                    output.write ("%s %i %i %f\n" % (isotope.string, isotope.z, isotope.a, results [isotope.string]))
        if record is not None:
            output.close ()

        keys = list (masses.keys ())
        keys.sort ()
        if elements:
            x = [[i] * len (masses [i]) for i in keys]
        else:
            x = [masses [i] for i in keys]
        y = [pFactors [i] for i in keys]
        label = [isos [i].getElementLabel () for i in keys]
        
        return x, y, label
            
    def plot (self, ax, removeIsotopes = None, imfLowerLimit = None, imfUpperLimit = None, record = None, names = True, elements = True, minDist = 1.5, relativeIso = None, moveLabels = None, **kwargs):
        x, y, label = self.getData (removeIsotopes = removeIsotopes, imfLowerLimit = imfLowerLimit, imfUpperLimit = imfUpperLimit, record = record, elements = elements, relativeIso = relativeIso, **kwargs)
        
        lines = []
        if moveLabels is None:
            moveLabels = {}
        marker = kwargs.pop ("marker", "o")
        
        if elements:
            lines.append (ax.plot (np.array (x).flat, np.array (y).flat, marker = marker, **kwargs) [0])
        
        annotationLocations = []
        
        for z in zip (x, y, label):
            if not elements:
                lines.append (ax.plot (z [0], z [1], marker = marker, **kwargs) [0])
            if names:
                i = np.argmax (z [1])
                current = (z [0] [i] - 1, z [1] [i] * 2)
                for label in moveLabels:
                    if re.sub ('\}', '', re.sub ("\\$", "", re.sub ('\\\\mathrm\{', '', z [2]))) == re.sub ('\}', '', re.sub ("\\$", "", re.sub ('\\\\mathrm\{', '', label))):
                        current = (current [0] + moveLabels [label] [0], current [1] * 10.0 ** (moveLabels [label] [1] / 10.))
                while len (annotationLocations) > 0 and (current [0] - annotationLocations [-1] [0]) ** 2 + (np.log10 (current [1]) - np.log10 (annotationLocations [-1] [1])) ** 2 * 100 < minDist ** 2:
                    current = (current [0], current [1] * 10.0 ** (minDist / 10))
                annotationLocations.append (current)

                ax.annotate (z [2], xy = (z [0] [i], z [1] [i]), xycoords = "data", xytext = (current [0], current [1]), textcoords = "data", arrowprops=dict (arrowstyle = "->"))

        ax.set_yscale ("log")
        
        ax.grid (True)

        ax.set_xlabel ("Atomic Mass")
        if relativeIso is None:
            ax.set_ylabel ("Production Factor")
        elif elements:
            ax.set_xlabel ("Atomic Number")
            ax.set_ylabel ("$\\frac{X}{%s}/\\left.\\frac{X}{%s}\\right|_{\\mathrm{solar}}$" % (Isotope (relativeIso).getElementLabel (inMath = True), Isotope (relativeIso).getElementLabel (inMath = True)))
        else:
            ax.set_ylabel ("$\\frac{X}{%s}/\\frac{X}{%s}_{\\mathrm{solar}}$" % (Isotope (relativeIso).getLabel (inMath = True), Isotope (relativeIso).getLabel (inMath = True)))

        return lines