import os

import periodictable
from kepler_utils.records.dump import Isotope

__location__ = os.path.realpath (os.path.join (os.getcwd (), os.path.dirname (__file__)))

class Abundances (object):
    def __init__ (self, abundanceDict, normalization = None):
        self.total = None
        self.quantity = {}
        self.isotopes = {}
        for isotope in abundanceDict:
            if isinstance (isotope, Isotope):
                iso = isotope.string
            else:
                iso = isotope
            
            isotope = Isotope (iso)
            
            if isotope.z not in self.isotopes:
                self.isotopes [isotope.z] = []
            self.isotopes [iso] = self.isotopes [isotope.z]
            self.isotopes [iso].append (iso)
                
            self.quantity [iso] = abundanceDict [iso]
            if self.total is None:
                try:
                    self.total = abundanceDict [iso].copy ()
                except AttributeError:
                    self.total = abundanceDict [iso]
            else:
                self.total += abundanceDict [iso]
        if normalization is not None:
            self.total = normalization
    
    def __getitem__ (self, isotope):
        if isinstance (isotope, Isotope):
            isotope = isotope.string
        return (self.quantity [isotope] / self.total)
        
    def __contains__ (self, isotope):
        if isinstance (isotope, Isotope):
            isotope = isotope.string
        return isotope in self.quantity
        
    def productionFactor (self, isotope, element = False, relativeIso = None):
        if isinstance (isotope, Isotope):
            isotope = isotope.string
        
        if element:
            result = sum ([self [iso] for iso in self.isotopes [isotope]]) / sum ([solar [iso] for iso in self.isotopes [isotope]])
        else:
            result = self [isotope] / solar [isotope]
        
        if relativeIso is None:
            return result
        else:
            return result / self.productionFactor (relativeIso, element = element)
        
file = open (os.path.join (__location__, "sollo03.dat"))
abundanceDict = {}
for line in file:
    if line [0] != ';':
        words = line.split ()
        abundanceDict [words [0]] = float (words [1])
solar = Abundances (abundanceDict, 1.0)
