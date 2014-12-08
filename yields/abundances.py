import os

import periodictable
from records.dump import Isotope

__location__ = os.path.realpath (os.path.join (os.getcwd (), os.path.dirname (__file__)))

class Abundances (object):
    def __init__ (self, abundanceDict, normalization = None):
        self.total = None
        self.quantity = {}
        for isotope in abundanceDict:
            if isinstance (isotope, Isotope):
                isotope = isotope.string
            self.quantity [isotope] = abundanceDict [isotope]
            if self.total is None:
                self.total = abundanceDict [isotope]
            else:
                self.total += abundanceDict [isotope]
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
        
    def productionFactor (self, isotope):
        return self [isotope] / solar [isotope]
        
file = open (os.path.join (__location__, "sollo03.dat"))
abundanceDict = {}
for line in file:
    if line [0] != ';':
        words = line.split ()
        abundanceDict [words [0]] = float (words [1])
solar = Abundances (abundanceDict, 1.0)
