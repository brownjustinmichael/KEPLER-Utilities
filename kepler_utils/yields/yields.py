import os
import collections
import astropy.units as u
import numpy as np
from kepler_utils.records.dump import Isotope
from pandas import DataFrame
import pandas as pd
import re

def makeDeluxe (text, index_label = "", caption = ""):
    text = text.replace ("tabular", "deluxetable*").replace ("E+", "E$+$").replace ("E-", "E$-$")
    begin = text.split (r"\toprule") [0]
    mid = text.split (r"\toprule") [1]
    end = mid.split (r"\midrule") [1].replace (r"\bottomrule", r"\enddata")
    mid = mid.split (r"\midrule") [0]
    mid = " & ".join ([r"\colhead{" + (head.rstrip ("\n\\ ").lstrip ("\n\\ ") if head.rstrip ("\n\\ ").lstrip ("\n\\ ") != "{}" else index_label) + "}" for head in re.split (r"&", mid)])
    return begin + (r"\tablecaption{" + caption + "}\n" if caption != "" else "") + "\\tablehead{" + mid + "}\n\\startdata" + end

def fromKeplerYield (filename, directory, table = 1):
    f = open (directory + filename, "r")
    i = 1
    prevLine = "None"
    line = f.readline ()
    for j in range (2):
        while line.split () != [] or prevLine.split () != []:
            prevLine = line
            line = f.readline ()
            i += 1
        prevLine = line
        line = f.readline ()
        i += 1
    for j in range (table):
        while line.split () == []:
            prevLine = line
            line = f.readline ()
            i += 1

        while line.split () != []:
            prevLine = line
            line = f.readline ()
            i += 1
    return i

class YieldReader (object):
    def __init__ (self, directory = "yields/wh07/", extension = ".yield", masses = None, winds = True, explosions = True, keplerYield = True, totalYieldName = "yieldall", windYieldName = "yieldwind", expYieldName = None, isotopeName = "isotope", table = 1):
        results = collections.OrderedDict ()
        file_dict = {}
        self.models = []
        self.masses = []
        self.winds = winds
        self.explosions = explosions
        mass_file = open (directory + "/masses", "r")
        for line in mass_file:
            if line == "\n":
                continue
            line = line.rstrip ("\n").split (" ")
            if len (line) > 1:
                file_dict [float (line [0])] = line [1]
            self.masses.append (float (line [0]) * u.solMass)
            
        while len (file_dict) > 0:
            mass_min = min (file_dict) * u.solMass
            filename = file_dict.pop (min (file_dict))
            if keplerYield:
                i = fromKeplerYield (filename, directory, table)
            else:
                i = 0
            results [mass_min] = np.genfromtxt (directory + filename, skip_header = i, names = True, dtype = None)
            
        self.masses = u.Quantity (masses if masses is not None else self.masses)
        
        self.yields = DataFrame ()
        self.isotopes = []
        i = 0
        for result in results:
            while i < len (self.masses) and self.masses [i] < result * 0.999999:
                self.yields = self.yields.append ({}, ignore_index = True)
                self.models.append (False)
                i += 1
                continue
            if self.masses [i] > result * 1.000001:
                continue
            yieldDF = {}
            for row in results [result]:
                if row [isotopeName] == "total" or row [isotopeName] == b"total":
                    break
                isotope = Isotope (row [isotopeName]).string
                if isotope not in self.isotopes:
                    self.isotopes.append (isotope)
                yieldDF [isotope] = 0.0
                if self.winds and self.explosions and totalYieldName is not None:
                    yieldDF [isotope] += float (row [totalYieldName])
                else:
                    if self.winds:
                        yieldDF [isotope] += float (row [windYieldName])
                    if self.explosions:
                        if expYieldName is None:
                            yieldDF [isotope] += float (row [totalYieldName]) - float (row [windYieldName])
                        else:
                            yieldDF [isotope] += row [expYieldName]
            self.yields = self.yields.append (yieldDF, ignore_index = True)
            self.models.append (True)
            i += 1
        for j in range (len (self.masses) - i):
            self.yields = self.yields.append ({}, ignore_index = True)
            self.models.append (False)
        self.yields = self.yields.fillna (0.0)
        self.isotopes = [Isotope (iso) for iso in self.isotopes]
        self.isotopes.sort ()
        self.models = np.array (self.models)
            
    def get_yield (self, isotope, massArray = None, tolerance = 0.0001):
        if isinstance (isotope, Isotope):
            isotope = isotope.string
        if isotope not in self.yields:
            if massArray is None:
                return u.Quantity ([0.0] * len (self.yields), u.solMass)
            return u.Quantity ([0.0] * len (massArray), u.solMass)
                
        if massArray is None:
            return u.Quantity (self.yields [isotope], u.solMass)
        return u.Quantity (self.yields [isotope].iloc [massArray], u.solMass)
        
    def get_masses (self):
        return self.masses
        
    def __add__ (self, other):
        return CompositeYieldReader ((self.masses, other.masses), (self.isotopes, other.isotopes), (self.yields, other.yields))
        
    def __mul__ (self, scalar):
        new_yields = self.yields * scalar
        return CompositeYieldReader ((self.masses,), (self.isotopes,), (new_yields,))
        
    __rmul__ = __mul__
    
    def __div__ (self, scalar):
        return self * (1.0 / scalar)

    def __getitem__ (self, i):
        if isinstance (i, slice):
            return CompositeYieldReader ((self.masses [i],), (self.isotopes [i],), (self.yields [i],))
        return CompositeYieldReader ((self.masses [i:i+1],), (self.isotopes [i:i+1],), (self.yields [i:i+1],))

class CompositeYieldReader (YieldReader):
    def __init__ (self, masses, isotopes, yields, scale = 1.0):
        self.masses = u.Quantity (np.array (np.concatenate (masses)), masses [0].unit)
        self.isotopes = []
        for isotopeArray in isotopes:
            for iso in isotopeArray:
                if iso not in self.isotopes:
                    self.isotopes.append (iso)
        self.isotopes.sort ()
        self.yields = DataFrame ()
        for dataframe in yields:
            self.yields = self.yields.append (dataframe)
        self.yields = self.yields.fillna (0.0)
        self.yields *= scale
        