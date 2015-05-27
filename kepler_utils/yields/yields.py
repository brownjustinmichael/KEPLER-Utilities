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

def fromKeplerYield (filename, table = 1):
    f = open (filename, "r")
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
    def __init__ (self, yields = None, scale = 1.0):
        self.masses = [] if yields is None else [ind [0] for ind in yields.index]
        self.isotopes = [] if yields is None else [Isotope (iso) for iso in yields.columns]
        self.isotopes.sort ()
        self.yields = DataFrame () if yields is None else yields
        self.yields = self.yields.fillna (0.0)
        self.yields *= scale

    @classmethod
    def from_file (cls, filename, mass, **kwargs):
        self = cls ()
        self.add_file (filename, mass, **kwargs)
        return self

    @classmethod
    def from_directory (cls, directory = "yields/wh07/", mass_file = "masses", **kwargs):
        self = cls ()
        mass_file = open (directory + "/" + mass_file, "r")
        for line in mass_file:
            if line == "\n":
                continue
            line = line.rstrip ("\n").split (" ")
            try:
                self.add_file (directory + "/" + line [1], float (line [0]) * u.solMass, **kwargs)
            except IndexError:
                self.yields = self.yields.append (DataFrame ([{"mass": float (line [0]) * u.solMass, "file": directory + "/"}]).set_index (["mass", "file"]))
                self.masses.append (float (line [0]) * u.solMass)
        return self

    @classmethod
    def combine (cls, yield_readers):
        self = cls ()
        self.masses = u.Quantity (np.array (np.concatenate ([yr.masses for yr in yield_readers])))
        for yr in yield_readers:
            isotopeArray = yr.isotopes
            for iso in isotopeArray:
                if iso not in self.isotopes:
                    self.isotopes.append (iso)
        self.isotopes.sort ()
        for yr in yield_readers:
            dataframe = yr.yields
            self.yields = self.yields.append (dataframe)
        self.yields = self.yields.fillna (0.0)
        return self

    def add_file (self, filename, mass, winds = True, explosions = True, keplerYield = True, totalYieldName = "yieldall", windYieldName = "yieldwind", expYieldName = None, isotopeName = "isotope", table = 1):
        self.masses.append (mass)
        if keplerYield:
            i = fromKeplerYield (filename, table)
        else:
            i = 0
        result = np.genfromtxt (filename, skip_header = i, names = True, dtype = None)

        yieldDF = {}
        yieldDF ["mass"] = mass
        yieldDF ["file"] = filename
        for row in result:
            if row [isotopeName] == "total" or row [isotopeName] == b"total":
                break
            isotope = Isotope (row [isotopeName])
            if isotope not in self.isotopes:
                self.isotopes.append (isotope)
            yieldDF [isotope.string] = 0.0
            if winds and explosions and totalYieldName is not None:
                yieldDF [isotope.string] += float (row [totalYieldName])
            else:
                if winds:
                    yieldDF [isotope.string] += float (row [windYieldName])
                if explosions:
                    if expYieldName is None:
                        yieldDF [isotope.string] += float (row [totalYieldName]) - float (row [windYieldName])
                    else:
                        yieldDF [isotope.string] += row [expYieldName]
        self.yields = self.yields.append (DataFrame ([yieldDF]).set_index (["mass", "file"]))
        self.yields = self.yields.fillna (0.0)
        self.isotopes.sort ()
  
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

    def get_keys (self):
        return [i for i in self.yields.index]
        
    def __add__ (self, other):
        return YieldReader (self.yields.add (other.yields, fill_value = 0.0))
        
    def __mul__ (self, scalar):
        return YieldReader (self.yields, scalar)
        
    __rmul__ = __mul__
    
    def __div__ (self, scalar):
        return self * (1.0 / scalar)

    def __getitem__ (self, i):
        if isinstance (i, slice):
            return YieldReader (self.yields [i])
        return YieldReader (self.yields [i:i+1])

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
        