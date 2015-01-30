import os
import collections
import astropy.units as u
import numpy as np
from kepler_utils.records.dump import Isotope
from pandas import DataFrame

class YieldReader (object):
    def __init__ (self, directory = "yields/wh07/", extension = ".yield", masses = None, windType = False):
        results = collections.OrderedDict ()
        mass_file = {}
        self.models = []
        self.windType = windType
        for file in os.listdir (directory):
            name = file.split (".") [0]
            for i in file.split (".") [1:-1]:
                name += "." + i
            found = False
            try:
                for i in range (1, len (name)):
                    if (not name [i].isdigit () and not name [i] == "."):
                        mod = name [i:]
                        if mod == "hi" or mod == "norm":
                            found = True
                            break
                        mass_file [(float (name [1:i]))] = file
                        found = True
                        break
                if not found:
                    mass_file [(float (name [1:]))] = file
            except ValueError:
                mass_file [0.0] = file

        while len (mass_file) > 0:
            mass_min = min (mass_file) * u.solMass
            filename = mass_file.pop (min (mass_file))
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
            while line.split () != []:
                prevLine = line
                line = f.readline ()
                i += 1
            results [mass_min] = np.genfromtxt (directory + filename, skip_header = i, names = True, usecols = [0, 1, 6], dtype = ("|S8", "f8", "f8"))
            
        self.masses = u.Quantity ([result for result in results])
        if masses is not None:
            self.masses = u.Quantity (masses)
        
        self.wind_yields = DataFrame ()
        self.yields = DataFrame ()
        self.isotopes = []
        i = 0
        for result in results:
            while i < len (self.masses) and self.masses [i] < result * 0.99:
                self.wind_yields = self.wind_yields.append ({}, ignore_index = True)
                self.yields = self.yields.append ({}, ignore_index = True)
                self.models.append (False)
                i += 1
                continue
            if self.masses [i] > result * 1.01:
                continue
            windDF = {}
            yieldDF = {}
            for isotope, exp_yield, wind_yield in results [result]:
                if isotope.decode () == "total":
                    break
                if isotope.decode () not in self.isotopes:
                    self.isotopes.append (isotope.decode ())
                windDF [isotope.decode ()] = wind_yield
                yieldDF [isotope.decode ()] = (exp_yield - wind_yield)
            self.wind_yields = self.wind_yields.append (windDF, ignore_index = True)
            self.yields = self.yields.append (yieldDF, ignore_index = True)
            self.models.append (True)
            i += 1
        for j in range (len (self.masses) - i):
            self.wind_yields = self.wind_yields.append ({}, ignore_index = True)
            self.yields = self.yields.append ({}, ignore_index = True)
            self.models.append (False)
        self.yields = self.yields.fillna (0.0)
        self.wind_yields = self.wind_yields.fillna (0.0)
        self.isotopes = [Isotope (iso) for iso in self.isotopes]
        self.isotopes.sort ()
        self.models = np.array (self.models)
        
        if self.windType:
            self.yields = self.wind_yields
            
    def get_yield (self, isotope, massArray = None, tolerance = 0.0001):
        if isinstance (isotope, Isotope):
            isotope = isotope.string
        if massArray is None:
            return u.Quantity (self.yields [isotope], u.solMass)
        return u.Quantity (self.yields [isotope].iloc [massArray], u.solMass)
            
    def get_wind (self, isotope, massArray = None, tolerance = 0.0001):
        if isinstance (isotope, Isotope):
            isotope = isotope.string
        if massArray is None:
            return u.Quantity (self.wind_yields [isotope], u.solMass)
        return u.Quantity (self.wind_yields [isotope].iloc [massArray], u.solMass)
        
    def get_masses (self):
        return self.masses
        
    def __add__ (self, other):
        return CompositeYieldReader ((self.masses, other.masses), (self.isotopes, other.isotopes), (self.yields, other.yields), (self.wind_yields, other.wind_yields))

class CompositeYieldReader (YieldReader):
    def __init__ (self, masses, isotopes, yields, wind_yields):
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
        
        self.wind_yields = DataFrame ()
        for dataframe in wind_yields:
            self.wind_yields = self.wind_yields.append (dataframe)
        self.wind_yields = self.wind_yields.fillna (0.0)
        