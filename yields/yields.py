import os
import collections
import astropy.units as u
import numpy as np
from records.dump import Isotope
from pandas import DataFrame


class YieldReader (object):
    def __init__ (self, directory = "yields/wh07/", extension = ".yield", masses = None):
        results = collections.OrderedDict ()
        mass_file = {}
        for file in os.listdir (directory):
            name = file.split (".") [0]
            for i in file.split (".") [1:-1]:
                name += "." + i
            found = False
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

        while len (mass_file) > 0:
            mass_min = min (mass_file) * u.solMass
            filename = mass_file.pop (min (mass_file))
            f = open (directory + filename, "r")
            for i in range (550):
                f.readline ()
            i = 550
            line = f.readline ()
            while line.split () != []:
                line = f.readline ()
                i += 1
            results [mass_min] = np.genfromtxt (directory + filename, skip_header = i, names = True, usecols = [0, 1, 6], dtype = ("|S8", "f8", "f8"))
            
        self.masses = u.Quantity ([result for result in results])
        if masses is not None:
            self.masses = masses
        
        self.wind_yields = DataFrame ()
        self.yields = DataFrame ()
        self.isotopes = []
        i = 0
        for result in results:
            while i < len (self.masses) and self.masses [i] <= result * 0.99:
                self.wind_yields = self.wind_yields.append ({}, ignore_index = True)
                self.yields = self.yields.append ({}, ignore_index = True)
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
            i += 1
        for j in range (len (self.masses) - i):
            self.wind_yields = self.wind_yields.append ({}, ignore_index = True)
            self.yields = self.yields.append ({}, ignore_index = True)
        self.yields = self.yields.fillna (0.0)
        self.wind_yields = self.wind_yields.fillna (0.0)
        self.isotopes = [Isotope (iso) for iso in self.isotopes]
        self.isotopes.sort ()
            
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