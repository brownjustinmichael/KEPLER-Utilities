import os
import collections
import astropy.units as u
import numpy as np

class YieldReader (object):
    def __init__ (self, directory = "yields/wh07/", extension = ".yield"):
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
        
        self.wind_yields = {}
        self.yields = {}
        for result in results:
            for isotope, exp_yield, wind_yield in results [result]:
                if isotope.decode () not in self.wind_yields:
                    self.wind_yields [isotope.decode ()] = []
                    self.yields [isotope.decode ()] = []
                self.wind_yields [isotope.decode ()].append (wind_yield)
                self.yields [isotope.decode ()].append (exp_yield - wind_yield)
        for isotope in self.wind_yields:
            self.wind_yields [isotope] = u.Quantity (self.wind_yields [isotope], u.solMass)
            self.yields [isotope] = u.Quantity (self.yields [isotope], u.solMass)
            
    def get_yield (self, isotope, massArray = None, tolerance = 0.0001):
        if massArray is None:
            return self.yields [index]
        results = []
        for mass in massArray:
            found = False
            for i in range (len (self.masses)):
                if self.masses [i] > mass - tolerance * mass.unit:
                    if self.masses [i] < mass + tolerance * mass.unit:
                        results.append (self.yields [isotope] [i])
                        found = True
                        break
                    break
            if not found:
                results.append (0.0 * u.solMass)
        return u.Quantity (results)
            
    def get_wind (self, isotope, massArray = None, tolerance = 0.0001):
        if massArray is None:
            return self.wind_yields [index]
        results = []
        for mass in massArray:
            found = False
            for i in range (len (self.masses)):
                if self.masses [i] > mass - tolerance * mass.unit:
                    if self.masses [i] < mass + tolerance * mass.unit:
                        results.append (self.wind_yields [isotope] [i])
                        found = True
                        break
                    break
            if not found:
                results.append (0.0 * u.solMass)
        return u.Quantity (results)
        
    def get_masses (self):
        return self.masses