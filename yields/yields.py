import os
import collections
import astropy.units as u
import numpy as np

class YieldReader (object):
    def __init__ (self, directory = "yields/wh07/", extension = ".yield"):
        results = collections.OrderedDict ()
        mass_file = {}
        for file in os.listdir (directory):
            name = file.replace (extension, "")
            for i in range (1, len (name)):
                if not name [i].isdigit ():
                    mod = name [i:]
                    if mod == "hi" or mod == "norm":
                        break
                    mass_file [(float (name [1:i]))] = file
                    break

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
        for result in results:
            for isotope, exp_yield, wind_yield in results [result]:
                if isotope.decode () not in self.wind_yields:
                    self.wind_yields [isotope.decode ()] = []
                self.wind_yields [isotope.decode ()].append (wind_yield)
        for isotope in self.wind_yields:
            self.wind_yields [isotope] = u.Quantity (self.wind_yields [isotope], u.solMass)
            
    def get_winds (self, index):
        return self.wind_yields [index]
        
    def get_masses (self):
        return self.masses