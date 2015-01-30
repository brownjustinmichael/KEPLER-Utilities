import numpy as np
import astropy.units as u
from .abundances as import Abundances

class IMFIntegrator (object):
    def __init__ (self, initial_masses, alpha = -2.35, numberArray = None, normalize = 1.0):
        self.masses = u.Quantity (initial_masses [:])
        if numberArray is None:
            # If we're here, a numberArray was not specified. Use a Salpeter IMF and check that the mass range is monotonic TODO allow for nonmonotonic mass arrays
            for i in range (1, len (initial_masses)):
                assert (initial_masses [i - 1] <= initial_masses [i])
            mdiff = np.zeros (len (self.masses))
            mdiff [:-1] = np.diff (initial_masses)
            self.freq = (self.masses**alpha * mdiff).value
        else:
            assert (len (initial_masses) == len (numberArray))
            self.freq = [i for i in numberArray]
        if normalize is not None:
            self.freq /= np.sum (self.freq)
            self.freq *= normalize
        
    def __call__ (self, array, mask = None, mask_frequency = False):
        if mask is None:
            mask = np.ones (len (self.masses), dtype = bool)
        freq = np.copy (self.freq)
        if mask_frequency:
            freq /= np.sum (self.freq [mask])
        return np.sum (array [mask] * freq [mask])
        
    def __add__ (self, other):
        masses = np.concatenate ((self.masses, other.masses))
        freq = np.concatenate ((self.freq, other.freq))
        return IMFIntegrator (masses, None, freq)
        
    def __mul__ (self, scalar):
        return IMFIntegrator (self.masses, None, self.freq, scalar)
        
    __rmul__ = __mul__
        
    def __div__ (self, scalar):
        return IMFIntegrator (self.masses, None, self.freq, 1.0 / scalar)
        
    def getAbundances (self, yieldReader, imfUpperLimit = None, imfLowerLimit = None):
        yd = {}
        mask = None
        if imfUpperLimit is not None:
            mask = yieldReader.get_masses () < imfUpperLimit * 1.01
        if imfLowerLimit is not None:
            if mask is None:
                mask = yieldReader.get_masses () > imfLowerLimit * 0.99
            else:
                mask = np.logical_and (mask, yieldReader.get_masses () > imfLowerLimit * 0.99)
           
        for isotope in yieldReader.isotopes:
            value = self (yieldReader.get_yield (isotope), mask = mask)
            yd [isotope.string] = value
           
        return Abundances (yd)