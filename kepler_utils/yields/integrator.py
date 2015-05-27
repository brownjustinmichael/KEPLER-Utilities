import numpy as np
import astropy.units as u
from .abundances import Abundances, solar

class Integrator (object):
    def __init__ (self, massArray, numberArray, normalize = 1.0, sort_by = None):
        self.sort_by = massArray if sort_by is None else sort_by
        self.freq = np.array ([x for y, x in sorted (zip (self.sort_by, numberArray), key = lambda pair: pair [0])])
        self.freq /= np.sum (self.freq)
        self.normalize = np.ones (len (self.freq))
        self.masses = u.Quantity (massArray)
        try:
            self.normalize *= float (normalize)
        except TypeError:
            self.normalize [:] = [x for y, x in sorted (zip (self.sort_by, normalize), key = lambda pair: pair [0])]
        self.sort_by.sort ()
                    
    def __call__ (self, array, mask = None, mask_frequency = False, imfLowerLimit = None, imfUpperLimit = None):
        if imfUpperLimit is not None:
            maskLower = np.logical_or (self.masses < imfUpperLimit * 1.01, self.masses == 0.0 * imfUpperLimit.unit)
            if mask is not None:
                mask = np.logical_and (mask, maskLower)
            else:
                mask = maskLower
        if imfLowerLimit is not None:
            maskUpper = np.logical_or (self.masses > imfLowerLimit * 0.99, self.masses == 0.0 * imfLowerLimit.unit)
            if mask is not None:
                mask = np.logical_and (mask, maskUpper)
            else:
                mask = maskUpper

        if mask is None:
            mask = np.ones (len (self), dtype = bool)
        freq = np.copy (self.freq) * self.normalize
        if mask_frequency:
            freq /= np.sum (self.freq [mask])
        return np.sum (array [mask] * freq [mask])
        
    def __add__ (self, other):
        freq = np.concatenate ((self.freq, other.freq))
        normalize = np.concatenate ((self.normalize, other.normalize))
        sort_by = np.concatenate ((self.sort_by, other.sort_by))
        return Integrator (freq, normalize, sort_by = sort_by)
        
    def __mul__ (self, scalar):
        return Integrator (self.freq, scalar * self.normalize, sort_by = self.sort_by)
        
    __rmul__ = __mul__
        
    def __div__ (self, scalar):
        return Integrator (self.freq, self.normalize / scalar, sort_by = self.sort_by)
        
    def __len__ (self):
        return len (self.freq)
        
    def getAbundances (self, yieldReader, imfUpperLimit = None, imfLowerLimit = None):
        yd = {}
        for isotope in yieldReader.isotopes:
            value = self (yieldReader.get_yield (isotope), imfLowerLimit = imfLowerLimit, imfUpperLimit = imfUpperLimit)
            yd [isotope.string] = value
           
        return Abundances (yd)
        
    def makeSolar (self, isotope1, isotope2, other, yr, other_yr = None):
        if other_yr is None:
            other_yr = yr
        ab1 = self.getAbundances (yr)
        ab2 = other.getAbundances (other_yr)
        
        # print (other_yr.get_yield ("fe56"))
        # print ((yr + other_yr).get_yield ("fe56"))
        
        x = ab1 [isotope1] / solar [isotope1] * ab1.total
        z = ab1 [isotope2] / solar [isotope2] * ab1.total
        
        y = ab2 [isotope1] / solar [isotope1] * ab2.total
        w = ab2 [isotope2] / solar [isotope2] * ab2.total
        
        a = float ((w - y) / (x + w - y - z))
        
        return a, a * self + (1.0 - a) * other, yr + other_yr

class IMFIntegrator (Integrator):
    def __init__ (self, massArray, alpha = -2.35, **kwargs):
        self.alpha = alpha
        self.original_normalize = kwargs.get ("normalize", 1.0)
        massArray = u.Quantity (massArray)
        
        if alpha == -1.0:
            raise ValueError ("Take your own logarithms.")
        
        values = (massArray ** (alpha + 1.0) / (alpha + 1.0)).value

        if len (massArray) == 1:
            numberArray = np.array ([1.0])
        else:
            numberArray = np.zeros (len (massArray))
            numberArray [:-1] = np.diff (values)
            numberArray [numberArray < 0.0] = 0.0
        
        # We assume here that the coverage is continuous
        
        total = ((max (massArray) ** (alpha + 1.0) - min (massArray) ** (alpha + 1.0)) / (alpha + 1.0) * kwargs.get ("normalize", 1.0)).value
        print (total)
        print ((1.0 + (((np.sum (numberArray) - total) / total) if total != 0.0 else 0.0)))

        super (IMFIntegrator, self).__init__ (u.Quantity (massArray [:]), numberArray, normalize = (1.0 + (((np.sum (numberArray) - total) / total) if total != 0.0 else 0.0)))
        
    @classmethod
    def from_yieldreader (cls, yieldreader, **kwargs):
        if "sort_by" not in kwargs:
            kwargs ["sort_by"] = yieldreader.get_keys ()
        return cls (yieldreader.masses, **kwargs)

    def __add__ (self, other):
        if isinstance (other, IMFIntegrator):
            if other.alpha != self.alpha or other.original_normalize != self.original_normalize:
                raise TypeError ("Can't add IMFs with different alphas or original normalization")
            masses = u.Quantity (np.array (np.concatenate ((self.masses, other.masses.to (self.masses.unit)))), self.masses.unit)
            sort_by = list (self.sort_by) + list (other.sort_by)
            return IMFIntegrator (masses, self.alpha, normalize = self.original_normalize, sort_by = sort_by)
        else:
            return super (IMFIntegrator, self).__add__ (other)
