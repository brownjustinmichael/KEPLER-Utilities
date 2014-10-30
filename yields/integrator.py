import numpy as np
import astropy.units as u

class IMFIntegrator (object):
    def __init__ (self, initial_masses, alpha = -2.35):
        self.masses = u.Quantity (initial_masses [:])
        self.mdiff = np.zeros (len (self.masses))
        self.mdiff [:-1] = np.diff (initial_masses)
        self.freq = (self.masses**alpha * self.mdiff).value
        self.freq /= np.sum (self.freq)
        
    def __call__ (self, array, mask = None, mask_frequency = False):
        if mask is None:
            mask = np.ones (len (self.masses), dtype = bool)
        freq = self.freq
        if mask_frequency:
            freq /= np.sum (self.freq [mask])
        return np.sum (array [mask] * freq [mask])