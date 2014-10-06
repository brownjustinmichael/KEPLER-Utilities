import numpy as np
import abc

import astropy.units as u
        
def calculate_he_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['h1'] > 0.1)]
    
def calculate_co_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['he4'] > 0.2)]

def calculate_ne_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['c12'] > 0.1)]

def calculate_si_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['o16'] > datadump ['si28'])]

def calculate_fe_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['si28'] > datadump ['fe56'])]

def calculate_tasbsg (cnv_record):
    radii = u.Quantity ([model [-1] for model in cnv_record ['rncoord']])
    return np.sum (cnv_record ['dt'] [np.logical_and (radii > 2e12 * u.cm, radii < 8e12 * u.cm)]).to (u.year)