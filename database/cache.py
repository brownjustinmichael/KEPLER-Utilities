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

def calculate_tc (datadump):
    return datadump ['tn'] [0] * u.K

def calculate_dc (datadump):
    return datadump ['dn'] [0] * u.g / u.cm ** 3

def calculate_pc (datadump):
    return datadump ['pn'] [0] * u.dyne

def calculate_hshell_peak (datadump):
    return np.sum ((datadump ['sn'] * datadump ['xm']) [np.logical_and (datadump ['h1'] > 0.00001, datadump ['sn'] > 0.)])

def calculate_hshell_abun (datadump):
    return np.sum ((datadump ['h1'] * datadump ['sn']) [datadump ['h1'] > 0.00001]) / np.sum (datadump ['sn'] [datadump ['h1'] > 0.00001])

def calculate_hshell_rate (datadump):
    return np.sum ((datadump ['h1'] * datadump ['n14'] * datadump ['dn'] * datadump ['xm']) [datadump ['h1'] > 0.00001])

def calculate_max_sn (datadump):
    return np.max (datadump ['sn'])

def calculate_tasbsg (cnv_record):
    radii = u.Quantity ([model [-1] for model in cnv_record ['rncoord']])
    return np.sum (cnv_record ['dt'] [np.logical_and (radii > 2e12 * u.cm, radii < 8e12 * u.cm)]).to (u.year)