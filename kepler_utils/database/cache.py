import numpy as np
import abc
import math

import astropy.units as u

def calculate_compactness_2_5 (datadump):
    return 2.5 / (datadump ['rn'] [np.argmax (datadump ['mass coordinate'] > 2.5 * u.solMass)] / (1000 * u.km)).to (1)
    
def calculate_co_ratio_in_core (datadump):
    return np.sum ((datadump ['xm'] * datadump ['c12']) [0:np.argmax (datadump ['he4'] > 0.2)]) / np.sum ((datadump ['xm'] * datadump ['o16']) [0:np.argmax (datadump ['he4'] > 0.2)])
        
def sum_iso (isotope, datadump):
    return np.sum (datadump ['xm'] * datadump [isotope])
    
def calculate_h_contained (datadump):
    return datadump ["mass coordinate"] [np.argmax (np.cumsum (datadump ['xm'] * datadump ["h1"]) > 1.0 * u.solMass)]
    
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
    
def core_overshoot (datadump):
    regions = []
    current = None
    previous = 0
    for i in range (len (datadump ['icon'])):
        if current == None:
            current = datadump ['icon'] [i] == "osht"
        if current != (datadump ['icon'] [i] == "osht"):
            if current:
                regions.append (list (range (previous, i + 1)))
            previous = i
            current = (datadump ['icon'] [i] == "osht")
    
    q = -1
    i = 0
    while i == 0 and q < len (regions):
        q += 1
        i = np.argmax (datadump ['difi'] [regions [0]] < datadump ['difi'] [regions [0] [0]] / math.e)
    ri = regions [0] [i]
    frac = 0.0
    for j in range (i):
        frac -= 0.5 * (datadump ['rn'] [ri + j + 1] - datadump ['rn'] [ri + j - 1]) / (0.5 * (datadump ['pn'] [ri + j] + datadump ['pn'] [ri + j - 1]) / (datadump ['pn'] [ri + j] - datadump ['pn'] [ri + j - 1]) * (datadump ['rn'] [ri + j] - datadump ['rn'] [ri + j - 1]))
    return frac

def env_overshoot (datadump):
    regions = []
    current = None
    previous = 0
    for i in range (len (datadump ['icon'])):
        if current == None:
            current = datadump ['icon'] [i] == "osht"
        if current != (datadump ['icon'] [i] == "osht"):
            if current:
                regions.append (list (range (previous, i + 1)))
            previous = i
            current = (datadump ['icon'] [i] == "osht")
    
    q = 0
    i = 0
    while (i == 0 and -(q) < len (regions)):
        q -= 1
        i = np.argmin ((datadump ['difi'] [regions [q]] < datadump ['difi'] [regions [q] [-1]] / math.e))
    ri = regions [q] [i]
    frac = 0.0
    for j in range (0, len (regions [q]) - i):
        frac -= 0.5 * (datadump ['rn'] [ri + j + 1] - datadump ['rn'] [ri + j - 1]) / (0.5 * (datadump ['pn'] [ri + j] + datadump ['pn'] [ri + j - 1]) / (datadump ['pn'] [ri + j] - datadump ['pn'] [ri + j - 1]) * (datadump ['rn'] [ri + j] - datadump ['rn'] [ri + j - 1]))
    return frac