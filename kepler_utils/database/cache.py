import numpy as np
import abc
import math

import astropy.units as u

def sn_hshell (datadump):
    return np.max (datadump ["sn"] [datadump ["h1"] > 0.01])

def sn_width (datadump):
    masses = datadump ["mass coordinate"] [np.where (np.logical_and (datadump ["sn"] > 0.1 * np.max (datadump ["sn"]), datadump ["h1"] > 0.01))]
    if len (masses) == 0:
        return datadump ["mass coordinate"] [0] - datadump ["mass coordinate"] [0]
    else:
        return masses [-1] - masses [0]

def rad_habun (datadump):
    return np.sum ((datadump ["xm"] * datadump ["h1"]) [np.argmax (datadump ["sn"] [100:]) + 100:np.argmax (datadump ["sn"] [100:]) + 200]) / (datadump ["mass coordinate"] [np.argmax (datadump ["sn"] [100:]) + 200] - datadump ["mass coordinate"] [np.argmax (datadump ["sn"] [100:]) + 100])

def compactness (datadump):
    return 2.5 / (datadump ['rn'] [np.argmax (datadump ['mass coordinate'] > 2.5 * u.solMass)] / (1000 * u.km)).to (1)
    
def co_ratio (datadump):
    return np.sum ((datadump ['xm'] * datadump ['c12']) [0:np.argmax (datadump ['he4'] > 0.2)]) / np.sum ((datadump ['xm'] * datadump ['o16']) [0:np.argmax (datadump ['he4'] > 0.2)])

def sum_iso (isotope, datadump):
    return np.sum (datadump ['xm'] * datadump [isotope])

def h_contained (datadump):
    return datadump ["mass coordinate"] [np.argmax (np.cumsum (datadump ['xm'] * datadump ["h1"]) > 1.0 * u.solMass)]
    
def he_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['h1'] > 0.1)]

def co_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['he4'] > 0.2)]

def ne_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['c12'] > 0.1)]

def si_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['o16'] > datadump ['si28'])]

def fe_core (datadump):
    return datadump ['mass coordinate'] [np.argmax (datadump ['si28'] > datadump ['fe56'])]

def tc (datadump):
    return datadump ['tn'] [0]

def dc (datadump):
    return datadump ['dn'] [0]

def pc (datadump):
    return datadump ['pn'] [0]

def h1c (datadump):
    return datadump ['h1'] [0]

def he4c (datadump):
    return datadump ['he4'] [0]

def hshell_peak (datadump):
    return np.sum ((datadump ['sn'] * datadump ['xm']) [np.logical_and (datadump ['h1'] > 0.00001, datadump ['sn'] > 0.)])

def hshell_abun (datadump):
    return np.sum ((datadump ['h1'] * datadump ['sn']) [datadump ['h1'] > 0.00001]) / np.sum (datadump ['sn'] [datadump ['h1'] > 0.00001])

def hshell_habun (datadump):
    return datadump ["h1"] [np.argmax (datadump ["xln"] > 0.97 * np.max (datadump ["xln"]))]

def hshell_rate (datadump):
    return np.sum ((datadump ['h1'] * datadump ['n14'] * datadump ['dn'] * datadump ['xm']) [datadump ['h1'] > 0.00001])

def max_sn (datadump):
    return np.max (datadump ['sn'])

def tasbsg (cnv_record):
    radii = u.Quantity ([model [-1] for model in cnv_record ['rncoord']])
    return np.sum (cnv_record ['dt'] [np.logical_and (radii > 2e12 * u.cm, radii < 8e12 * u.cm)]).to (u.year)

def ledd_shell (datadump):
    return datadump ["ledd"] [np.argmax (datadump ["xln"] > np.max (datadump ["xln"]))]

def core_overshoot (datadump):
    i = np.argmax ((datadump ["icon"] == "osht"))
    if not np.all (np.logical_or (datadump ["icon"] == "conv", datadump ["icon"] == "semi") [:i]):
        return 0.0
    j = np.min ([np.argmax ((datadump ["icon"] != "osht") [i:]), i + np.argmax (datadump ["difi"] [i:] < datadump ["difi"] [i] / np.e)]) + 1
    print (j - i)
    return (datadump ["rn"] [j] - datadump ["rn"] [i]) / (-0.5 * (datadump ["pn"] [i + 1] + datadump ["pn"] [i + 1]) * (datadump ["rn"] [i + 1] - datadump ["rn"] [i]) / (datadump ["pn"] [i + 1] - datadump ["pn"] [i]))

def env_overshoot (datadump):
    i = np.argmax ((datadump ["icon"] == "conv") [::-1])
    i = np.argmax ((datadump ["icon"] == "osht") [i::-1])
    j = np.max ([np.argmax ((datadump ["icon"] != "osht") [i::-1]), np.argmax (i - (datadump ["difi"] [i::-1] < datadump ["difi"] [i] / np.e) [i::-1])]) - 1
    return (datadump ["rn"] [i] - datadump ["rn"] [j]) / (-0.5 * (datadump ["pn"] [i + 1] + datadump ["pn"] [i + 1]) * (datadump ["rn"] [i + 1] - datadump ["rn"] [i]) / (datadump ["pn"] [i + 1] - datadump ["pn"] [i]))
    
def L_ratio (datadump):
    return np.max (datadump ["xln"]) / datadump ["xln"] [-1]

def fasRSG (cnvdata):
    h1 = cnvdata ["abun_cnv"] [:,1]
    total = (cnvdata [-1] ["timesec"] - cnvdata [np.argmax (h1 < 0.001)] ["timesec"]) * u.s
    radii = u.Quantity ([i [-1] for i in cnvdata ["rncoord"]])
    return np.sum (cnvdata ["dt"] [radii > 1.e13 * u.cm]) / total

def fasBL (cnvdata):
    radii = u.Quantity ([i [-1] for i in cnvdata ["rncoord"]])
    for i, r in enumerate (radii):
        if r > 1.e13 * u.cm:
            break

    h1 = cnvdata ["abun_cnv"] [:,1]
    total = cnvdata [-1] ["timesec"] - cnvdata [np.argmax (h1 < 0.001)] ["timesec"]

    j = i + 100 + np.argmax (radii [i + 100:] < 1.e13 * u.cm)

    if j == i + 100:
        return 0.0

    for i, r in enumerate (radii [j:]):
        if r > 1.e13 * u.cm:
            break

    return (cnvdata [i + j] ["timesec"] - cnvdata [j] ["timesec"]) / total  
