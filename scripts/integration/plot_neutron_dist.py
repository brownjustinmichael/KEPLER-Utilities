import kepler_utils.yields.janka as janka
import scipy.interpolate
import kepler_utils.yields.yields as yld
import kepler_utils.yields.integrator as integrator
import astropy.units as u
import numpy as np
import astropy.constants as consts
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from scipy.interpolate import interp1d

janka_dir = "/Users/justinbrown/Dropbox/Research/Stan/Nucleosynthesis/links14/"
janka_calibrations = ["w15.0", "w18.0", "w20.0", "n20.0"]

lower_sn_limit = 8
remnant_approximation = "si"

def bary_to_grav (x):
    return (x * (1 - 3.0 / 5.0 * consts.G * x / (12 * u.km * consts.c ** 2))).to (u.solMass)
    
def grav_to_bary (x):
    return (10 * u.km * consts.c ** 2 / consts.G * (1 - np.sqrt (1 - consts.G * x / 5. / consts.c ** 2 / u.km))).to (u.solMass)

fontP = FontProperties()
fontP.set_size('small')

observed = np.array ([1.3381, 1.2489, 1.3332, 1.3452, 1.32, 1.24, 1.365, 1.248, 1.4414, 1.3867, 1.358, 1.354, 1.438, 1.27])

for calibration in janka_calibrations:
    jp = None
    jps = []
    janka_files = ["results_%s_revised.txt" % calibration]
    for f in janka_files:
        jps.append (janka.JankaParser.readFrom (janka_dir + f))
        if jp is None:
            jp = jps [-1]
        else:
            jp += jps [-1]
            
    janka_imf = integrator.IMFIntegrator (jp [":M"])
    
    print ("Calibration: ", calibration)
    
    wdata = np.genfromtxt ("wcores.dat", names = True)
    w_imf = integrator.IMFIntegrator (wdata ["Minit"])
    w_remnants = wdata [remnant_approximation] * u.solMass

    baryonic_masses = u.Quantity (list (w_remnants) + list (jp [":with_fallback:M_mass_cut_after_fb"]))
    # baryonic_masses = jp [":with_fallback:M_mass_cut_after_fb"]
    
    baryonic_masses [baryonic_masses > 2.5 * u.solMass] = 0. * u.solMass
    gravitational_masses = bary_to_grav (baryonic_masses)
    
    x = np.random.rand (1000000)
    x = ((1 - x) * max ((w_imf + janka_imf).masses) ** (-1.35) + x * min ((w_imf + janka_imf).masses) ** (-1.35)) ** (1. / -1.35)
    
    fig, ax = plt.subplots (1, 1)

    ax.set_xlabel ("Gravitational Mass")
    ax.set_ylabel ("PDF")
    axtop = ax.twiny ()
    
    interpolator = interp1d ((w_imf + janka_imf).masses, gravitational_masses)
    inters = [interpolator (x [np.logical_and (low < x, x < high)]) for low, high in zip ((9, 10, 13, 15, 18), (10, 13, 15, 18, 120))]
    num, bins, patches = ax.hist (observed, 9, normed = True, color = "gray", label = "Schwab 2010", linestyle = "dashed")
    num, bins, patches = ax.hist ([inter [inter > 1.18] for inter in inters], 40, normed = True, stacked = True, label = ["%d $M_{\odot}$ < M < %d $M_{\odot}$" % (low, high) for low, high in zip ((9, 10, 13, 15, 18), (10, 13, 15, 18, 120))])
    
    ticks = np.arange (1, 3, 0.1) * u.solMass
    
    axtop.set_xticks (bary_to_grav (ticks).value)
    axtop.set_xlim (ax.get_xlim ())
    axtop.set_xticklabels (ticks.value)
    axtop.set_xlabel ("Baryonic Mass")
    
    plt.tight_layout ()
    ax.legend (prop = fontP)
    plt.savefig (calibration + "_neutron_dist.pdf")