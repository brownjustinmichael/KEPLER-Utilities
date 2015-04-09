import kepler_utils.yields.janka as janka
import scipy.interpolate
import kepler_utils.yields.yields as yld
import kepler_utils.yields.integrator as integrator
import astropy.units as u
import numpy as np
import astropy.constants as consts
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

janka_dir = "/Users/justinbrown/Dropbox/Research/Stan/links14/"
janka_calibrations = ["w15.0", "w18.0", "w20.0", "n20.0"]

lower_sn_limit = 8
remnant_approximation = "si"

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
    integrated_baryonic_mass = (w_imf + janka_imf) (baryonic_masses, mask = baryonic_masses > 0.0 * u.solMass, mask_frequency = True)
    gravitational_masses = baryonic_masses * (1 - 3.0 / 5.0 * consts.G * baryonic_masses / (12 * u.km * consts.c ** 2))
    integrated_gravitational_mass = (w_imf + janka_imf) (gravitational_masses, mask = gravitational_masses > 0.0 * u.solMass, mask_frequency = True)
    print ("M_remnant", integrated_baryonic_mass, integrated_gravitational_mass)
    
    x = np.random.rand (1000000)
    x = ((1 - x) * max ((w_imf + janka_imf).masses) ** (-1.35) + x * min ((w_imf + janka_imf).masses) ** (-1.35)) ** (1. / -1.35)
    interpolator = interp1d ((w_imf + janka_imf).masses, baryonic_masses)
    inter = interpolator (x)
    
    fig, (ax1, ax2) = plt.subplots (2, 1)
    
    num, bins, patches = ax1.hist (inter [inter > 1.3], 50, label = calibration)
    ax1.set_xlabel ("Baryonic Mass")
    ax1.legend ()
    
    ax2.set_xlabel ("Gravitational Mass")
    
    interpolator = interp1d ((w_imf + janka_imf).masses, gravitational_masses)
    inter = interpolator (x)
    num, bins, patches = ax2.hist (inter [inter > 1.1], 50, color = "green")
    
    plt.tight_layout ()
    plt.savefig (calibration + "_neutron_dist.pdf")