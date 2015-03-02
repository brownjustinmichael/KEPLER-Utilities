import kepler_utils.yields.janka as janka
import scipy.interpolate
import kepler_utils.yields.yields as yld
import kepler_utils.yields.integrator as integrator
import astropy.units as u
import numpy as np
import astropy.constants as consts
import matplotlib.pyplot as plt



janka_dir = "/Users/justinbrown/Dropbox/Research/Stan/links14/"
janka_calibrations = ["w15.0", "w18.0", "w20.0", "n20.0"]

lower_sn_limit = 8
remnant_approximation = "si"

for calibration in janka_calibrations:
    jp = None
    jps = []
    janka_files = ["results_low_mass_%s.txt" % calibration, "results_%s.txt" % calibration, "results_heavy_%s.txt" % calibration]
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
    
    baryonic_masses [baryonic_masses > 2.5 * u.solMass] = 0. * u.solMass
    integrated_baryonic_mass = (w_imf + janka_imf) (baryonic_masses, mask = baryonic_masses > 0.0 * u.solMass, mask_frequency = True)
    gravitational_masses = baryonic_masses * (1 - 3.0 / 5.0 * consts.G * baryonic_masses / (12 * u.km * consts.c ** 2))
    integrated_gravitational_mass = (w_imf + janka_imf) (gravitational_masses, mask = gravitational_masses > 0.0 * u.solMass, mask_frequency = True)
    print ("M_remnant", integrated_baryonic_mass, integrated_gravitational_mass)
    
    numbers = (w_imf + janka_imf).freq [baryonic_masses > 0.0 * u.solMass]
    
    fig, (ax1, ax2) = plt.subplots (2, 1)
    
    tohist = [x for sublist in [[mass] * int (number) for mass, number in zip (baryonic_masses [baryonic_masses > 0.0 * u.solMass], numbers * 1000)] for x in sublist]

    num, bins, patches = ax1.hist (u.Quantity (tohist), 50, label = calibration)
    ax1.set_xlabel ("Baryonic Mass")
    ax1.legend ()
    
    tohist = [x for sublist in [[mass] * int (number) for mass, number in zip (gravitational_masses [gravitational_masses > 0.0 * u.solMass], numbers * 1000)] for x in sublist]
    ax2.set_xlabel ("Gravitational Mass")

    num, bins, patches = ax2.hist (u.Quantity (tohist), 50, color = "green")
    
    plt.tight_layout ()
    plt.savefig (calibration + "_neutron_dist.pdf")