from sys import argv

from kepler_utils.yields.yields import YieldReader, Isotope, makeDeluxe
from kepler_utils.yields.integrator import IMFIntegrator, Integrator
import astropy.units as u
import numpy as np

from kepler_utils.plots.yields import YieldPlot

if len (argv) > 1:
    masses = argv [1:]
else:
    masses = [25]

w18yr = YieldReader (directory = "yields/yields_w18/", masses = np.concatenate ([np.arange (9.5, 12.0, 0.5), np.arange (12.0, 30.0, 0.1), np.array ([30., 31., 32., 33., 35., 40., 45., 50., 55., 60., 70., 80., 100., 120.])]) * u.solMass, winds = False)
whyr = YieldReader (explosions = False)

table = None
indices = []
for mass in masses:
    index = np.argmax (w18yr.masses >= float (mass) * u.solMass)
    ejecta = w18yr.yields.iloc [index : index + 1].reset_index () [[str (iso) for iso in w18yr.isotopes]]
    index = np.argmax (whyr.masses >= float (mass) * u.solMass)
    winds = whyr.yields.iloc [index : index + 1].reset_index () [[str (iso) for iso in whyr.isotopes]]
    
    for df in [ejecta, winds]:
        if table is None:
            table = df
        else:
            table = table.append (df, ignore_index = True)
    indices += ["Ejecta", "Winds"]

isos = [Isotope (col) for col in table.columns]
isos.sort ()
table = table [[str (iso) for iso in isos]]
table.columns = [iso.getLabel () for iso in isos]
table = table.fillna (0.0)
table.index = indices

print (makeDeluxe (table.transpose ().to_latex (float_format = '{:,.2E}'.format, escape = False), index_label = "Isotope", caption = r"Yields from $\sim$15 \Msun\ and $\sim$25 \Msun\ in \Msun"))