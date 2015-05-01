import math

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.mlab
import matplotlib.collections
import matplotlib.colors
import matplotlib.scale
import matplotlib.transforms

from kepler_utils.records.cnv import CNVFile
import kepler_utils.plots.shiftlog

msun = 1.988435e33

converter = matplotlib.colors.ColorConverter ()

def returnTrue (value):
    return True
    
def returnValue (value, model):
    return value
    
def returnZero (value, model):
    return 0

def energyBin (value, model):
    if value > 0:
        return float (10. ** (value + model ['mingain']))
    elif value < 0:
        return float (-10. ** (-value + model ['minloss']))
    else:
        return float (0.0)

transparent_cmap = matplotlib.colors.LinearSegmentedColormap ("transparent", {'red':((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), 'green':((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), 'blue':((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), 'alpha':((0.0, 0.0001, 0.0001), (1.0, 0.0001, 0.0001))})

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.legend import Legend
import matplotlib.patches
import matplotlib.legend_handler

class PolyCollectionHandler (matplotlib.legend_handler.HandlerBase):
    def __init__(self, patch_func=None, **kw):
        matplotlib.legend_handler.HandlerBase.__init__(self, **kw)

        self._patch_func = patch_func
        
    def create_artists (self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):
        x0, y0 = xdescent, ydescent
        width, height = width, height
        linewidth = orig_handle.get_linewidth () [0] if len (orig_handle.get_linewidth ()) != 0 else None
        patch = matplotlib.patches.Rectangle([x0, y0], width, height, hatch=orig_handle.get_hatch (), lw=linewidth,
                                   transform=trans)
        cmap = orig_handle.get_cmap ()
        patch.set_facecolor (orig_handle.get_cmap () (0.5))
        patch.set_edgecolor (orig_handle.get_edgecolor () [0])
        return [patch]

Legend.update_default_handler_map({matplotlib.collections.PolyCollection: PolyCollectionHandler()})
            
class KippenhahnPlot(object):
    """This object takes a matplotlib axis and attaches a Kippenhahn plot to it associated with the given CNVFile"""
    def __init__(self, axis, cnv_file, useModels = False):
        super(KippenhahnPlot, self).__init__()
        self.axis = axis
        self.cnv_file = cnv_file
        
        self.times = (u.Quantity ([model ['timesec'] for model in self.cnv_file], 's').to ('year')).value
        
        if (not useModels):
            self.x = self.times
        else:
            self.x = list (range (len (self.times)))
            
        self.xmin = min (self.x)
        self.xmax = max (self.x)
        
        self.mmin = self.cnv_file [0] ['xmcoord'] [0] / msun
        self.mmax = self.cnv_file [0] ['xmcoord'] [-1] / msun
        
        self.tend = self.xmax
        i = -1
        while (self.tend == self.xmax):
            self.tend += self.x [-1] - self.x [-1 + i]
            i -= 1
            
        # Set the plot labels
        if (not useModels):
            self.axis.set_xlabel ("t (years)")
        else:
            self.axis.set_xlabel ("Model Number")
        self.axis.set_ylabel ("m (solar masses)")
        
    def get_tmax (self):
        return self.xmax
        
    def plotMax (self, index, scale = 1.0, **kwargs):
        return self.axis.plot (self.x, scale * np.array ([max (model [index]) for model in self.cnv_file]), **kwargs)
        
    def _distance (self, value1, value2, logspace = True):
        if logspace:
            return abs (np.log10 ((self.tend - value2) / (self.tend - value1)))
        else:
            return abs (value2 - value1)
        
    def plotRectangles (self, indices, values, top_indexed, logspace = True, points = 400, extent = None, set_conditions = (returnTrue,), get_value = None, **kwargs):
        hatch = kwargs.pop ('hatch', "")
        colors = converter.to_rgba_array (kwargs.pop ('color', ""))
        label = kwargs.pop ('label', None)
        
        if extent == None:
            extent = range (0, len (self.cnv_file))
            
        if get_value == None:
            get_value = []
            for i in range (len (set_conditions)):
                get_value.append (returnValue)
                
        total = self._distance (self.x [list (extent) [0]], self.x [list (extent) [-1]], logspace)
        
        polygons = []
        z = []
        for condition in set_conditions:
            polygons.append ([])
            z.append ([])
        previous_i = None
        for i in extent:
            if previous_i == None:
                previous_i = i
                continue
            model = self.cnv_file [i]
            left = self.x [previous_i]
            right = self.x [i]
            if points is not None and (self._distance (left, right, logspace) < total / points):
                continue
            previous_i = i
            if (left == right):
                continue
            for j in range (len (model [indices])):
                index = model [indices] [j] - 1
                if top_indexed:
                    top = model ['xmcoord'] [index] / msun
                    if j == 0:
                        bottom = model ['xmcoord'] [0] / msun
                    else:
                        bottom = model ['xmcoord'] [model [indices] [j - 1]] / msun
                if not top_indexed:
                    bottom = model ['xmcoord'] [index] / msun
                    if j == len (model [indices]) - 1:
                        top = model ['xmcoord'] [model ['ncoord'] - 1] / msun
                    else:
                        top = model ['xmcoord'] [model [indices] [j + 1]] / msun
                for k in range (len (set_conditions)):
                    if set_conditions [k] (model [values] [j]):
                        polygons [k].append ([(left, bottom), (left, top), (right, top), (right, bottom)])
                        z [k].append (get_value [k] (model [values] [j], model))
        collections = []
        converter
        for i in range (len (polygons)):
            collections.append (matplotlib.collections.PolyCollection (polygons [i], array = np.array (z [i]), **kwargs))
            self.axis.add_collection (collections [-1], autolim = True)
            if isinstance (hatch, str):
                collections [-1].set_hatch (hatch)
            elif len (hatch) > 0:
                collections [-1].set_hatch (hatch [i % len (hatch)])
            if len (colors) > 0:
                collections [-1].set_color (colors [i % len (colors)])
            if isinstance (label, str):
                collections [-1].set_label (label)
            elif label is not None and len (label) > 0:
                collections [-1].set_label (label [i % len (label)])
        return collections
        
    def plotEnergy (self, points = 400, logspace = True, extent = None, **kwargs):
        # Calculate the maximum log10 value of energy gain/loss throughout the star to evenly balance the SymLog colorbar 
        vmax = max ([max (model ['nuc']) if len (model ['nuc']) != 0 else 0 for model in self.cnv_file])

        cmap = kwargs.pop ('cmap', plt.get_cmap ('coolwarm_r'))
        linewidth = kwargs.pop ('linewidth', 0.0)
        norm = kwargs.pop ('norm', matplotlib.colors.SymLogNorm (0.1, vmin = -10.**vmax, vmax = 10.**vmax))
        edgecolor = kwargs.pop ('edgecolor', 'white')

        # Generate the energy plot from indices 'inuc' and values 'nuc'
        return self.plotRectangles ('inuc', 'nuc', False, get_value = (energyBin,), cmap = cmap, linewidth = linewidth, norm = norm, edgecolor = edgecolor, logspace = logspace, points = points, extent = extent)
        
    def addEnergyColorBar (self, energy):
        vmax = max ([max (model ['nuc']) if len (model ['nuc']) != 0 else 0 for model in self.cnv_file])

        # Add the color bar to the plot
        cb = plt.colorbar (energy)
        # Specify the ticks for the colorbar and make their string representation with powers of 10
        ticks = [10. ** i for i in range (-1, vmax + 1, vmax // 5 if vmax >= 5 else 1)]
        ticks += [-10. ** i for i in range (-1, vmax + 1, vmax // 5 if vmax >= 5 else 1)]
        cb.set_ticks(ticks)
        if (matplotlib.rcParams ['text.usetex']):
            cb.set_ticklabels([r'${10^{%d}}$' % (t) for t in range (-1, vmax + 1, vmax // 5 if vmax >= 5 else 1)] + [r'${-10^{%d}}$' % (t) for t in range (-1, vmax + 1, vmax // 5 if vmax >= 5 else 1)])
        else:
            cb.set_ticklabels(['$\\mathdefault{10^{%d}}$' % (t) for t in range (-1, vmax + 1, vmax // 5 if vmax >= 5 else 1)] + ['$\\mathdefault{-10^{%d}}$' % (t) for t in range (-1, vmax + 1, vmax // 5 if vmax >= 5 else 1)])
    
        cb.set_label ("energy generation loss/gain (erg/g/s)")
        return cb
        
    def plotConvection (self, points = 400, logspace = True, extent = None, **kwargs):
        # Generate the mixing plot
        conditions = (lambda x: x == 'C', lambda x: x == 'O', lambda x: x == 'S', lambda x: x == 'T')
        returns = ((lambda x, model: 0),) * 4
        
        edgecolor = kwargs.pop ('edgecolor', 'white')
        cmap = kwargs.pop ('cmap', transparent_cmap)
        linewidth = kwargs.pop ('linewidth', 0.0)
        hatch = kwargs.pop ('hatch', ["////", "\\\\\\", "xxxx", "."])
        color = kwargs.pop ('color', ('black', 'gray', 'r', 'b'))
        label = kwargs.pop ('label', ("Convective", "Overshooting", "Semi-Convective", "Thermohaline"))
        
        return self.plotRectangles ('iconv', 'yzip', True, set_conditions = conditions, get_value = returns, edgecolor = edgecolor, cmap = cmap, linewidth = linewidth, logspace = logspace, points = points, hatch = hatch, color = color, label = label, extent = extent)
        
def jTDPlot (cnv_record, logspace = True, extent = None, points = 400):
    if not isinstance (cnv_record, CNVFile):
        cnv_record = CNVFile (cnv_record)

    # Generate a matplotlib figure with one subplot
    fig = plt.figure (figsize = (18,10))
    ax = plt.subplot (111)

    # Initialize the KippenhahnPlot object with the given axis and record file
    kippplot = KippenhahnPlot (ax, cnv_record)

    energy, = kippplot.plotEnergy (points = points, logspace = logspace, extent = extent)
    cb = kippplot.addEnergyColorBar (energy)

    kippplot.plotConvection (points = points, logspace = logspace, extent = extent)

    # Generate the outer edge of the star
    mass, = kippplot.plotMax ('xmcoord', 1.0 / msun, color = 'black', label = "Total Mass")

    # Set the axis to be logarithmically scaled around the end of the star's life, with exponents increasing backward
    if (logspace):
        ax.set_xscale ('shiftlog', base = 10, zero = kippplot.tend, sign = -1.0)

    # Set the plot limits and show the grid and legend
    if extent == None:
        ax.set_xlim (kippplot.xmin, kippplot.xmax - 10.**-5)
    else:
        ax.set_xlim (kippplot.x [list (extent) [0]], kippplot.x [list (extent) [-1]])
    ax.set_ylim (kippplot.mmin, kippplot.mmax)
    plt.grid ()
    plt.legend ()
    
    radii = u.Quantity ([model [-1] for model in cnv_record ['rncoord']])
    
    return fig, ax