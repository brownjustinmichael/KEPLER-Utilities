#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
from sqlalchemy.sql import func

from kepler_utils.database.database import DumpFileEntry

class CorePlot (object):
    """
    ax is a MatPlotLib axis object
    query is a SQLAlchemy instance from the KEPLER database
        This works best if query has no columns. If there is no binning, query can have columns. If there is binning, query can only have aggregate columns.
    """
    
    def __init__ (self, query, ax = None):
        self.ax = ax
        if self.ax is None:
            fig, self.ax = plt.subplots (1, 1)
        self.query = query
        
    def plotScatter (self, xkey, ykey, skey = None, ckey = None, binkey = None, bins = 10, binMin = None, binMax = None, clsType = DumpFileEntry, errorbars = False, resize = False, xscale = 1.0, yscale = 1.0, **kwargs):
        query = self.query
        
        # If binning is intended (binkey is not None), group the data into bins
        if binkey is not None:
            # Calculate the min and max of the bin if binMin and binMax are not specified
            if binMax is None:
                query = clsType.cacheQuery (binkey, query, function = "max", label = "max")
            if binMin is None:
                query = clsType.cacheQuery (binkey, query, function = "min", label = "min")
            if query != self.query:
                res = query.first ()
                binMax = res.max if binMax is None else binMax
                binMin = res.min if binMin is None else binMin
                if binMax is None or binMin is None:
                    raise ValueError ("Nothing matches query")
            # Add the binning procedure to the query
            # TODO It might be nice if this could bin by non-numerical data too
            if isinstance (binkey, str):
                binkey = getattr (DumpFileEntry, binkey)
            query = self.query.add_column (func.floor (binkey / (binMax - binMin) * bins + binMin).label ("bin")).group_by ("bin").filter (binkey >= binMin, binkey <= binMax)
            
        # Add the keys to the query, using the cacheQuery method
        for key, name in zip ([xkey, ykey, skey, ckey], ["xkey", "ykey", "skey", "ckey"]):
            if key is not None:
                query = clsType.cacheQuery (key, query, function = None if binkey is None else "avg", label = name)
                if binkey is not None and errorbars:
                    query = clsType.cacheQuery (key, query, function = "stddev", label = name + "std")
                    
        # Execute the query
        results = query.all ()
        
        # Plot the data
        if skey is not None:
            kwargs ["s"] = np.array ([res.skey for res in results])
            if resize:
                minimum = np.min (kwargs ["s"])
                maximum = np.max (kwargs ["s"])
                if minimum != maximum:
                    kwargs ["s"] = (kwargs ["s"] - minimum) / (maximum - minimum) * 100 + 20
                else:
                    kwargs ["s"] = 200
        if ckey is not None:
            kwargs ["c"] = [res.ckey for res in results]
        
        if binkey is None or not errorbars:
            return self.ax.scatter ([res.xkey * xscale for res in results], [res.ykey * yscale for res in results], **kwargs)
        else:
            return self.ax.errorbar ([res.xkey * xscale for res in results], [res.ykey * yscale for res in results], xerr = [res.xkeystd * xscale for res in results], yerr = [res.ykeystd * yscale for res in results], linestyle = "None", **kwargs)
