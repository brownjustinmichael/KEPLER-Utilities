#!/usr/bin/env python

from sqlalchemy.sql import func

from kepler_utils.database.database import DumpFileEntry

class CorePlot (object):
    """ax is a MatPlotLib axis object, query is a SQLAlchemy instance from the KEPLER database"""
    
    def __init__ (self, ax, query):
        self.ax = ax
        self.query = query
        
    def plotScatter (self, xkey, ykey, skey = None, ckey = None, binkey = None, bins = 10, binMin = None, binMax = None, clsType = DumpFileEntry, **kwargs):
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
            # Add the binning procedure to the query
            # TODO It might be nice if this could bin by non-numerical data too
            if isinstance (binkey, str):
                binkey = getattr (DumpFileEntry, binkey)
            query = self.query.add_column (func.floor (binkey / (binMax - binMin) * bins + binMin).label ("bin")).group_by ("bin")
            
        # Add the keys to the query, using the cacheQuery method
        query = clsType.cacheQuery (xkey, query, function = None if binkey is None else "avg", label = "xkey")
        query = clsType.cacheQuery (ykey, query, function = None if binkey is None else "avg", label = "ykey")
        if skey is not None:
            query = clsType.cacheQuery (skey, query, function = None if binkey is None else "avg", label = "skey")
        if skey is not None:
            query = clsType.cacheQuery (ckey, query, function = None if binkey is None else "avg", label = "ckey")
                    
        # Execute the query
        results = query.all ()
        
        # Plot the data
        if skey is not None:
            kwargs ["s"] = [res.skey for res in results]
        if ckey is not None:
            kwargs ["c"] = [res.ckey for res in results]
        
        self.ax.scatter ([res.xkey for res in results], [res.ykey for res in results], **kwargs)
