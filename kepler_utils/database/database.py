import os
import glob
import datetime
import inspect
import abc
import hashlib
import types

import numpy as np
import astropy.units as u
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base, declared_attr

from kepler_utils.records.dump import pfile, qfile, Dump, DataDump
from kepler_utils.records.cnv import CNVFile

# Setup the declarative base class for inheritance
Base = declarative_base ()

# Start the database engine with a database called dumpfiles.db in the current directory
engine = sqlalchemy.create_engine ('sqlite:///' + os.path.join (os.getcwd (), 'dumpfiles.db'), echo = False)

# Create a session class that's bound to the above engine
Session = sqlalchemy.orm.sessionmaker (bind = engine)

# A conversion from string to sqlalchemy types
types = {'float' : sqlalchemy.Float, 'integer' : sqlalchemy.Integer}

class setup_parameters (object):
    """
    This decorator is designed to take a parameter file with the following tab-separated format:
    ParamNumber	ParamName	(float or integer)	defaultValue	unit	Description with lots of words	versionNumberAdded
    
    It will add an appropriate column to a SQL table for each parameter in the file
    
    It will also add (or add to) a list called parameters such that it contains all the imported parameters
    """
    def __init__ (self, file):
        """
        Save the parameter file name
        """
        self.file = file
    
    def __call__ (self, cls):
        """
        Load the parameters from the saved file name
        """
        # Open the file and create the parameter list
        parameter_file = open (self.file, 'r')
        parameters = []
        
        # Iterate through the lines in the file
        for line in parameter_file:
            words = line.split ('\t')
            if words [0] [0] == '#':
                # Skip commented lines
                continue
                
            # Add the new parameter to the parameters list
            parameters.append (words [1])
            
            # Add the parameter attribute to the class
            setattr (cls, words [1], sqlalchemy.Column (words [1], types [words [2]]))
            
        try:
            # Try adding the parameter list to the old one in the class
            old_params = getattr (cls, 'parameters')
            setattr (cls, 'parameters', old_params + parameters)
        except AttributeError:
            # If there is no old parameter list, add a new one
            setattr (cls, 'parameters', parameters)
            
        return cls

class Cache (object):
    """
    A mixin object for a cache entry object to inherit, providing the relevant columns and accessors
    """
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key = True)
    name = sqlalchemy.Column (sqlalchemy.String)
    code = sqlalchemy.Column (sqlalchemy.String)
    
    value = sqlalchemy.Column (sqlalchemy.Float)
    unit = sqlalchemy.Column (sqlalchemy.String)

    def get_value (self):
        """
        Return the current cache value (with appropriate AstroPy units)
        """
        if self.value is not None:
            return self.value * u.Unit (self.unit)
        else:
            return np.nan * u.Unit (self.unit)
        
    @classmethod
    def CacheFactory (cls, name, entrycls, tablename, **kwargs):
        """
        Create a cache object inside of cls that is itself a SQLAlchemy entry object
        """
        newclass = type (name, (cls, Base), {"__tablename__": tablename, \
            "file": sqlalchemy.Column (sqlalchemy.String, sqlalchemy.ForeignKey (entrycls.file)), \
            "entry": sqlalchemy.orm.relationship (entrycls, **kwargs), \
            "__table_args__": (sqlalchemy.UniqueConstraint('name', 'file', name='uix_1'),)})
        setattr (entrycls, "Cache", newclass)
        return newclass
    
    @classmethod
    def caching (cls, entrycls):
        """
        This decorator function is designed to add a caching ability to any subclass of Base
        """
        # Create the cache entry and table using CacheFactory
        entrycls.Cache = cls.CacheFactory (entrycls.__name__ + "Cache", entrycls, entrycls.__name__.lower ().replace ("entry", "") + "cache", backref = sqlalchemy.orm.backref ('cached'))
    
        # Define the caching method for cls, which will check whether the value has been cached previously
        def cache (self, session, cache_name, function, **kwargs):
            # Create a string version of the function for code comparison
            try:
                code = ''.join (inspect.getsourcelines (function) [0])
            except TypeError:
                code = 'unable to read'
    
            # Check if the cached function already exists in the current cache
            for cache in self.cached:
                if cache.name == cache_name and cache.code == code:
                    # If it exists, return the value
                    return cache.get_value ()
                elif cache.name == cache_name and cache.code != code:
                    # If the code has changed, delete the cache
                    session.delete (cache)
                    session.commit ()
            
            # Evaluate the cached function
            kwargs ["cache"] = kwargs.pop ("cache", True)
            if function == None:
                raise ValueError ("Unable to find value in database.")
            result = u.Quantity (function (self.get_data (**kwargs)))
            self.cached.append (self.Cache (name = cache_name, value = result.value, unit = str (result.unit), code = code))
    
            # Commit the result and return it
            session.commit ()
            return result
        
        setattr (entrycls, "cache", cache)
        return entrycls
        
sim_tags = sqlalchemy.Table ('sim_tags', Base.metadata, sqlalchemy.Column ('sim_id', sqlalchemy.Integer, sqlalchemy.ForeignKey ('simulations.id')), sqlalchemy.Column ('tag_id', sqlalchemy.Integer, sqlalchemy.ForeignKey ('tags.id')))

@setup_parameters (pfile)
class SimulationEntry (Base):
    """
    A database entry for simulations, which will each hold a number of file entry objects
    """
    __tablename__ = 'simulations'
    
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key=True)
    runid = sqlalchemy.Column(sqlalchemy.LargeBinary, nullable = True)
    name = sqlalchemy.Column (sqlalchemy.String)
    path = sqlalchemy.Column (sqlalchemy.String)
    loaded = sqlalchemy.Column (sqlalchemy.Boolean, default = False)
    complete = sqlalchemy.Column (sqlalchemy.Boolean, default = False)
    template_name = sqlalchemy.Column (sqlalchemy.String, nullable = True)
    template_hash = sqlalchemy.Column (sqlalchemy.String, nullable = True)
    
    tags = sqlalchemy.orm.relationship ('Tag', secondary = sim_tags, backref = 'sims')
    
    __table_args__ = (sqlalchemy.UniqueConstraint ('path', 'name', name = '_path_name_uc'),)
    
    def __repr__ (self):
        return "<SimulationEntry(name='%s', id=%s)>" % (self.name, str (self.id))
        
    def getStateDump (self, state):
        """
        Get a dump file in the simulation with the given string state, e.g. "presn", "hdep", "2000", etc.
        """
        for dump in self.dumpfiles:
            if dump.state == state:
                return dump
        raise IndexError ("State %s not found in simulation" % state)
        
    def copy_parameters (self, entry):
        """
        Copy the parameters from entry into the simulation
        """
        for param in self.parameters:
            setattr (self, param, getattr (entry, param))
            
    def check_parameters (self, entry):
        """
        Compare the simulation parameters to those from entry. If any are different, turn them to None in the simulation entry.
        
        The reasoning behind this is that the simulation has a set of parameters that define it, but because of the nature of KEPLER, these parameters can change during the course of a run. Any parameters that change would not be good parameters to query for the simulation. (Sadly, any parameters can change at any point, so we can't just remove these from the parameter list.)
        """
        for param in self.parameters:
            if getattr (self, param) != getattr (entry, param):
                setattr (self, param, None)
                
    def tag (self, session, tagName):
        tag = Tag.get (session, tagName)
        if tag not in self.tags:
            self.tags.append (tag)
        session.commit ()

class FileEntry (object):
    """
    This base mixin is a superclass for a file entry object, and it contains all the pieces that will be needed for any generic file entry
    """
    file = sqlalchemy.Column (sqlalchemy.String, primary_key=True)
    date = sqlalchemy.Column (sqlalchemy.DateTime)
    
    @declared_attr
    def simulation_id(cls):
        return sqlalchemy.Column ('simulation_id', sqlalchemy.ForeignKey('simulations.id'))
    
    @declared_attr
    def simulation (cls):
        return sqlalchemy.orm.relationship ("SimulationEntry", backref = sqlalchemy.orm.backref (cls.__name__.lower ().replace ("entry", "") + "s", cascade="all, delete-orphan"))
        
    def __init__ (self, **kwargs):
        super (FileEntry, self).__init__ (**kwargs)
        self.dataobject = None
        
    @sqlalchemy.orm.reconstructor
    def init_on_load (self):
        self.dataobject = None
    
    @abc.abstractmethod
    def addToSimulation (self, simulationEntry):
        pass
        
    @classmethod
    @abc.abstractmethod
    def genFromFile (cls, file):
        pass
        
    @classmethod
    def update_database (cls, session, file_name, tags = None, template_name = None, goal_state = 'presn', log_info = False):
        """
        Check whether a file should be added to the database, and do so
        """
        real_tags = []
        if isinstance (tags, str):
            tags = [tags]
        if tags is not None:
            for tag in tags:
                real_tags.append (Tag.get (session, tag))
        
        template_hash = None
        if template_name is not None:
            template_hash = hashlib.md5 (open (template_name).read ().encode ()).hexdigest ()
        
        # Get the absolute path to the file
        file = os.path.abspath (file_name)
        if log_info:
            print ("Checking for file " + file + " in database")
        
        try:
            # Check that there is not a newer or equivalent file in the database already
            oldentry = session.query (cls).filter_by (file = file).one ()
            if oldentry.date >= datetime.datetime.fromtimestamp (os.path.getmtime(file)):
                if log_info:
                    print ("Newer or equivalent file in database: skipping")
                for tag in real_tags:
                    if tag not in oldentry.simulation.tags:
                        oldentry.simulation.tags.append (tag)
                # oldentry.simulation.template_name = template_name
                # oldentry.simulation.template_hash = template_hash
                return
            else:
                # If there's an older file, be sure to delete it
                session.delete (oldentry)
        except sqlalchemy.orm.exc.NoResultFound:
            if log_info:
                print ("File not found in database: adding")
            
        # If we've arrived here, we'll need to add a new entry to the database.
        # Load the data from the file
        try:
            entry, runid, name = cls.genFromFile (file)
        except Exception as e:
            print (e)
            raise e
        
        try:
            # print ("Checking for a corresponding simulation in the database...")
            simulation = session.query (SimulationEntry).filter_by (name = name).filter_by (path = os.path.dirname (file)).one ()
            if runid is not None and simulation.runid != runid:
                print ("Warning: Simulation runid does not match run, resetting")
                # raise RuntimeError ("Simulation at %s with name %s does not match file %s" % (os.path.dirname (file), name, file_name))
                session.delete (simulation)
                session.commit ()
                raise sqlalchemy.orm.exc.NoResultFound
        except sqlalchemy.orm.exc.NoResultFound:
            if log_info:
                print ("No appropriate simulation found: creating")
            simulation = SimulationEntry (runid = runid, name = name, path = os.path.dirname (file), template_name = template_name, template_hash = template_hash)
            session.add (simulation)
        for tag in real_tags:
            if tag not in simulation.tags:
                simulation.tags.append (tag)
        try:
            simulation.getStateDump (goal_state)
            simulation.complete = True
        except IndexError:
            pass
        entry.addToSimulation (simulation)
        session.add (entry)
        
    @classmethod
    def scan_for_updates (cls, directory, glob_string = '*', tags = None, template_name = None, log_info = False):
        """
        Scan the directory for files matching glob_string and add them to the database
        """
        # Create a session
        session = Session ()
        
        # Go through the database, searching for entries that match the current class and check whether the corresponding file has been deleted
        for entry in session.query (cls).all ():
            if not os.path.isfile (entry.file):
                # If the file has been deleted, removed the entry from the database
                if log_info:
                    print ("File " + entry.file + " has been deleted: removing from database")
                session.delete (entry)
        
        # If any files match the glob string in the current directory, send them to update_database
        for file in glob.glob (os.path.join (directory, glob_string)):
            cls.update_database (session, file, tags, template_name = template_name, log_info = log_info)
            
        # If any files match the glob string in any subdirectories, send them to update_database
        for root, dirs, files in os.walk (directory):
            for direct in dirs:
                direct = os.path.join (root, direct)
                if log_info:
                    print ("Searching in " + direct + " for " + os.path.join (direct, glob_string))
                for file in glob.glob (os.path.join (direct, glob_string)):
                    cls.update_database (session, file, tags, template_name = template_name, log_info = log_info)
        
        try:
            session.commit ()
        except:
            session.rollback ()
            raise TypeError ("Failed to commit")
        finally:
            session.close ()

@setup_parameters (pfile)
@setup_parameters (qfile)
@Cache.caching
class DumpFileEntry (FileEntry, Base):
    __tablename__ = 'dumpfiles'
    
    timestep = sqlalchemy.Column (sqlalchemy.Integer)
    state = sqlalchemy.Column (sqlalchemy.String)
        
    def __repr__(self):
        return "<Dump (name='%s', timestep='%s', state=%s)>" % (self.simulation.name, self.timestep, self.state)
    
    def __init__ (self, **kwargs):
        super (DumpFileEntry, self).__init__ (**kwargs)
       
    def addToSimulation (self, simulationEntry):
        self.simulation = simulationEntry
        if not simulationEntry.loaded:
            simulationEntry.copy_parameters (self)
            simulationEntry.loaded = True
        simulationEntry.check_parameters (self)
        
    @classmethod
    def genFromFile (cls, file):
        d = Dump (file)
        entry = DumpFileEntry (file = file, date = datetime.datetime.fromtimestamp (os.path.getmtime(file)), timestep = d.ncyc, state = d.getState ())
        entry.set_parameters (d, pfile)
        entry.set_parameters (d, qfile)
        return entry, d.getRunId (), d.namep
        
    def set_parameters (self, dump, file):
        parameter_file = open (file, 'r')
        for line in parameter_file:
            words = line.split ('\t')
            if words [0] [0] == '#':
                return
            try:
                setattr (self, words [1], dump.parameters [words [1]].value)
            except KeyError as e:
                pass
        
    def get_data (self, cache = True, **kwargs):
        if self.dataobject is not None:
            return self.dataobject
        dataobject = DataDump (self.file, **kwargs)
        if not cache:
            return dataobject
        if self.dataobject is None:
            self.dataobject = dataobject
        return self.dataobject
        
    def flush (self):
        self.dataobject = None

@Cache.caching
class CNVFileEntry (FileEntry, Base):
    __tablename__ = 'cnvfiles'
    
    name = sqlalchemy.Column (sqlalchemy.String)
    
    def __init__ (self, **kwargs):
        super (CNVFileEntry, self).__init__ (**kwargs)
        
    def __repr__(self):
       return "<CNV File (name='%s')>" % (self.simulation.name)
        
    @classmethod
    def genFromFile (cls, file):
        name = os.path.basename (os.path.splitext (file) [0])
        entry = CNVFileEntry (file = file, name = name, date = datetime.datetime.fromtimestamp (os.path.getmtime(file)))
        return entry, None, name
        
    def get_data (self, cache = True, **kwargs):
        if self.dataobject is not None:
            return self.dataobject
        dataobject = CNVFile (self.file, **kwargs)
        if not cache:
            return dataobject
        if self.dataobject == None:
            self.dataobject = dataobject
        return self.dataobject
        
    def flush (self):
        self.dataobject = None
        
    def addToSimulation (self, simulationEntry):
        self.simulation = simulationEntry

class Tag (Base):
    __tablename__ = 'tags'
    
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key = True)
    tag = sqlalchemy.Column (sqlalchemy.String, unique=True)
    
    def __repr__ (self):
        return "<Tag %s>" % self.tag
    
    @classmethod
    def get (cls, session, tag):
        try:
            return session.query (cls).filter (cls.tag == tag).one ()
        except sqlalchemy.orm.exc.NoResultFound:
            tag_obj = cls (tag = tag)
            session.add (tag_obj)
            session.commit ()
            return tag_obj

def basicQuery (session):
    return session.query (SimulationEntry, DumpFileEntry).join (DumpFileEntry)

def cache (session, sims, funcs, states = ("presn")):
    dumps = [[sim.getStateDump (state) for state in states] for sim in sims]
    
    results = {}
    for name in funcs:
        results [name] = []
        for run in dumps:
            results [name].append ([])
        
    for i, run in enumerate (dumps):
        for state, dump in zip (states, run):
            for name in funcs:
                results [name] [i].append (dump.cache (session, name + "_" + state, funcs [name]))
            dump.flush ()
        for name in results:
            results [name] [i] = u.Quantity (results [name] [i])
            
    for name in results:
        results [name] = u.Quantity (results [name])
            
    return results
            
def cnv_cache (session, sims, funcs):
    if isinstance (sims, sqlalchemy.orm.query.Query):
        sims = sims.all ()
        
    if not isinstance (sims, SimulationEntry):
        sims = [sim [0] for sim in sims]
        
    cnvs = [sim.cnvfiles [0] for sim in sims]
    
    results = {}
    for name in funcs:
        results [name] = []
        
    for cnv in cnvs:
        for name in funcs:
            results [name].append (cnv.cache (session, name, funcs [name], verbose = True))
        cnv.flush ()
        
    for name in results:
        results [name] = u.Quantity (results [name])
        
    return results

Base.metadata.create_all (engine)
