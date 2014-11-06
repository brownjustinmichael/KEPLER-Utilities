import os
import glob
import datetime
import inspect
import abc
import hashlib

import astropy.units as u
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base, declared_attr

import records.dump
import records.cnv

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
        return self.value * u.Unit (self.unit)
        
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
        def cache (self, session, cache_name, function, cache_data = True):
            # Create a string version of the function for code comparison
            code = ''.join (inspect.getsourcelines (function) [0])
    
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
            result = u.Quantity (function (self.get_data (cache_data)))
            self.cached.append (self.Cache (name = cache_name, value = result.value, unit = str (result.unit), code = code))
    
            # Commit the result and return it
            session.commit ()
            return result
        
        setattr (entrycls, "cache", cache)
        return entrycls
        
sim_tags = sqlalchemy.Table ('sim_tags', Base.metadata, sqlalchemy.Column ('sim_id', sqlalchemy.Integer, sqlalchemy.ForeignKey ('simulations.id')), sqlalchemy.Column ('tag_id', sqlalchemy.Integer, sqlalchemy.ForeignKey ('tags.id')))

@setup_parameters (records.dump.pfile)
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
    
    def __repr__ (self):
        return "<SimulationEntry(name='%s', id=%s)>" % (self.name, str (self.id))
        
    def get_state_dump (self, state):
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

# @sqlalchemy.event.listens_for(Session, 'after_flush')
# def delete_sim_orphans(session, ctx):
#     """
#     A simulation can only exist if it contains at least one file. After flushing the session, delete any orphaned simulations
#     """
#     session.query(SimulationEntry).\
#         filter(~SimulationEntry.dumpfiles.any()).\
#         filter(~SimulationEntry.cnvfiles.any()).\
#         delete(synchronize_session=False)

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
        return sqlalchemy.orm.relationship ("SimulationEntry", backref = cls.__name__.lower ().replace ("entry", "") + "s")
        
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
        entry, runid, name = cls.genFromFile (file)
        
        try:
            # print ("Checking for a corresponding simulation in the database...")
            if runid is None:
                simulation = session.query (SimulationEntry).filter_by (name = name).one ()
            else:
                simulation = session.query (SimulationEntry).filter_by (runid = runid).one ()
        except sqlalchemy.orm.exc.MultipleResultsFound:
            if log_info:
                print ("Multiple possible simulations in database, skipping")
            return
        except sqlalchemy.orm.exc.NoResultFound:
            if log_info:
                print ("No appropriate simulation found: creating")
            simulation = SimulationEntry (runid = runid, name = name, path = os.path.dirname (file), template_name = template_name, template_hash = template_hash)
            session.add (simulation)
        for tag in real_tags:
            if tag not in simulation.tags:
                simulation.tags.append (tag)
        try:
            simulation.get_state_dump (goal_state)
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
            cls.update_database (session, file, tags, template_name = template_name)
            
        # If any files match the glob string in any subdirectories, send them to update_database
        for root, dirs, files in os.walk (directory):
            for direct in dirs:
                direct = os.path.join (root, direct)
                if log_info:
                    print ("Searching in " + direct + " for " + os.path.join (direct, glob_string))
                for file in glob.glob (os.path.join (direct, glob_string)):
                    cls.update_database (session, file, tags, template_name = template_name)
        
        try:
            session.commit ()
        except:
            session.rollback ()
            raise TypeError ("Failed to commit")
        finally:
            session.close ()

@setup_parameters (records.dump.pfile)
@setup_parameters (records.dump.qfile)
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
        d = records.dump.Dump (file)
        entry = DumpFileEntry (file = file, date = datetime.datetime.fromtimestamp (os.path.getmtime(file)), timestep = d.ncyc, state = d.getState ())
        entry.set_parameters (d, records.dump.pfile)
        entry.set_parameters (d, records.dump.qfile)
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
        
    def get_data (self, cache = True):
        dataobject = records.dump.DataDump (self.file)
        if not cache:
            return dataobject
        if self.dataobject is None:
            self.dataobject = dataobject
        return self.dataobject

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
        
    def get_data (self, cache = True):
        dataobject = records.cnv.CNVFile (self.file)
        if not cache:
            return dataobject
        if self.dataobject == None:
            self.dataobject = dataobject
        return self.dataobject
        
    def addToSimulation (self, simulationEntry):
        self.simulation = simulationEntry

class Tag (Base):
    __tablename__ = 'tags'
    
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key = True)
    tag = sqlalchemy.Column (sqlalchemy.String, unique=True)
    
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

Base.metadata.create_all (engine)
