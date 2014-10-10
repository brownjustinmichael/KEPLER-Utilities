import os
import glob
import datetime
import inspect

import astropy.units as u
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base, DeferredReflection

import records.dump
import records.cnv

__location__ = os.path.realpath (os.path.join (os.getcwd (), os.path.dirname (__file__)))

engine = sqlalchemy.create_engine ('sqlite:///' + os.path.join (os.getcwd (), 'dumpfiles.db'), echo = False)
# Base = declarative_base ()
Base = declarative_base (DeferredReflection)
Session = sqlalchemy.orm.sessionmaker (bind = engine)

class setup_parameters:
    def __init__ (self, file):
        self.file = file
    
    def __call__ (self, cls):
        if not cls.loaded:
            parameter_file = open (self.file, 'r')
            parameters = []
            for line in parameter_file:
                words = line.split ('\t')
                if words [0] [0] == '#':
                    continue
                parameters.append (words [1])
                if words [2] == 'float':
                    setattr (cls, words [1], sqlalchemy.Column (words [1], sqlalchemy.Float))
                elif words [2] == 'integer':
                    setattr (cls, words [1], sqlalchemy.Column (words [1], sqlalchemy.Integer))
                else:
                    assert (False)
            try:
                old_params = getattr (cls, 'parameters')
                setattr (cls, 'parameters', old_params + parameters)
            except AttributeError:
                setattr (cls, 'parameters', parameters)
        return cls

@setup_parameters (records.dump.pfile)
class SimulationEntry (Base):
    __tablename__ = 'simulations'
    loaded = False
    
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key=True)
    runid = sqlalchemy.Column(sqlalchemy.LargeBinary)
    name = sqlalchemy.Column (sqlalchemy.String)
    path = sqlalchemy.Column (sqlalchemy.String)
    
    def __repr__ (self):
        return "<SimulationEntry(name='%s', id=%s)>" % (self.name, str (self.id))
        
    def get_state_dump (self, state):
        for dump in self.dumpfiles:
            if dump.state == state:
                return dump
        raise IndexError ("State %s not found in simulation" % state)
        
    def copy_parameters (self, entry):
        for param in self.parameters:
            setattr (self, param, getattr (entry, param))
            
    def check_parameters (self, entry):
        for param in self.parameters:
            if getattr (self, param) != getattr (entry, param):
                setattr (self, param, None)

@setup_parameters (records.dump.pfile)
@setup_parameters (records.dump.qfile)
class DumpFileEntry (Base):
    loaded = False

    __tablename__ = 'dumpfiles'
    file = sqlalchemy.Column (sqlalchemy.String, primary_key=True)

    runid = sqlalchemy.Column (sqlalchemy.LargeBinary, sqlalchemy.ForeignKey ('simulations.runid'))
    date = sqlalchemy.Column (sqlalchemy.DateTime)
    timestep = sqlalchemy.Column (sqlalchemy.Integer)
    state = sqlalchemy.Column (sqlalchemy.String)
    
    simulation = sqlalchemy.orm.relationship ("SimulationEntry", backref = sqlalchemy.orm.backref ('dumpfiles'))

    @sqlalchemy.orm.reconstructor
    def init_on_load (self):
        self.dataobject = None

    def __repr__(self):
       return "<Dump (name='%s', timestep='%s', state=%s)>" % (
                            self.simulation.name, self.timestep, self.state)
       
    @classmethod
    def scan_for_updates (cls, directory, glob_string = '*'):
        glob_string += '#*'
        session = Session ()
        for entry in session.query (DumpFileEntry).all ():
            if not os.path.isfile (entry.file):
                print ("File " + entry.file + " has been deleted: removing from database")
                session.delete (entry)
                session.commit ()
        for simentry in session.query (SimulationEntry).all ():
            if len (simentry.dumpfiles) == 0 and simentry.cnvfile is None:
                print ("Simulation " + simentry.name + " is empty: removing from database")
                session.delete (simentry)
                session.commit ()
        for file in glob.glob (os.path.join (directory, glob_string)):
            cls.update_database (session, file)
        for root, dirs, files in os.walk (directory):
            for direct in dirs:
                direct = os.path.join (root, direct)
                print ("Searching in " + direct + " for " + os.path.join (direct, glob_string))
                for file in glob.glob (os.path.join (direct, glob_string)):
                    cls.update_database (session, file)
        session.commit ()
    
    @classmethod
    def update_database (cls, session, file_name):
        file = os.path.abspath (file_name)
        print ("Checking for file " + file + " in database")
        try:
            oldentry = session.query (DumpFileEntry).filter_by (file = file).one ()
            if oldentry.date >= datetime.datetime.fromtimestamp (os.path.getmtime(file)):
                print ("Newer or equivalent file in database: skipping")
                return
            else:
                session.delete (oldentry)
        except sqlalchemy.orm.exc.NoResultFound:
            print ("File not found in database: adding")
            pass
        d = records.dump.Dump (file)
        entry = DumpFileEntry (file = file, date = datetime.datetime.fromtimestamp (os.path.getmtime(file)), timestep = d.ncyc, state = d.getState ())
        cls.set_parameters (entry, d, records.dump.pfile)
        cls.set_parameters (entry, d, records.dump.qfile)
                        
        try:
            print ("Checking for a corresponding simulation in the database...")
            simulation = session.query (SimulationEntry).filter_by (runid = d.getRunId ()).one ()
            simulation.check_parameters (entry)
        except sqlalchemy.orm.exc.NoResultFound:
            print ("No appropriate simulation found: creating")
            simulation = SimulationEntry (runid = d.getRunId (), name = d.namep, path = os.path.dirname (file))
            simulation.copy_parameters (entry)
            session.add (simulation)
        simulation.dumpfiles.append (entry)
        session.add (entry)
        session.commit ()
        
        
    @classmethod
    def set_parameters (cls, entry, dump, file):
        parameter_file = open (file, 'r')
        for line in parameter_file:
            words = line.split ('\t')
            if words [0] [0] == '#':
                return
            try:
                setattr (entry, words [1], dump.parameters [words [1]].value)
            except KeyError as e:
                pass
        
    def get_data (self):
        if self.dataobject is None:
            self.dataobject = records.dump.DataDump (self.file)
        return self.dataobject
        
    def cache (self, session, cache_name, function):
        code = ''.join (inspect.getsourcelines (function) [0])
        for dumpcache in self.dumpcache:
            if dumpcache.name == cache_name and dumpcache.code == code:
                return dumpcache.get_value ()
            elif dumpcache.name == cache_name and dumpcache.code != code:
                session.delete (dumpcache)
                session.commit ()
        result = u.Quantity (function (self.get_data ()))
        self.dumpcache.append (DumpCache (name = cache_name, value = result.value, unit = str (result.unit), code = code))
        session.commit ()
        return result
        
class CNVFileEntry (Base):
    __tablename__ = 'cnvfiles'

    file = sqlalchemy.Column (sqlalchemy.String, primary_key=True)
    
    name = sqlalchemy.Column (sqlalchemy.String)
    runid = sqlalchemy.Column (sqlalchemy.String, sqlalchemy.ForeignKey ('simulations.id'))
    simulation = sqlalchemy.orm.relationship ("SimulationEntry", backref = sqlalchemy.orm.backref ('cnvfile', order_by = file, uselist = False))
    date = sqlalchemy.Column (sqlalchemy.DateTime)

    @sqlalchemy.orm.reconstructor
    def init_on_load (self):
        self.dataobject = None

    def __repr__(self):
       return "<CNV File (name='%s', timestep='%s', state=%s)>" % (self.simulation.name)
        
    def get_data (self):
        if self.dataobject == None:
            self.dataobject = records.cnv.CNVFile (self.file)
        return self.dataobject
        
    def cache (self, session, cache_name, function):
        code = ''.join (inspect.getsourcelines (function) [0])
        for cnvcache in self.cnvcache:
            if cnvcache.name == cache_name and cnvcache.code == code:
                return cnvcache.get_value ()
            elif cnvcache.name == cache_name and cnvcache.code != code:
                session.delete (cnvcache)
                session.commit ()
        result = u.Quantity (function (self.get_data ()))
        self.cnvcache.append (CNVCache (name = cache_name, value = result.value, unit = str (result.unit), code = code))
        session.commit ()
        return result
        
    @classmethod
    def scan_for_updates (cls, directory, glob_string = "*"):
        glob_string += ".cnv"
        session = Session ()
        for entry in session.query (CNVFileEntry).all ():
            if not os.path.isfile (entry.file):
                print ("File " + entry.file + " has been deleted: removing from database")
                session.delete (entry)
                session.commit ()
        for simentry in session.query (SimulationEntry).all ():
            if len (simentry.dumpfiles) == 0 and simentry.cnvfile is None:
                print ("Simulation " + simentry.name + " is empty: removing from database")
                session.delete (simentry)
                session.commit ()
        for file in glob.glob (os.path.join (directory, glob_string)):
            cls.update_database (session, file)
        for root, dirs, files in os.walk (directory):
            for direct in dirs:
                direct = os.path.join (root, direct)
                print ("Searching in " + direct + " for " + os.path.join (direct, glob_string))
                for file in glob.glob (os.path.join (direct, glob_string)):
                    cls.update_database (session, file)
        session.commit ()
        
    @classmethod
    def update_database (cls, session, file_name):
        name = os.path.basename (os.path.splitext (file_name) [0])
        file = os.path.abspath (file_name)
        print ("Checking for file " + file + " in database")
        try:
            oldentry = session.query (CNVFileEntry).filter_by (file = file).one ()
            if oldentry.date >= datetime.datetime.fromtimestamp (os.path.getmtime(file)):
                print ("Newer or equivalent file in database: skipping")
                return
            else:
                session.delete (oldentry)
        except sqlalchemy.orm.exc.NoResultFound:
            print ("File not found in database: adding")
            pass
        entry = CNVFileEntry (file = file, name = name, date = datetime.datetime.fromtimestamp (os.path.getmtime(file)))
        try:
            print ("Checking for a corresponding simulation in the database...")
            try:
                simulation = session.query (SimulationEntry).filter_by (path = os.path.dirname (file), name = name).one ()
            except sqlalchemy.orm.exc.MultipleResultsFound:
                print ("Multiple possible simulations, skipping file")
                return
        except sqlalchemy.orm.exc.NoResultFound:
            print ("No appropriate simulation found: creating")
            simulation = SimulationEntry (name = name, path = os.path.dirname (file))
            session.add (simulation)
        simulation.cnvfile = entry
        session.add (entry)
        session.commit ()
        
class DumpCache (Base):
    """docstring for DumpCache """
    __tablename__ = "dumpcache"
    
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key = True)
    name = sqlalchemy.Column (sqlalchemy.String)
    dumpfile = sqlalchemy.Column (sqlalchemy.String, sqlalchemy.ForeignKey (DumpFileEntry.file))
    dump = sqlalchemy.orm.relationship ("DumpFileEntry", backref = sqlalchemy.orm.backref ('dumpcache', order_by = name))
    code = sqlalchemy.Column (sqlalchemy.String)
    
    __table_args__ = (sqlalchemy.UniqueConstraint('name', 'dumpfile', name='uix_1')),
    
    value = sqlalchemy.Column (sqlalchemy.Float)
    unit = sqlalchemy.Column (sqlalchemy.String)
    
    def get_value (self):
        return self.value * u.Unit (self.unit)
        
class CNVCache (Base):
    """docstring for DumpCache """
    __tablename__ = "cnvcache"
    
    id = sqlalchemy.Column (sqlalchemy.Integer, primary_key = True)
    name = sqlalchemy.Column (sqlalchemy.String)
    cnvfile = sqlalchemy.Column (sqlalchemy.String, sqlalchemy.ForeignKey (CNVFileEntry.file))
    cnv = sqlalchemy.orm.relationship ("CNVFileEntry", backref = sqlalchemy.orm.backref ('cnvcache', order_by = name))
    code = sqlalchemy.Column (sqlalchemy.String)
    
    __table_args__ = (sqlalchemy.UniqueConstraint('name', 'cnvfile', name='uix_1')),
    
    value = sqlalchemy.Column (sqlalchemy.Float)
    unit = sqlalchemy.Column (sqlalchemy.String)
    
    def get_value (self):
        return self.value * u.Unit (self.unit)
        
        
# Base.prepare (engine)

def basicQuery (session):
    return session.query (SimulationEntry, DumpFileEntry).join (DumpFileEntry)

Base.metadata.create_all (engine)
