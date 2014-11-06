import os
import subprocess
import glob
import hashlib
import time

import database.database as db

class Generator (object):
    """docstring for Generator """
    def __init__(self, original):
        super(Generator , self).__init__()
        
        infile = open (original, 'r')
        self.content = infile.read ()
        infile.close ()
        
        self.original = original
        
    def generate (self, new_file, force = False, query = False, **kwargs):
        parameters = kwargs
        # parameters = {'scpower' : scpower, 'osfactor' : osfactor}
        
        for param in parameters:
            if param [0] == 'p' and param [1:].isnumeric:
                parameters [param [1:]] = parameters.pop (param)
        
        if ('scpower' in parameters or 'scpower' in parameters) and 'woodscon' not in parameters:
            parameters ['woodscon'] = 1.0
            
        if ('osfactor' in parameters or 'osherwig' in parameters) and 'brumoson' not in parameters:
            parameters ['brumoson'] = 1.0
            
        if ('osherwig' in parameters):
            parameters ['osfactor'] = 0.0
        
        if query:
            session = db.Session ()
            query = session.query (db.SimulationEntry).filter (db.SimulationEntry.template_name == self.original).filter (db.SimulationEntry.template_hash == hashlib.md5 (open (self.original).read ().encode ()).hexdigest ()).filter (db.SimulationEntry.complete == True)
            for param in parameters:
                if type (parameters [param]) == float:
                    query = query.filter (getattr (db.SimulationEntry, param) > parameters [param] * 0.99)
                    query = query.filter (getattr (db.SimulationEntry, param) < parameters [param] * 1.01)
                else:
                    query = query.filter (getattr (db.SimulationEntry, param) == parameters [param])
            if query.count () != 0:
                raise TypeError ("Simulation already exists")
            session.close ()
                
        if os.path.isfile (new_file):
            if not force:
                raise NameError ("File already exists")
            
        file = open (new_file, 'w')
        file.write (self.content + '\n')
        
        for param in parameters:
            file.write ("p %s %s\n" % (param, str (parameters [param])))
            
        file.close ()
        
class Simulation (object):
    def __init__ (self, name, generator, run_location = '.', command = './kepler', force = False):
        if len (name) > 8:
            raise NameError ("KEPLER can only handle 8 character names")
            
        if len (glob.glob (os.path.join (run_location, name + ".cnv"))) > 0 or len (glob.glob (os.path.join (run_location, name + "#*"))) > 0:
            if not force:
                raise NameError ("Dump files would clash with existing files")
        
        self.generator = generator
        self.force = force
        self.command = command
        self.run_location = run_location
        self.name = name
        
    def run (self, query = False, **kwargs):
        self.generator.generate (os.path.join (self.run_location, self.name + 'g'), self.force, query = query, **kwargs)
        return subprocess.Popen ([self.command, self.name, self.name + 'g'], cwd = self.run_location)
        
    def rebase (self, tags = None):
        while True:
            try:
                db.DumpFileEntry.scan_for_updates (self.run_location, self.name + "#*", tags, template_name = self.generator.original)
                db.CNVFileEntry.scan_for_updates (self.run_location, self.name + ".cnv", tags, template_name = self.generator.original)
                return
            except:
                print ("There was an issue. Retrying in 5 seconds...")
                time.sleep (5)
        
# Simulation ("s15o.1s1", Generator (osfactor = 0.1, scpower = 1.0), run_location = 'generator', command = '/Users/justinbrown/Codes/kepler/run/kepler')