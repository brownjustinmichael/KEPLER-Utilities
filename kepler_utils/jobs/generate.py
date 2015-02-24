import os
import subprocess
import glob
import hashlib
import time
import datetime

from kepler_utils.database.database import SimulationEntry, DumpFileEntry, CNVFileEntry, Session

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
            session = Session ()
            query = session.query (SimulationEntry).filter (SimulationEntry.template_name == self.original).filter (SimulationEntry.template_hash == hashlib.md5 (open (self.original).read ().encode ()).hexdigest ()).filter (SimulationEntry.complete == True)
            for param in parameters:
                if type (parameters [param]) == float:
                    query = query.filter (getattr (SimulationEntry, param) > parameters [param] * 0.99)
                    query = query.filter (getattr (SimulationEntry, param) < parameters [param] * 1.01)
                else:
                    query = query.filter (getattr (SimulationEntry, param) == parameters [param])
            if query.count () != 0:
                print ("Simulation already exists")
                raise TypeError ("Simulation already exists")
                
            query = session.query (SimulationEntry).filter (SimulationEntry.template_name == self.original).filter (SimulationEntry.template_hash == hashlib.md5 (open (self.original).read ().encode ()).hexdigest ()).filter (SimulationEntry.complete == False)
            for param in parameters:
                if type (parameters [param]) == float:
                    query = query.filter (getattr (SimulationEntry, param) > parameters [param] * 0.99)
                    query = query.filter (getattr (SimulationEntry, param) < parameters [param] * 1.01)
                else:
                    query = query.filter (getattr (SimulationEntry, param) == parameters [param])
            for sim in query.all ():
                for entry in sim.dumpfiles:
                    session.delete (entry)
                for entry in sim.cnvfiles:
                    session.delete (entry)
                session.delete (sim)
                session.commit ()
                
            session.close ()
        
        if os.path.isfile (new_file):
            if not force:
                print ("Generator file already exists")
                raise NameError ("File already exists")
            
        file = open (new_file, 'w')
        file.write (self.content + '\n')
        
        for param in parameters:
            file.write ("p %s %s\n" % (param, str (parameters [param])))
            
        file.close ()
        
class Simulation (object):
    def __init__ (self, name, generator, run_location = '.', command = './kepler', force = False):
        if len (name) > 8:
            print ("KEPLER can only handle 8 character names")
            raise NameError ("KEPLER can only handle 8 character names")
        
        self.generator = generator
        self.force = force
        self.command = command
        self.run_location = run_location
        self.name = name
        
    def run (self, query = False, **kwargs):
        self.generator.generate (os.path.join (self.run_location, self.name + 'g'), force = True, query = query, **kwargs)
        if len (glob.glob (os.path.join (self.run_location, self.name + ".cnv"))) > 0 or len (glob.glob (os.path.join (self.run_location, self.name + "#*"))) > 0:
            if not self.force:
                print ("Dump files would clash with existing files.")
                raise NameError ("Dump files would clash with existing files")
            else:
                os.remove (os.path.join (self.run_location, self.name + ".cnv"))
                for file in glob.glob (os.path.join (self.run_location, self.name + "#*")):
                    os.remove (file)
        return subprocess.Popen ([self.command, self.name, self.name + 'g'], cwd = self.run_location)
        
    def rebase (self, tags = None):
        session = Session ()
        while True:
            try:
                DumpFileEntry.scan_for_updates (self.run_location, self.name + "#*", tags, template_name = self.generator.original, log_info = True)
                CNVFileEntry.scan_for_updates (self.run_location, self.name + ".cnv", tags, template_name = self.generator.original)
                return
            except Exception as e:
                # print ("There was an issue. Retrying in 5 seconds...")
                print (e)
                print (type (e))
                # time.sleep (5)
                raise (e)
                
    def getLatest (self):
        mostRecent = None
        timeStamp = 0.0
        for file in glob.glob (os.path.join (self.run_location, self.name + "#*")):
            name = file.split ('/') [-1].split ('#') [-1]
            newTime = datetime.datetime.fromtimestamp (os.path.getmtime(file))
            if mostRecent is None or timeStamp < newTime:
                mostRecent = name
                timeStamp = newTime
        return mostRecent
        
# Simulation ("s15o.1s1", Generator (osfactor = 0.1, scpower = 1.0), run_location = 'generator', command = '/Users/justinbrown/Codes/kepler/run/kepler')