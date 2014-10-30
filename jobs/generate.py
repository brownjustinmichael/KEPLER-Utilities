import os
import subprocess
import glob

import database.database

class Generator (object):
    """docstring for Generator """
    def __init__(self, original):
        super(Generator , self).__init__()
        
        infile = open (original, 'r')
        self.content = infile.read ()
        infile.close ()
        
        self.original = original
        
    def generate (self, new_file, force = False, **kwargs):
        parameters = kwargs
        # parameters = {'scpower' : scpower, 'osfactor' : osfactor}
        
        for param in parameters:
            if param [0] == 'p' and param [1:].isnumeric:
                parameters [param [1:]] = parameters.pop (param)
        
        if ('scpower' in parameters or 'scpower' in parameters) and 'woodscon' not in parameters:
            parameters ['woodscon'] = 1.0
            
        if ('osfactor' in parameters or 'osherwig' in parameters) and 'brumoson' not in parameters:
            parameters ['brumoson'] = 1.0
        
        print ("Generating")
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
        print ("DONE")
        
    def run (self, **kwargs):
        self.generator.generate (os.path.join (self.run_location, self.name + 'g'), self.force, **kwargs)
        return subprocess.Popen ([self.command, self.name, self.name + 'g'], cwd = self.run_location)
        
    def rebase (self):
        database.database.DumpFileEntry.scan_for_updates (self.run_location, self.name)
        database.database.CNVFileEntry.scan_for_updates (self.run_location, self.name)
        
# Simulation ("s15o.1s1", Generator (osfactor = 0.1, scpower = 1.0), run_location = 'generator', command = '/Users/justinbrown/Codes/kepler/run/kepler')