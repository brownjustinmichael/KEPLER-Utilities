import argparse
import atexit
import matplotlib.pyplot as plt

class PlotArgumentParser (argparse.ArgumentParser):
    """A subclass of argument parser for plotting routines to conveniently add command line options"""
    
    def __init__ (self, inputFile = True, *args, **kwargs):
        super(PlotArgumentParser, self).__init__(*args, **kwargs)
        if inputFile:
            self.add_argument ('input_file')
        self.add_argument ('--output', default = None)
        self.add_argument ('--style', default = None)
        
        self.namespace = None
        
    def parse_args (self, *args, **kwargs):
        self.namespace = super (PlotArgumentParser, self).parse_args ()
        
        if self.namespace.style is not None:
            plt.style.use (self.namespace.style)
            
        atexit.register (self.outputOrShow)
            
        return self.namespace
    
    def outputOrShow (self):
        try:
            plt.tight_layout ()
        except ValueError:
            print ("Issue with tight layout.")
        
        if (self.namespace.output is not None):
            plt.savefig (self.namespace.output)
        else:
            # Plot the result
            plt.show ()

