import argparse
import atexit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

puor = cm.get_cmap ("PuOr")
brbg = cm.get_cmap ("BrBG")

colors = []

for i in range (puor.N):
    if i <128:
        colors.append (puor (i))
    else:
        colors.append (brbg (i))

BGOr = ListedColormap (colors, 'BGOr')
plt.register_cmap (cmap = BGOr)

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
            
        # atexit.register (self.outputOrShow)
            
        return self.namespace
    
    def exit (self, fig = None):
        if fig is not None:
            fig.patch.set_alpha (0.0)

        try:
            plt.tight_layout ()
        except ValueError:
            print ("Issue with tight layout.")
        
        if (self.namespace.output is not None):
            if self.namespace.output.split (".") [-1] == "png":
                plt.savefig (self.namespace.output, dpi = 1000)
            else:
                plt.savefig (self.namespace.output)
        else:
            # Plot the result
            plt.show ()

