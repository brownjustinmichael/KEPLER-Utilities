import os
import time
import re

import fortranfile
import numpy
import struct
import pandas
import astropy.units as u

__location__ = os.path.realpath (os.path.join (os.getcwd (), os.path.dirname (__file__)))

pfile = os.path.realpath (os.path.join (__location__, "parameters.dat"))
qfile = os.path.realpath (os.path.join (__location__, "qparameters.dat"))

class Isotope(object):
    """
    A python representation of an isotope
    """
    def __init__(self, string, a, z, **kwargs):
        super (Isotope, self).__init__()
        self.string = string
        self.a = a
        self.z = z
        self.data = {}
        self.label = self.makeLabel (self.string, self.a)
        for key in kwargs:
            self.data [key] = kwargs [key]
            
    def __str__ (self):
        return self.string
        
    def getLabel (self):
        """
        Return a LaTeX version of the string representation of the isotope
        """
        return self.label
        
    @staticmethod
    def makeLabel (string, a):
        """
        Return a LaTeX version of the string representation of an isotope
        """
        label = ''.join([i for i in string if not i.isdigit()])
        if string == 'nt1':
            label = '$n$'
        elif string == 'pn1':
            label = '$p$'
        elif string == 'ye':
            label = '$Y_{e}$'
        elif string == "'fe'":
            label = "$\\mathrm{\\lq{Fe}\\rq}$"
        else:
            label = '$^{' + str (int (a)) + '}\mathrm{' + label.capitalize () + '}$'
        return label

class Dump (object):
    """
    This base class contains most of the dump file information except the actual data
    """
    def __init__(self, fileName):
        super(Dump, self).__init__()
        self.parameters = {}
        self.file = fortranfile.FortranFile (fileName)
        
        self.state = fileName.split ('#') [-1]
        print ("Loading Dump from " + fileName)
        
        # Read dump file contents
        self._read ()
        
        self.file.close ()
        
    def getRunId (self):
        """
        Return a unique identifier for the run
        """
        return self.uuidrun
    
    def getCycle (self):
        """
        Return the timestep number of the dump file
        """
        return self.ncyc
    
    def getState (self):
        """
        Return the string 'state' of the dump file, e.g. 'cdep', 'presn', '2000'
        """
        return self.state
        
    def _readInt (self, kind = 8):
        """
        Read a single integer from self.file, kind is the fortran kind of integer to read (2, 4, 8 possible)
        """
        if kind == 16:
            return struct.unpack ('>' + str (1) + 'q', self.file.read (8)) [0]
        if kind == 8:
            return struct.unpack ('>' + str (1) + 'i', self.file.read (4)) [0]
        if kind == 4:
            return struct.unpack ('>' + str (1) + 'h', self.file.read (2)) [0]
        if kind == 2:
            return struct.unpack ('>' + str (1) + 'b', self.file.read (1)) [0]
        raise TypeError ("Unrecognized integer kind.")

    def _readInts (self, num, kind = 8):
        """
        Read num integers from self.file, kind is the fortran kind of integer to read (2, 4, 8 possible)
        Return a numpy array of integers
        """
        if kind == 16:
            return numpy.array (struct.unpack ('>' + str (num) + 'q', self.file.read (num * 8)))
        if kind == 8:
            return numpy.array (struct.unpack ('>' + str (num) + 'i', self.file.read (num * 4)))
        if kind == 4:
            return numpy.array (struct.unpack ('>' + str (num) + 'h', self.file.read (num * 2)))
        if kind == 2:
            return numpy.array (struct.unpack ('>' + str (num) + 'b', self.file.read (num)))
        else:
            raise TypeError ("Unrecognized integer kind.")

    def _readDouble (self, kind = 8):
        """
        Read a single double from self.file, kind is the fortran kind of real to read (8 possible)
        """
        if kind == 8:
            return struct.unpack ('>' + str (1) + 'd', self.file.read (8)) [0]
        if kind == 4:
            return struct.unpack ('>' + str (1) + 'f', self.file.read (4)) [0]

    def _readDoubles (self, num, kind = 8):
        """
        Read num doubles from self.file, kind is the fortran kind of real to read (8 possible)
        Return a numpy array of doubles
        """
        return numpy.array (struct.unpack ('>' + str (num) + 'd', self.file.read (num * kind)))

    def _readChars (self, num, asByteString = False):
        """
        Read num characters from self.file; if asByteString, return a byte string, else return a python string
        """
        chars = self.file.read (num)
        if asByteString:
            return chars
        chars = re.sub(r'[^ -~]', '', chars.decode ('utf-8'))
        return chars.rstrip ()
        
    def _readStrings (self, num, charsPer, asByteString = False):
        """
        Read num strings from self.file, if asByteString, return a list of byte strings, else return a list of python strings
        """
        return [self._readChars (charsPer, asByteString) for i in range (num)]
        
    def _skip (self, num):
        """
        Skip num bytes in self.file
        """
        self.file.read (num)
        
    def _read (self):
        """
        Read in only the critical information about the dump file:
            -Header information
            -Parameters
            -Qparameters
            -UUID Information
        """
        
        # Read off the first few bytes
        test = self._readInt ()
        
        # Check the version; this is only designed to read dump files from versions as early as 170011
        self.version = self._readInt ()
        print ("Version is " + str (self.version))
        # assert (self.version >= 170010)
        
        # Read the header information
        self._readHeader ()
        
        # Read the parameters
        self._readParameters (self.parameters, self.numParameters, os.path.join (__location__, 'parameters.dat'))
        self._readParameters (self.parameters, self.numDerivedParameters, os.path.join (__location__, 'qparameters.dat'), unitCol = 3, defaultCol = None)
        
        if self.version < 170011:
            self.parameters ['iold'] += 8
        
        #Skip the small arrays
        self._skip ((self.lenSmallArrays + 1) * 8)

        # Read the program name
        self.namep = self._readChars (8)
        
        # Skip the small character arrays, zones and ppn
        self._skip ((self.lenSmallCharArrays - 1 + (self.jmsave + 1 + 1) * min (self.numCenteredZoneArrays + self.numInterfacialZoneArrays, 33) + self.numNetworkIons * (self.jmsave - 1) + 1) * 8)
        
        # Magnetic fields should be included here
        assert (self.parameters ['magnet'] <= 0.0)
        
        # Skip viscosity and diffusion arrays
        self._skip ((self.jmsave + 2) * 2 * 8)
        
        # TODO Include flame data reader
        assert (self.parameters ['sharp1'] <= 0.0)
        
        # TODO Include WIMP data reader
        assert (self.parameters ['wimp'] <= 0.0)
        
        # Skip viscosity and diffusion arrays
        self._skip ((self.jmsave - 1) * 8)
        
        # Read in UUIDs
        if self.version >= 170004:
            self.uuidrun, self.uuidcycle, self.uuiddump, self.uuidprev, self.uuidprog, self.uuidexec = self._readStrings (6, 16, True)
            self.nuuidhist = self._readInt (16)
            self.uuidhist = [self._readStrings (self.nuuidhist, 16, True) for i in range (6)]
        else:
            self.uuidrun = None
        
    def _readHeader (self):
        """
        Read in the header of the dump file and set the appropriate variables
        """
        self.ncyc = self._readInt ()
        self.lenHeader, self.maxZones, self.maxBurnZones, self.nburn, self.iratioz, self.nvar = self._readInts (6)
        self.nheadz, self.nhead, self.numParameters, self.nparm, self.numDerivedParameters, self.nqparm = self._readInts (6)
        self.numNetworks, self.maxNetworkIons, self.numTotalIons, self.numNetworksb, self.maxNetworkIonsb, self.numTotalIonsb = self._readInts (6)
        self.nreacz, self.maxTimestepControllers, self.ndt, self.maxPistons, self.maxYeInitializations, self.nsubz, self.nsub = self._readInts (7)
        self.numInterfacialZoneArrays, self.numCenteredZoneArrays, self.nzoneb, self.lenSmallArrays, self.lenSmallCharArrays = self._readInts (5)
        self.numNetworkIons, self.numBurnIons, self.nreac, self.jmsave, self.lencom, self.lencomc = self._readInts (6)
        self.nedtcom, self.ndatqz, self.ngridz, self.nylibz, self.nyoffst = self._readInts (5)
        self.lenshed, self.lenqhed, self.nzedz, self.ncsavdz = self._readInts (4)
        things = self._readInts (self.iratioz * self.lenHeader - 47)
        self._readDouble ()
        self._readDouble ()
            
    def _readParameters (self, dic, numParameters, tableFile, comment = '#', unitCol = 4, defaultCol = 3):
        """
        Read in the parameters from self.file using tableFile as a guide; comment is the comment character for tableFile, unitCol is the column of the unit string in tableFile
        Puts parameters in dic as astropy Quantity objects, complete with units
        """
        # Open the parameter file
        file = open (tableFile, 'r')
        i = 0
        for line in file:
            words = line.split ('\t')
            if len (words) == 0 or words [0] == comment:
                # Check that this isn't a comment line
                continue
            elif len (words) > 6 and words [6] != '\n' and words [6] != '' and int (words [6].rstrip ('\n')) > self.version:
                # Check that this parameter wasn't added after the current version
                print ("Skipping parameter " + words [1] + " because this version of KEPLER predates it")
                continue
            else:
                i += 1
                # Check that we haven't exceeded the number of parameters in the dump file
                if i > numParameters:
                    # If we have a default, use the default
                    print ("WARNING: DEFAULTING for " + words [1])
                    if defaultCol != None:
                        if words [2] == 'float':
                            dic [words [1]] = float (words [3])
                        elif words [2] == 'integer':
                            dic [words [1]] = int (words [3])
                    else:
                        continue
                else:
                    #Check the type of the parameter and read it in
                    if words [2] == 'float':
                        dic [words [1]] = self._readDouble ()
                    
                    elif words [2] == 'integer':
                        dic [words [1]] = self._readInt ()
                    
                        #If it's an integer, read off the remainder of the 8 byte chunk
                        self._readInt ()
                    
                    else:
                        # Unrecognized type, raise an error
                        raise TypeError (words [2])
            # Read in units
            if words [unitCol] != '-':
                dic [words [1]] = u.Quantity (dic [words [1]], words [unitCol], dtype = type (dic [words [1]]))
            else:
                dic [words [1]] = u.Quantity (dic [words [1]], "1", dtype = type (dic [words [1]]))
                
        # Check that we've read all the parameters in the dump file; if there are extras, give them names of 'parm###', type double, and no units
        for j in range (i + 1, numParameters + 1):
            dic ['parm' + str (j)] = u.Quantity (self._readDouble (), '1', dtype = float)
            
        # Read off the remainder of the chunk
        self._readDoubles (2)

class DataDump (Dump):
    """
    This class takes KEPLER dump file as its argument and produces an indexable object containing the star data
    """
    def __init__(self, fileName, loadBurn = False):
        super(DataDump, self).__init__ (fileName)
        self.data = {}
        self.units = {}
        self.stardata = {}
        self.file = fortranfile.FortranFile (fileName)
        
        # Read dump file contents
        self._read_full (loadBurn)
        
        self.file.close ()
    
    def __getitem__ (self, index):
        """
        Get the index from the pandas dataframe object, adding in units when possible
        Returns an astropy Quantity object
        """
        # If possible, return quantity object
        try:
            return u.Quantity (self.df [index])
        except Exception as e:
            print (type (e))
            return self.df [index]
    
    def getIsotope (self, string):
        """
        Given a string, find the equivalent isotope object in the network
        """
        if string not in [str (i) for i in self.ions]:
            raise IndexError (string)
        for key in self.df:
            if str (key) == string:
                return key
        raise IndexError (string)

    def getIsotopes (self):
        """
        Return an iterator to iterate over the isotopes
        """
        return iter (self.ions)
    
    def getCore (self, string):
        """
        Returns an astropy quantity of the core mass associated with the isotope string
        """
        isos = [str (iso) for iso in self.getIsotopes ()]
        index = len (self [string]) - 1 - numpy.argmax (numpy.diff (self [string] < self.df [isos].max (1)) [::-1])
        if index == len (self [string]) - 1:
            return u.Quantity (0.0, 'g')
        return self.df ['mass coordinate'] [index]
    
    def _read_full (self, loadBurn = False):
        """
        Read the dump file into the DataDump object; if loadBurn, read the full burn array, else read only the approx network isotopes
        """
        self._readInt ()
        
        self.version = self._readInt ()

        # Skip the header and parameters; these are read in during the superclass constructor
        self._skip ((self.iratioz * self.lenHeader + 3) * 4)
        self._skip ((self.numParameters + 2 + self.numDerivedParameters + 2) * 8)
        
        # Read in the 'small' arrays
        self._readSmallArrays (self.data, self.units)
        self._readSmallCharArrays (self.data)
        
        # Read in the interfactial and centered zone arrays
        self._readZones (self.stardata, self.units, self.jmsave, self.numInterfacialZoneArrays, self.numCenteredZoneArrays, 33)
        
        # Read in the abundance data for the small network
        ppn = numpy.array ([self._readDoubles (self.numNetworkIons) for i in range (self.jmsave - 1)])
        self._readDouble ()
        
        # Magnetic fields should be included here
        assert (self.parameters ['magnet'] <= 0.0)
        
        # Read in the viscosity and diffusion arrays
        self.stardata ['angdgeff'] = numpy.zeros (self.jmsave + 1)
        self.stardata ['angdgeff'] [:-1] = u.Quantity (self._readDoubles (self.jmsave + 1) [:-1], u.cm ** 2 / u.s)
        self.stardata ['difieff'] = numpy.zeros (self.jmsave + 1)
        self.stardata ['difieff'] [:-1] = u.Quantity (self._readDoubles (self.jmsave + 1) [:-1], u.cm ** 2 / u.s)
        
        # Flame data should be included here
        assert (self.parameters ['sharp1'] <= 0.0)
        
        # WIMP data should be included here
        assert (self.parameters ['wimp'] <= 0.0)
        
        # Read in the convective energy contribution
        self.stardata ['sadv'] = numpy.zeros (self.jmsave + 1)
        self.stardata ['sadv'] [:-1] = self._readDoubles (self.jmsave + 1) [:-1]
        
        # Skip UUID information
        if self.version >= 170004:
            self._skip (6 * 16 + 8 + 16 * 6 * self.nuuidhist)
            self._readDoubles (1)
        
        # TODO Implement log data reader
        if (self.version >= 168450):
            self.nlog = self._readInt ()
            self._readInt ()
            assert (self.nlog == 0)
            self._readDouble ()
        
        # Read energy dissipation due to shear
        self.stardata ['sv'] = numpy.zeros (self.jmsave + 1)
        self.stardata ['sv'] [1:-2] = self._readDoubles (self.jmsave - 2) 
        
        # Read out user-defined parameters
        # TODO Implement user-defined parameter reader
        self.data ['noparm'] = self._readInt ()
        assert (self.data ['noparm'] == 0)
        self._readDoubles (1)
        if (self.parameters ['imaxb'] <= 0 or loadBurn == False):
            loadBurn = False
        else:
            self._readInt ()
            self._readDoubles (1)
        
            # Read burn zone data
            self._readBurnZones (self.stardata, self.units, self.jmsave, self.parameters ['nsaveb'])
        
            # Read in the abundance data for the large network
            ppnb = numpy.array ([self._readDoubles (self.numBurnIons) for i in range (self.jmsave - 1)])
            self._readDoubles (1)
        
            # Read in burn ions informational arrays
            self.maxBurnIons = self._readInt ()
            self._readInts (2)
            self.data ['nabmax'] = self._readInts (self.maxBurnIons + 2) [:-2]
            self.data ['nzbmax'] = self._readInts (self.maxBurnIons + 2) [:-2]
            self.data ['nibmax'] = self._readInts (self.numBurnIons + 2) [:-2]
            self.data ['ionbmax'] = self._readStrings (self.maxBurnIons, 8)
            self._readChars (8, True)
            self.data ['burnamax'] = self._readDoubles (self.maxBurnIons + 1) [:-1]
            self.data ['burnmmax'] = self._readDoubles (self.maxBurnIons + 1) [:-1]
            self.data ['ibcmax'] = self._readInts (self.maxBurnIons + 2) [:-3]
        
            # Read in surface composition
            self.data ['compsurfb'] = self._readDoubles (self.numBurnIons)
        
        # Load ppn data into abundance arrays
        if not loadBurn:
            self._loadAbundances (self.stardata, ppn, self.data ['zionn'], self.numNetworkIons, self.stardata ['netnum'], self.jmsave, self.data ['ions'], self.data ['aion'], self.data ['zion'], self.numTotalIons, self.numNetworks)
        else:
            self._loadAbundances (self.stardata, ppnb, self.data ['zionnb'], self.numBurnIons, self.stardata ['netnumb'], self.jmsave, self.data ['ionbmax'], self.data ['nabmax'], self.data ['nzbmax'], self.numBurnIons, self.numNetworksb)
        
        # Convert stardata into a pandas dataframe object
        starItems = [(key, list (self.stardata [key] [1:-1])) for key in self.stardata]
        starItems.append (("mass coordinate", list (self.parameters ['totm'] - self.stardata ['ym'] [1:-1])))
        self.df = pandas.DataFrame.from_items (starItems)
        # self.df = pandas.DataFrame (self.stardata)
        
    def _readSmallArrays (self, dic, unitDic):
        """
        Read in the small arrays:
            -Piston Arrays
            -Ion informational arrays
            -Burn ion informational arrays
            -Timestep controllers
            -Rate arrays
            -Surface composition
            -Wind arrays
        """
        # Read piston arrays
        dic ['tpist'] = self._readDoubles (self.maxPistons) * u.s
        dic ['rpist'] = self._readDoubles (self.maxPistons) * u.cm
        dic ['yemass'] = self._readDoubles (self.maxYeInitializations) * u.g
        dic ['yeq0'] = self._readDoubles (self.maxYeInitializations) * u.mol / u.g

        # Read ion informational arrays
        dic ['aion'] = self._readDoubles (self.numTotalIons) [:-1]
        dic ['zion'] = self._readDoubles (self.numTotalIons) [:-1]
        dic ['znumi'] = self._readInts (self.numNetworks * 2) [:self.numNetworks]
        dic ['zionn'] = numpy.array ([self._readInts (self.maxNetworkIons) [:-1] for i in range (2 * self.numNetworks)]) [:self.numNetworks]
        
        # Read burn ion informational arrays
        dic ['aionb'] = self._readDoubles (self.numBurnIons) [:-1]
        dic ['zionb'] = self._readDoubles (self.numBurnIons) [:-1]
        dic ['znumib'] = self._readInts (self.numNetworksb * 2) [:self.numNetworksb]
        dic ['zionnb'] = numpy.array ([self._readInts (self.numBurnIons) for i in range (2 * self.numNetworksb)]) [:self.numNetworksb]
        
        # Read timestep controllers
        dic ['dtc'] = self._readDoubles (self.maxTimestepControllers) * u.s
        dic ['zjdtc'] = self._readDoubles (self.maxTimestepControllers)
        
        # Read subroutine timing
        dic ['timeused'] = numpy.array ([self._readDoubles (3) for i in range (self.nsubz + 1)]) * u.s

        # Read nuclear rate arrays
        dic ['totalr'] = self._readDoubles (self.nreacz) * u.mol
        dic ['ratr'] = self._readDoubles (self.nreacz) * u.mol / u.s
        dic ['qval'] = self._readDoubles (self.nreacz) * u.MeV
        dic ['zjrate'] = self._readDoubles (self.nreacz) * u.Unit ('1')
        dic ['rrx'] = self._readDoubles (self.nreacz) * u.mol / u.s / u.g
        
        # Read surface composition
        dic ['compsurf'] = self._readDoubles (self.maxNetworkIons)
        
        # Read q data
        # TODO What is this?
        dic ['zlocqz'] = self._readDoubles (self.ndatqz)
        dic ['zlocqz0'] = self._readDoubles (self.ndatqz)
        dic ['ratzdump'] = self._readDoubles (self.ndatqz)
        dic ['rationdez'] = self._readDoubles (self.ndatqz)
        dic ['ratioadz'] = self._readDoubles (self.ndatqz)
        
        # Read in user-specified zonal edit arrays
        dic ['zndatzed'] = self._readDoubles (self.nzedz)
        dic ['zncyczed'] = self._readDoubles (self.nzedz)
        dic ['zedmass1'] = self._readDoubles (self.nzedz)
        unitDic ['zedmass1'] = self.parameters ['scalem']
        dic ['zedmass2'] = self._readDoubles (self.nzedz)
        unitDic ['zedmass2'] = self.parameters ['scalem']
        
        # Read in wind arrays
        dic ['wind'] = self._readDoubles (self.numTotalIons + 1) [:-1]
        dic ['windb'] = self._readDoubles (self.numBurnIons + 1) [:-1]
        
    def _readSmallCharArrays (self, dic):
        # Read in names, flags, id words
        dic ['namep'] = self._readChars (8) # Name of program
        dic ['namec'] = self._readChars (16) # Name of command
        dic ['iflag8'] = self._readChars (8) # Sequence letter for naming graphics files
        dic ['iqbrnflg'] = self._readChars (8) # Communications flag for subroutine SDOTQ
        dic ['craybox'] = self._readChars (8) # Name of users computer output box
        dic ['idword'] = self._readChars (8) # ID flag, usually user name
        
        # Read in storage/status info
        dic ['nxdirect'] = self._readChars (16) # Name of storage directory
        dic ['lastrun'] = self._readChars (16) # Date of last run
        dic ['lastmod'] = self._readChars (16) # Supposedly date of last code compilation, but actually garbage
        
        # Read in ion labels
        dic ['ions'] = self._readStrings (self.numTotalIons, 8)
        dic ['ionsb'] = self._readStrings (self.numBurnIons, 8)
        
        # Read in timestep controller labels
        dic ['idtcsym'] = self._readStrings (self.maxTimestepControllers, 8)
        
        # Read in nuclear reaction symbols
        dic ['isymr'] = self._readStrings (self.nreacz, 8)
        
        # Read in post-processor file names and arrays
        dic ['nameqq'], dic ['nameqlib'], dic ['nameolds'], dic ['namenews'] = self._readStrings (4, 16)
        dic ['namedatq'] = self._readStrings (self.ndatqz, 8)
        dic ['labldatq'] = self._readStrings (self.ndatqz, 48)
        
        # Read in edit variable arrays and saved command strings
        dic ['namedzed'] = [self._readStrings (self.nzedz, 8) for i in range (10)]
        dic ['savdcmd0'] = self._readStrings (30, 80)
        
        # Read in 'look' variable arrays and post-processor file names
        dic ['namedatl'] = self._readStrings (self.ndatqz, 8)
        dic ['labldatl'] = self._readStrings (self.ndatqz, 48)
        dic ['nameqql0'], dic ['nameqql1'], dic ['nameqlbl'] = self._readStrings (3, 16)
        
        # Miscellaneous other symbols
        dic ['nsdirect'] = self._readChars (48)
        dic ['isosym'] = self._readStrings (50, 8)
        dic ['isoicon'] = self._readStrings (50, 16)
        dic ['savedcmd'] = self._readStrings (self.ncsavdz, 80)
        dic ['datapath'] = self._readChars (80)
        
        # Read in convection icons
        self.stardata ['icon'] = self._readStrings (self.jmsave + 1, 8)
        self._readDouble ()
        
    def _readZones (self, dic, unitsDic, jm, numInterfacialZoneArrays, numCenteredZoneArrays, nzonemax):
        """
        Read the interfacial zones then the centered zones
        """
        # Read the interfacial zones using the following zone names and units
        zoneiNames = ['ym', 'rn', 'rd', 'un', 'xln', 'qln', 'qld', 'difi', 'vconvect', 'oslen', 'adindex']
        zoneiUnits = [u.g, u.cm, u.cm, u.cm / u.s, u.erg / u.s, u.erg / u.s, u.erg / u.s, u.cm ** 2 / u.s]
    
        for i in range (min (numInterfacialZoneArrays, 33)):
            # If there is no zone name for this zone, use 'zonei###'
            if i < len (zoneiNames):
                dic [zoneiNames [i]] = self._readDoubles (jm + 2) [0:-1]
            else:
                dic ['zonei' + str (i)] = self._readDoubles (jm + 2) [0:-1]
            
            # Read in units
            if i < len (zoneiUnits):
                # unitsDic [zoneiNames [i]] = zoneiUnits [i]
                dic [zoneiNames [i]] = u.Quantity (dic [zoneiNames [i]], zoneiUnits [i], dtype = dic [zoneiNames [i]].dtype)
            
        # Read the centered zones using the following zone names, units, and types
        zonecNames = ['netnum', 'xm', 'dn', 'tn', 'td', 'en', 'pn', 'zn', 'etan', 'sn', 'snn', 'abar', 'zbar', 'xkn', 'xnei', 'stot', 'angj','angdg', 'angddsi', 'angdshi', 'angdssi', 'angdez', 'angdgsf', 'dsold', 'tsold', 'snold', 'snbd', 'snbt', 'abarold', 'abarnbd', 'abarnbt', 'ypbtime', 'ynbtime']
        zonecUnits = [1, u.g, u.g / u.cm ** 3, u.K, u.K, u.erg / u.g, u.erg / u.cm ** 3, u.erg, 1, u.erg / u.g / u.s, u.erg / u.g / u.s, u.g / u.mol, 1, u.cm ** 2 / u.g, 1.0 / u.cm ** 3, u.k, u.g / u.cm ** 3, u.K, u.erg / u.g / u.s, u.erg * u.cm ** 3 / u.s / u.g ** 2, u.erg / u.g / u.s / u.K, 1, u.cm ** 3 / u.g, 1.0 / u.K, u.mol / u.g / u.s, u.mol / u.g / u.s]
        zonecIsint = [True] + [False] * (len (zonecNames) - 1)
    
        for i in range (min (numCenteredZoneArrays, nzonemax - numInterfacialZoneArrays)):
            #If there is no zone name or type, use 'zonec###', assume double and unitless
            if i < len (zonecNames):
                name = zonecNames [i]
                isint = zonecIsint [i]
            else:
                name = 'zonec' + str (i)
                isint = False
            
            if isint:
                dic [name] = self._readInts (2 * (jm + 2)) [:jm + 1]
            else:
                dic [name] = self._readDoubles (jm + 2) [0:-1]
            
            # Read in units
            if i < len (zonecUnits):
                # unitsDic [name] = zonecUnits [i]
                dic [name] = u.Quantity (dic [name], zonecUnits [i], dtype = dic [name].dtype)
            
    
    def _readBurnZones (self, dic, unitDic, jm, nsaveb):
        """
        Read in the burn zone arrays
        """
        # Read using the following names, units, and types
        zonebNames = ['netnumb', 'zlimab', 'zlimzb', 'zlimcb', 'timen', 'dtimen', 'dnold', 'tnold', 'ymb', 'sburn', 'etab', 'pbuf']
        zonebUnits = [1, 1, 1, 1, u.s, u.s, u.g / u.cm ** 3, u.K, u.g, u.erg / u.g / u.s, u.mol / u.g, u.mol / u.g]
        zonebIsint = [True] * 4 + [False] * (len (zonebNames) - 4)
    
        for i in range (nsaveb):
            # If name isn't in name list, use 'zoneb###', assume double and unitless
            if i < len (zonebNames):
                name = zonebNames [i]
                isint = zonebIsint [i]
            else:
                name = 'zoneb' + str (i)
                isint = False
            
            if isint:
                dic [name] = self._readInts (2 * (jm + 2)) [:self.jmsave + 1]
            else:
                dic [name] = self._readDoubles (jm + 2) [:-1]
            
            # Read in units
            if i < len (zonebUnits):
                # unitDic [name] = zonebUnits [i]
                dic [name] = u.Quantity (dic [name], zonebUnits [i], dtype = dic [name].dtype)
        
    def _loadAbundances (self, dic, ppn, ionn, imax, netnum, jmsave, ions, aion, zion, numTotalIons, nn):
        """
        Load abundance arrays from PPN data
        """
        self.ions = []
        test = numpy.zeros ((numTotalIons, jmsave + 1))
        
        for j in range (jmsave - 1):
            net  = netnum [j + 1] - 1
            test [ionn [net] - 1,j + 1] = ppn [j] * aion [ionn [net] - 1]
            
        for i in range (numTotalIons):
            if ions [i] != '':
                self.ions.append (Isotope (ions [i], aion [i], zion [i]))
                dic [ions [i]] = test [i]

# For testing
# d = DataDump ("s15o0s0#presn")
