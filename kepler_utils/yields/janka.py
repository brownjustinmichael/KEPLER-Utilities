import collections
import astropy.units as u

class JankaParser (object):
    keys = {"Neutrino-driven explosion result without calculation of fallback:\n" : "without_fallback", "With fallback:\n" : "with_fallback", "Trajectories:\n" : "trajectories", "90% trajectory:\n" : "90%", "95% trajectory:\n" : "95%", "special trajectory:\n" : "special"}
    units = {"M_sun" : u.solMass, "foe" : 10.**51 * u.erg}
    
    def __init__ (self, calibration, records):
        self.calibration = calibration
        self.records = collections.OrderedDict (sorted (records.items (), key = lambda x: x [1] ['M']))
            
    @classmethod
    def readFrom (cls, filename):
        fileobj = open (filename, "r")
        calibration = fileobj.readline ().split () [1]
        fileobj.readline ()
        fileobj.readline ()
        
        records = collections.OrderedDict ()
        
        try:
            while (True):
                tup = cls.readRecord (fileobj)
                records [tup [0]] = tup [1]
        except ValueError:
            pass
            
        return cls (calibration = calibration, records = records)
        
    @staticmethod
    def readRecord (fileobj):
        record = {}
        line = fileobj.readline ()
        if line == "":
            raise ValueError

        name = line.split () [0]
        record ["M"] = float (name [1:]) * u.solMass
        record ["metallicity"] = name [0]
        fileobj.readline ()
        curdict = ()
        
        while (True):
            pos = fileobj.tell ()

            line = fileobj.readline ()
            posmid = fileobj.tell ()
            nextline = fileobj.readline ()

            if nextline == "------------\n":
                fileobj.seek (pos)
                return (name, record)
            fileobj.seek (posmid)
            if line == "\n":
                if (curdict != ()):
                    curdict = curdict [:-1]
                continue
            thisdict = record
            for index in curdict:
                thisdict = thisdict [index]
            words = line.split ()
            if words == []:
                return (name, record)
            if line in JankaParser.keys:
                thisdict [JankaParser.keys [line]] = {}
                curdict = curdict + (JankaParser.keys [line],)
                continue
            if line [-2] == ":":
                thisdict [line [:-2]] = {}
                curdict = curdict + (line [:-2])
                continue
            if words [1] == "=":
                thisdict [words [0]] = float (words [2]) * JankaParser.units [words [3]]
                continue
            raise TypeError ()
            
    def __iter__ (self):
        return iter (self.records)
        
    def __in__ (self, record):
        return record in self.records
        
    def __getitem__ (self, index):
        indices = index.split (":")
        results = [self.records]
        for index in indices:
            if index == "":
                newresults = []
                for result in results:
                    for inindex in result:
                        newresults.append (result [inindex])
                results = newresults
            else:
                results = [result [index] if result != 0.0 and index in result else {} for result in results]
        resultunit = None
        for result in results:
            if not isinstance (result, dict):
                try:
                    resultunit = result.unit
                    break
                except Exception as e:
                    print (e)
        if resultunit is not None:
            for i in range (len (results)):
                if isinstance (results [i], dict):
                    results [i] = u.Quantity (0.0 * resultunit)
            results = u.Quantity (results)
        return results
        
    def __add__ (self, other):
        for selfRecord in self:
            if selfRecord in other:
                raise TypeError ("Repeated mass")
        if self.calibration != other.calibration:
            raise TypeError ("Different calibrators used")
        records = collections.OrderedDict ()
        for selfRecord in self:
            records [selfRecord] = self.records [selfRecord]
        for otherRecord in other:
            records [otherRecord] = other.records [otherRecord]
        return JankaParser (calibration = self.calibration, records = records)