import fortranfile
import numpy
import struct
import astropy.units as u

def selected_kind (kind_length):
    if kind_length == 4:
        return 2
    if kind_length == 2:
        return 1
    else:
        raise TypeError ("Unrecognized kind length %i" % kind_length)
        
class CNVFile:
    """
    This class reads a KEPLER cnv file into memory.
    
    It requires a file_name parameter, and can be given a maximum model number to read as well.
    
    This object stores the cnv file as a list of dictionaries, each holding all the cnv file keywords for a given record. Indexing the object indexes that list.
    """
    def __init__ (self, file_name, max_model = None):
        self.models = []
        self.units = {}
        self.units ['dt'] = u.s
        self.units ['timesec'] = u.s
        self.units ['xmcoord'] = u.g
        self.units ['rncoord'] = u.cm
        self.units ['xlum_cnv'] = u.erg / u.s
        self.file = fortranfile.FortranFile (file_name)
        i = 0
        tbase = 0.0
        while (True):
#             print ("Attempting iteration %i" % i)
            try:
                self.read_int ()
                version = self.read_int ()
            except struct.error:
                print ("End of file.")
                return
            
#             print ("Found version, filling variables...")

            if version != 10500:
                raise TypeError ("Can't handle version type %i" % self.models ['version'])
            ncyc = self.read_int ()
            timesec = self.read_double ()
            dt = self.read_double ()
                        
            if (i > 1):
                if timesec + tbase < self [-1] ['timesec']:
                    tbase = self [-1] ['timesec']
                    print ("Time reset at %d." % i)
                    # return
                
            self.models.append ({})
            timesec += tbase 
            self [i] ['ncyc'] = ncyc
            self [i] ['timesec'] = timesec
            self [i] ['dt'] = dt

            self [i] ['nconv'] = self.read_int ()
            self [i] ['nnuc'] = self.read_int ()
            self [i] ['nnuk'] = self.read_int ()
            self [i] ['nneu'] = self.read_int ()
            self [i] ['nnucd'] = self.read_int ()
            self [i] ['nnukd'] = self.read_int ()
            self [i] ['nneud'] = self.read_int ()
            self [i] ['ncoord'] = self.read_int ()
            self [i] ['idx_kind_len'] = self.read_int ()
            self [i] ['nuc_kind_len'] = self.read_int ()

            self [i] ['nuc'] = self.read_ints (self [-1] ['nnuc'], self [-1] ['nuc_kind_len'])
            self [i] ['nuk'] = self.read_ints (self [-1] ['nnuk'], self [-1] ['nuc_kind_len'])
            self [i] ['neu'] = self.read_ints (self [-1] ['nneu'], self [-1] ['nuc_kind_len'])
            self [i] ['nucd'] = self.read_ints (self [-1] ['nnucd'], self [-1] ['nuc_kind_len'])
            self [i] ['nukd'] = self.read_ints (self [-1] ['nnukd'], self [-1] ['nuc_kind_len'])
            self [i] ['neud'] = self.read_ints (self [-1] ['nneud'], self [-1] ['nuc_kind_len'])
            self [i] ['yzip'] = self.read_chars (self [-1] ['nconv'])
            self [i] ['xmcoord'] = self.read_doubles (self [-1] ['ncoord'])
            self [i] ['rncoord'] = self.read_doubles (self [-1] ['ncoord'])

            self [i] ['ladv'] = self.read_int ()
            self [i] ['iadv'] = self.read_ints (self [-1] ['ladv'], self [-1] ['idx_kind_len'])
            self [i] ['dmadv'] = self.read_doubles (self [-1] ['ladv'])
            self [i] ['dvadv'] = self.read_doubles (self [-1] ['ladv'])

            self [i] ['inuc'] = self.read_ints (self [-1] ['nnuc'], self [-1] ['idx_kind_len'])
            self [i] ['inuk'] = self.read_ints (self [-1] ['nnuk'], self [-1] ['idx_kind_len'])
            self [i] ['ineu'] = self.read_ints (self [-1] ['nneu'], self [-1] ['idx_kind_len'])
            self [i] ['inucd'] = self.read_ints (self [-1] ['nnucd'], self [-1] ['idx_kind_len'])
            self [i] ['inukd'] = self.read_ints (self [-1] ['nnukd'], self [-1] ['idx_kind_len'])
            self [i] ['ineud'] = self.read_ints (self [-1] ['nneud'], self [-1] ['idx_kind_len'])
            self [i] ['iconv'] = self.read_ints (self [-1] ['nconv'], self [-1] ['idx_kind_len'])

            self [i] ['levcnv'] = self.read_int ()
            self [i] ['minloss'] = self.read_int ()
            self [i] ['mingain'] = self.read_int ()
            self [i] ['minnucl'] = self.read_int ()
            self [i] ['minnucg'] = self.read_int ()
            self [i] ['minneul'] = self.read_int ()
            self [i] ['minneug'] = self.read_int ()

            self [i] ['minlossd'] = self.read_int ()
            self [i] ['mingaind'] = self.read_int ()
            self [i] ['minnucld'] = self.read_int ()
            self [i] ['minnucgd'] = self.read_int ()
            self [i] ['minneuld'] = self.read_int ()
            self [i] ['minneugd'] = self.read_int ()

            self [i] ['tc_cnv'] = self.read_double ()
            self [i] ['dc_cnv'] = self.read_double ()
            self [i] ['pc_cnv'] = self.read_double ()
            self [i] ['ec_cnv'] = self.read_double ()
            self [i] ['sc_cnv'] = self.read_double ()
            self [i] ['ye_cnv'] = self.read_double ()
            self [i] ['ab_cnv'] = self.read_double ()

            self [i] ['et_cnv'] = self.read_double ()
            self [i] ['sn_cnv'] = self.read_double ()
            self [i] ['su_cnv'] = self.read_double ()
            self [i] ['g1_cnv'] = self.read_double ()
            self [i] ['g2_cnv'] = self.read_double ()
            self [i] ['s1_cnv'] = self.read_double ()
            self [i] ['s2_cnv'] = self.read_double ()

            self [i] ['aw_cnv'] = self.read_double ()
            self [i] ['summ0'] = self.read_double ()
            self [i] ['radius0'] = self.read_double ()
            self [i] ['an_cnv'] = self.read_double ()

            self [i] ['abun_cnv'] = self.read_doubles (20)

            self [i] ['eni_cnv'] = self.read_double ()
            self [i] ['enk_cnv'] = self.read_double ()
            self [i] ['enp_cnv'] = self.read_double ()
            self [i] ['ent_cnv'] = self.read_double ()
            self [i] ['epro_cnv'] = self.read_double ()
            self [i] ['enn_cnv'] = self.read_double ()
            self [i] ['enr_cnv'] = self.read_double ()
            self [i] ['ensc_cnv'] = self.read_double ()
            self [i] ['enes_cnv'] = self.read_double ()
            self [i] ['enc_cnv'] = self.read_double ()
            self [i] ['enpist_cnv'] = self.read_double ()
            self [i] ['enid_cnv'] = self.read_double ()
            self [i] ['enkd_cnv'] = self.read_double ()
            self [i] ['enpd_cnv'] = self.read_double ()
            self [i] ['entd_cnv'] = self.read_double ()
            self [i] ['eprod_cnv'] = self.read_double ()
            self [i] ['xlumn_cnv'] = self.read_double ()
            self [i] ['enrd_cnv'] = self.read_double ()
            self [i] ['enscd_cnv'] = self.read_double ()
            self [i] ['enesd_cnv'] = self.read_double ()
            self [i] ['encd_cnv'] = self.read_double ()
            self [i] ['enpistd_cnv'] = self.read_double ()
            self [i] ['xlum_cnv'] = self.read_double ()
            self [i] ['xlum0_cnv'] = self.read_double ()
            self [i] ['entloss_cnv'] = self.read_double ()
            self [i] ['eniloss_cnv'] = self.read_double ()
            self [i] ['enkloss_cnv'] = self.read_double ()
            self [i] ['enploss_cnv'] = self.read_double ()
            self [i] ['enrloss_cnv'] = self.read_double ()
            self [i] ['angit_cnv'] = self.read_double ()
            self [i] ['anglt_cnv'] = self.read_double ()
            self [i] ['xmacc_cnv'] = self.read_double ()
            self.read_int ()
            i += 1
            if max_model != None and i > max_model:
                return
        
    def read_int (self, kind_length = 8):
        if kind_length == 8:
            return struct.unpack ('>' + str (1) + 'i', self.file.read (4)) [0]
        if kind_length == 4:
            return struct.unpack ('>' + str (1) + 'h', self.file.read (2)) [0]
        if kind_length == 2:
            return struct.unpack ('>' + str (1) + 'b', self.file.read (1)) [0]
        raise TypeError ("Unrecognized integer kind.")

    def read_double (self, kind = 8):
        return struct.unpack ('>' + str (1) + 'd', self.file.read (8)) [0]
    
    def read_ints (self, num, kind_length = 8):
        if kind_length == 8:
            return list (struct.unpack ('>' + str (num) + 'i', self.file.read (num * 4)))
        if kind_length == 4:
            return list (struct.unpack ('>' + str (num) + 'h', self.file.read (num * 2)))
        if kind_length == 2:
            return list (struct.unpack ('>' + str (num) + 'b', self.file.read (num)))
        else:
            raise TypeError ("Unrecognized integer kind.")
    
    def read_doubles (self, num, kind = 8):
        x = self.file.read (num * kind)
        return list (struct.unpack ('>' + str (num) + 'd', x))
    
    def read_chars (self, num):
        return self.file.read (num).decode ('utf-8')
    
    def __getitem__ (self, index):
        if isinstance (index, str):
            if index in self.units:
                units = self.units [index]
            else:
                units = "1"
            x = [u.Quantity (model [index], units) for model in self.models]
            try:
                x = u.Quantity (x)
            except TypeError:
                pass
            return x
        return self.models [index]
    
    def __iter__ (self):
        return iter (self.models)
        
    def __len__ (self):
        return len (self.models)
        