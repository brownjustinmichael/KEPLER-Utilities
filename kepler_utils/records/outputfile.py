import fortranfile
import numpy as np
import struct
import re

class OutputFile:
    """
    This class reads a fortran output file into memory.

    :type file_name: :class:`str`
    :param file_name: The location of the CNV file to load
    """
    
    def __init__ (self, file_name, max_model = None, verbose = False):
        self.file_name = file_name
        self.reload ()

    def reload (self):
        self.file = fortranfile.FortranFile (self.file_name)

    def _read_int (self, kind_length = 8):
        """
        Read an integer from self.file

        :type kind_length: :class:`int`
        :param kind_length: The byte length of the integer
        :rtype: :class:`int`
        :return: The integer representatation of the data
        """
        return self._read_ints (1, kind_length) [0]

    def _read_double (self, kind_length = 8):
        """
        Read a real number from self.file

        :type kind_length: :class:`int`
        "param kind_length: The byte length of the real number
        :rtype: :class:`float`
        :return: The float representatation of the read data
        """
        return self._read_doubles (1, kind_length) [0]
    
    def _read_ints (self, num, kind_length = 8):
        """
        Read a list of integers from self.file

        :type num: :class:`int`
        :param num: The number of integers to read
        :type kind_length: :class:`int`
        :param kind_length: The byte length of the integers
        :rtype: :class:`list` of :class:`int`
        :return: The list of integers representatating of the read data
        """
        if kind_length == 16:
            return np.array (struct.unpack ('>' + str (num) + 'q', self.file.read (num * 8)), dtype = int)
        if kind_length == 8:
            return np.array (struct.unpack ('>' + str (num) + 'i', self.file.read (num * 4)), dtype = int)
        if kind_length == 4:
            return np.array (struct.unpack ('>' + str (num) + 'h', self.file.read (num * 2)), dtype = int)
        if kind_length == 2:
            return np.array (struct.unpack ('>' + str (num) + 'b', self.file.read (num)), dtype = int)
        raise TypeError ("Unrecognized integer kind.")
    
    def _read_doubles (self, num, kind = 8):
        """
        Read a list of real numbers from self.file

        :type num: :class:`int`
        :param num: The number of real numbers to read
        :type kind_length: :class:`int`
        :param kind_length: The byte length of the real numbers
        :rtype: :class:`list` of :class:`float`
        :return: The string representatation of the read data
        """
        if kind == 8:
            return np.array (struct.unpack ('>' + str (num) + 'd', self.file.read (num * 8)))
        if kind == 4:
            return np.array (struct.unpack ('>' + str (num) + 'f', self.file.read (num * 4)))
        raise TypeError ("Unrecognized double kind.")

    def _read_chars (self, num, as_byte_string = False):
        """
        Read a string from self.file

        :type num: :class:`int`
        :param num: The number of integers to read
        :rtype: :class:`str`
        :return: The string representatation of the read data
        """
        chars = self.file.read (num)
        if as_byte_string:
            return chars
        chars = re.sub(r'[^ -~]', '', chars.decode ('utf-8'))
        return chars.rstrip ()

    def _read_strings (self, num, chars_per, as_byte_string = False):
        """
        Read num strings from self.file, if asByteString, return a list of byte strings, else return a list of python strings

        :type num: :class:`int`
        :param num: The number of strings to load
        :type charsPer
        """
        return [self._read_chars (chars_per, as_byte_string) for i in range (num)]
        
    def _skip (self, num):
        """
        Skip num bytes in self.file

        :type num: :class:`int`
        :param num: The number of bytes to skip
        """
        self.file.read (num)
        
    