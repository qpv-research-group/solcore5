"""This module gives access to the vast library of optical constant data, made freely available by the SOPRA-SA
optoelectronics company founded in 1948. For further detail on the data and SOPRA-SA see the legacy website:
http://www.sspectra.com/sopra.html"""

import numpy as np
import os, sys
import re
from natsort import natsorted
from configparser import ConfigParser
from solcore.science_tracker import science_reference
from solcore import config, SOLCORE_ROOT

SOPRA_PATH = os.path.abspath(config['Others']['sopra'].replace('SOLCORE_ROOT', SOLCORE_ROOT))
compounds_path = os.path.join(SOPRA_PATH, "compounds.txt")
compounds_info = ConfigParser()
compounds_info.read(compounds_path)


# Defining the SOPRA_DB class variable
class sopra_database:
    """Import the SOPRA_DB module from the solcore.material_system package and get started by selecting a material from
    the extensive list that SOPRA-SA compiled;

    >>> GaAs = sopra_database('GaAs')

    Once imported a number of useful methods can be called to return n, k and alpha data for the desired material.
    """

    science_reference("All optical constant data made avaialble by SOPRA-SA",
                      "http://www.sspectra.com/sopra.html")

    def __init__(self, Material):
        # Define filepath to the SOPRA database for file import...
        self.__SOPRA_PATH = SOPRA_PATH

        # Load in SOPRA_DB.csv database file
        DB = np.genfromtxt(os.path.join(self.__SOPRA_PATH, "SOPRA_DB_Updated.csv"), delimiter=",", dtype=str)

        self.__fname = None
        for fname, symbol, range, info in DB:

            if re.fullmatch(Material.upper(), fname) is not None:
                self.__fname = fname

                self.path = os.path.join(self.__SOPRA_PATH, self.__fname + ".MAT")

                # self.info contains all detail loaded from SOPRA_DB,csv file...
                self.info = {"Material": symbol,
                             "Wavelength (nm)": range,
                             "File Info": info,
                             "File Path": self.path}

        # If the material name is incorrect then the material attribute is not written to
        if self.__fname is None:

            print("SOPRA_DB :: ERROR :: Material not found in SOPRA_DB... Check materials list...")
            print("Similar Matches ::")
            for fname, symbol, range, info in DB:

                if re.match(Material.upper(), fname) is not None:
                    print(fname)

            # If the exception is caught, exit the program as nothing else useful can be done...
            # sys.exit()
            raise SOPRAError("Material not found in SOPRA database: {}".format(Material))

    @staticmethod
    def material_list():
        """ SOPRA_DB.material_list() :: Loads a list (.pdf file) of all available SOPRA materials. """

        print("Opening List of Available Materials in the SOPRA database")

        # Need different treatment depending on computer OS.
        if sys.platform == 'darwin':
            # Find spaces in the filename and add a \ before (for unix based systems)
            directory = SOPRA_PATH.split(" ")

            new_path = directory[0]
            for i in range(1, len(directory), 1):
                new_path = new_path + "\ " + directory[i]

            os.system("open " + os.path.join(new_path, "List_Of_Files_Updated_PDF.pdf"))

        elif sys.platform == 'linux':
            # Find spaces in the filename and add a \ before (for unix based systems)
            directory = SOPRA_PATH.split(" ")

            new_path = directory[0]
            for i in range(1, len(directory), 1):
                new_path = new_path + "\ " + directory[i]

            os.system("xdg-open " + os.path.join(new_path, "List_Of_Files_Updated_PDF.pdf"))

        elif sys.platform == 'win32':
            # Find spaces in the filename and add a \ before (for unix based systems)

            os.system("start " + os.path.join(SOPRA_PATH, "List_Of_Files_Updated_PDF.pdf"))

    def load_n(self, Lambda=None):
        """ SOPRA_DB.load_n(Lambda) :: Load refractive index (n) data of the requested material.
            Optional argument Lambda allows user to specify a custom wavelength range. data will be interpolated into
                this range before output.

            Returns: Tuple of (Wavelength, n)"""

        try:
            os.stat(self.path)
        except FileNotFoundError:
            print('load_n :: WARNING :: There is no individual data file for, ' + self.material + ".")
            print('This material may be part of a set of varying composition, check the materials list...')
            sys.exit()

        # Load in data from file...
        Wav, n = np.genfromtxt(self.path, delimiter="*", skip_header=3, skip_footer=3, usecols=(2, 3), unpack=True)

        if Lambda is not None:
            # Interpolate in range specified by Lambda...
            n_interp = np.interp(Lambda, Wav, n)

            return (Lambda, n_interp)

        else:
            return (Wav, n)

    def load_k(self, Lambda=None):
        """ SOPRA_DB.load_k(Lambda) :: Load refractive index (n) data of the requested material.
            Optional argument Lambda allows user to specify a custom wavelength range. data will be interpolated into
                this range before output.

            Returns: Tuple of (Wavelength, k) """

        try:
            os.stat(self.path)
        except FileNotFoundError:
            print('load_k :: WARNING :: There is no individual data file for, ' + self.material + ".")
            print('This material may be part of a set of varying composition, check the materials list...')
            sys.exit()

        # Load in data from file...
        Wav, k = np.genfromtxt(self.path, delimiter="*", skip_header=3, skip_footer=3, usecols=(2, 4), unpack=True)

        if Lambda is not None:
            # Interpolate in range specified by Lambda...
            k_interp = np.interp(Lambda, Wav, k)

            return (Lambda, k_interp)
        else:
            return (Wav, k)

    def load_alpha(self, Lambda=None):
        """ SOPRA_DB.load_alpha(Lambda) :: Load refractive index (n) data of the requested material.
            Optional argument Lambda allows user to specify a custom wavelength range. data will be interpolated into
                this range before output.

            Returns: Tuple of (Wavelength, alpha) """

        Wav, k = self.load_k(Lambda=Lambda)

        return (Wav, ((4 * np.pi) / (Wav * 1E-9)) * k)

    def load_temperature(self, Lambda, T=300):
        """ SOPRA_DB.load_temperature(T, Lambda) :: Loads n and k data for a set of materials with temperature dependent
                data sets
            Optional argument T defaults to 300K
            Required argument Lambda specifies a wavelength range and the data is interpolated to fit. This is a
                required argument here as not all data sets in a group are the same length (will be fixed in a
                subsequent update).

            Returns: Tuple of (Wavelength, n, k) """

        T_degC = T - 273.15  # Convert from Kelvin to degC (units given in the data)...

        # Navigate to the correct folder that contains temperature dependent data...
        path = os.path.join(self.__SOPRA_PATH, self.__fname + "_T")

        try:
            os.stat(path)

        except FileNotFoundError:
            print("load_temperature :: WARNING :: Material folder does not exists... Check materials list...")
            # If material folder is not found exit program as nothing more useful can be done...
            sys.exit()

        # if folder exists, read in files from folder...
        Folder = natsorted(os.listdir(path))

        DATA = []
        TEMP = []

        for files in Folder:

            # .DS_Store is a metadata file used in Mac OS X to store various file/ folder info. Ignoring...
            if ".DS_Store" not in files:

                # extract temperature from filename...
                Num = re.findall("[-+]?\d+[\.]?\d*", files)
                Num.append("0")
                TEMP.append(float(Num[0]))

                Wav, n, k = np.genfromtxt(os.path.join(path, files), delimiter="*",
                                          skip_header=3, skip_footer=3, usecols=(2, 3, 4), unpack=True)

                if Lambda is not None:
                    # Interpolate if the Lambda argument is specified, if not pass loaded Wav, n and k...
                    n_interp = np.interp(Lambda, Wav, n)
                    k_interp = np.interp(Lambda, Wav, k)

                    DATA.append((Lambda, n_interp, k_interp, float(Num[0])))

                else:
                    DATA.append((Wav, n, k, float(Num[0])))

        # Check and see if the entered temperature is within the range of data...
        if T_degC <= min(TEMP):
            print("load_temperature :: WARNING :: Desired Temperature < than the minimum (%6.1f K)" % (
                min(TEMP) + 273.15))
            print("Returned interpolated data will be that at Tmin = %6.1f K" % (min(TEMP) + 273.15))

        elif T_degC >= max(TEMP):
            print("load_temperature :: WARNING :: Desired Temperature > than the maximum (%6.1f K)" % (
                max(TEMP) + 273.15))
            print("Returned interpolated data will be that at Tmax = %6.1f K" % (max(TEMP) + 273.15))

        # use linear interpolation to interpolate the data at the desired temperature...
        n_interp_data = []
        k_interp_data = []

        # In range of all wavelenghs...
        for i in range(0, len(DATA[0][0])):

            T_list = []
            k_at_T = []
            n_at_T = []

            # At each wavelenth, build a list of n and k at each T...
            for X, n, k, temp in DATA:
                T_list.append(temp)
                k_at_T.append(k[i])
                n_at_T.append(n[i])

            # Interpolate the point corresponding to T and build the new data array...
            n_interp_data.append(np.interp(T - 273.19, T_list, n_at_T))
            k_interp_data.append(np.interp(T - 273.19, T_list, k_at_T))

        # Return the Wavelength vector and the new n and k data...
        return (DATA[0][0], n_interp_data, k_interp_data)

    def load_composition(self, Lambda, **kwargs):
        """ SOPRA_DB.load_temperature(T, Lambda) :: Loads n and k data for a set of materials with varying composition.
            Required argument Lambda specifies a wavelength range and the data is interpolated to fit. This is a
                required argument here as not all data sets in a group are the same length (will be fixed in a
                subsequent update).
            Keyword argument :: Specify the factional material and fraction of the desired alloy.

            Returns: Tuple of (Wavelength, n, k) """
        # Use of keyword args allows the user to specify the required material fraction for neatness...
        for material in kwargs:
            mat_fraction = material
            frac = kwargs[material]

        # Navigate to the correct folder that contains temperature dependent data...
        path = os.path.join(self.__SOPRA_PATH, self.__fname + "_" + mat_fraction.upper())
        try:
            os.stat(path)

        except FileNotFoundError:
            print("load_composition :: WARNING :: Material folder does not exists... Check materials list or check" +
                  " that composition material is correct...")
            # If material folder is not found exit program as nothing more useful can be done...
            sys.exit()

        # if folder exists, read in files from folder...
        Folder = natsorted(os.listdir(path))

        DATA = []
        COMP = []

        for files in Folder:

            # .DS_Store is a metadata file used in Mac OS X to store various file/ folder info. Ignoring...
            if ".DS_Store" not in files:

                # extract temperature from filename...
                Num = re.findall("[-+]?\d+[\.]?\d*", files)
                Num.append("0")
                COMP.append(float(Num[0]))

                Wav, n, k = np.genfromtxt(os.path.join(path, files), delimiter="*",
                                          skip_header=3, skip_footer=3, usecols=(2, 3, 4), unpack=True)

                if Lambda is not None:
                    # Interpolate if the Lambda argument is specified, if not pass loaded Wav, n and k...
                    n_interp = np.interp(Lambda, Wav, n)
                    k_interp = np.interp(Lambda, Wav, k)

                    DATA.append((Lambda, n_interp, k_interp, float(Num[0])))

                else:
                    DATA.append((Wav, n, k, float(Num[0])))

        # Check and see if the entered temperature is within the range of data...
        if frac <= min(COMP):
            print("load_composition :: WARNING :: Desired composition < than the minimum (%6.1f %%)" % min(COMP))
            print("Returned interpolated data will be that at %6.1f %%" % min(COMP))

        elif frac >= max(COMP):
            print("load_composition :: WARNING :: Desired composition > than the maximum (%6.1f %%)" % max(COMP))
            print("Returned interpolated data will be that at %6.1f %%" % max(COMP))

        # use linear interpolation to interpolate the data at the desired temperature...
        n_interp_data = []
        k_interp_data = []

        # In range of all wavelenghs...
        for i in range(0, len(DATA[0][0])):

            x_list = []
            k_at_C = []
            n_at_C = []

            # At each wavelenth, build a list of n and k at each T...
            for X, n, k, x in DATA:
                x_list.append(x)
                k_at_C.append(k[i])
                n_at_C.append(n[i])

            # Interpolate the point corresponding to T and build the new data array...
            n_interp_data.append(np.interp(frac, x_list, n_at_C))
            k_interp_data.append(np.interp(frac, x_list, k_at_C))

        # Return the Wavelength vector and the new n and k data...
        return (DATA[0][0], n_interp_data, k_interp_data)


class SOPRAError(Exception):
    def __init__(self, message):
        self.message = message
