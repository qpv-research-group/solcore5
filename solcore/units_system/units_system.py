from collections import defaultdict
import re  # tools to manage regular expresions
import numpy as np
from typing import Optional, Callable

import solcore
from solcore.source_managed_class import SourceManagedClass
from solcore.constants import *


class UnitError(Exception):
    def __init__(self, msg):
        BaseException.__init__(self, msg)


class WrongDimensionError(Exception):
    def __init__(self, msg):
        BaseException.__init__(self, msg)


def generateConversionDictForSISuffix(suffix, centi=False, deci=False, non_base_si_factor=1):
    prefixes = "Y,Z,E,P,T,G,M,k,,m,u,n,p,f,a,z,y".split(",")
    exponents = list(range(8, -9, -1))

    if centi:
        prefixes.append("c")
        exponents.append(-2. / 3.)

    if deci:
        prefixes.append("d")
        exponents.append(-1. / 3.)

    unitNames = ["%s%s" % (prefix, suffix) for prefix in prefixes]
    conversion = [1000. ** exponent * non_base_si_factor for exponent in exponents]

    return dict(zip(unitNames, conversion))


class UnitsSystem(SourceManagedClass):
    """ Contains all the functions related with the conversion of units. While defined inside this class, most of these
    functions are available outside it, being decorated with 'breakout' (see 'Singleton')"""

    def __init__(self, sources: Optional[Callable] = None):
        super().__init__({k: sources(k) for k in sources()})
        self.separate_value_and_unit_RE = re.compile(u"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)(?:[ \t]*(.*))?")
        self.split_units_RE = re.compile(u"(?:([^ \+\-\^\.0-9]+)[\^]?([\-\+]?[^ \-\+]*)?)")
        self.siConversions = {}
        self.dimensions = defaultdict(dict)
        self.read()

    def read(self, source=None, value=None):
        """ Reads the units file and creates a database with all units and conversion factors. """
        if source is not None:
            super().read(source, value)

        for dimension in self.database.sections():
            units = self.database.options(dimension)
            for unit in units:
                if "META_GenerateConversions" in unit:
                    expression = self.database.get(dimension, unit)
                    si_base_unit = expression.split()[0]
                    centi = "centi" in expression
                    deci = "deci" in expression
                    non_base_si_factor = self.safe_eval(expression.split()[1]) if "altbase" in expression else 1
                    dimension_conversions = generateConversionDictForSISuffix(
                        si_base_unit,
                        centi=centi,
                        deci=deci,
                        non_base_si_factor=non_base_si_factor
                    )
                    self.siConversions.update(dimension_conversions)
                    self.dimensions[dimension].update(dimension_conversions)
                    continue
                string_expression = self.database.get(dimension, unit)
                self.siConversions[unit] = self.safe_eval(string_expression)
                self.dimensions[dimension][unit] = self.siConversions[unit]

    def safe_eval(self, string_expression):
        return eval(string_expression, {"__builtins__": {}}, {"constants": solcore.constants})

    def siUnits(self, value, unit):
        """ Convert value from unit to equivalent si-unit

        >>> print(siUnits(1,"mm")) # yields meters
        0.001
        >>> print(siUnits(1,"um")) # yields meters
        1e-06

        :param value: the value to convert
        :param unit: the units of the value
        :return: the value expresed in SI units
        """
        if unit is None or value is None:
            return value

        units_list = self.split_units_RE.findall(unit)
        for unit, power in units_list:
            power = float(power) if power != '' else 1
            value = value * np.power((self.siConversions[unit]),
                                     power)  ### caution, *= is WRONG because it modifies original obj. DO NOT WANT

        return value

    def asUnit(self, value, unit):
        """ Converts from si unit to other unit. It is the reversed of siUnits function

        >>> print(asUnit(1, "mA")) # print 1A in mA.
        1000.0

        :param value: the value to convert, assumed in SI units
        :param unit: the new units
        :param dimension: the value expressed in the new units.
        :return:
        """
        if unit is None or value is None:
            return value

        units_list = self.split_units_RE.findall(unit)
        for unit, power in units_list:
            power = float(power) if power != '' else 1
            value = value / (self.siConversions[unit]) ** power  ### caution, /= is WRONG because it modifies original obj. DO NOT WANT

        return value

    def si(self, *args):
        """ Utility function that forwards to either siUnit or siUnitFromString"""
        if type(args[0]) == str:
            return self.siUnitFromString(*args)
        return self.siUnits(*args)

    def siUnitFromString(self, string):
        """ Converts a string of a number with units into si units of that quantity

        >>> print(si("5 mm s-1")) # output in m/s
        0.005
        >>> print(si("5e-0mm-2")) # output in m2
        5000000.0
        >>> print(si("5"))
        5.0

        :param string: the string to convert
        :return: the value in SI units
        """
        # if unit is None or value is None:
        #     return value

        matchObj = self.separate_value_and_unit_RE.match(string)
        value, unit = matchObj.groups()
        value = float(value)
        units_list = self.split_units_RE.findall(unit)
        for unit, power in units_list:
            power = float(power) if power != '' else 1
            value *= (self.siConversions[unit]) ** power
        return value

    def convert(self, value, from_unit, to_unit):
        """ Converts between comparable units, does NOT check if units are comparable.

        >>> print(convert(1, "nm", "mm"))
        1e-06
        >>> print(convert(1, "um", "nm"))
        1000.0
        >>> print(convert(1, "cm s-1", "km h-1"))
        0.036

        :param value: the value ot convert
        :param from_unit: the original unit
        :param to_unit: the final unit
        :return: the value expressed in the final unit
        """
        return self.asUnit(self.siUnits(value, from_unit), to_unit)

    def eVnm(self, value):
        """ Bi-directional conversion between nm and eV.

        >>> print('%.3f'%eVnm(1000))
        1.240
        >>> print('%i'%round(eVnm(1)))
        1240

        :param value: a number with units [nm] or [eV].
        :return: either the conversion [nm] --> [eV], or [eV] --> [nm]
        """
        factor = self.asUnit(h, "eV") * self.asUnit(c, "nm")
        return factor / value

    def nmJ(self, value):
        """ Bi-directional conversion between nm and J.

        >>> print(nmJ(1000))
        1.9864452126e-19
        >>> print(nmJ(2e-18))
        99.3222606298

        :param value: a number with units [nm] or [J].
        :return: either the conversion [nm] --> [J], or [J] --> [nm]
        """
        factor = h * c
        return factor / self.siUnits(value, "nm")

    def mJ(self, value):
        """ Bi-directional conversion between m and J.

        >>> print(mJ(1000))
        1.986445212595144e-25
        >>> print(mJ(2e-18))
        9.93222606297572e-08

        :param value: a number with units [m] or [J].
        :return: either the conversion [m] --> [J], or [J] --> [m]
        """
        factor = h * c
        return factor / value

    def nmHz(self, value):
        """ Bi-directional conversion between nm and Hz.

        :param value: a number with units [nm] or [Hz].
        :return: Either a number which is the conversion [nm] --> [Hz] or [Hz] --> [nm]
        """
        factor = self.asUnit(c, "nm s-1")
        return factor / value

    def spectral_conversion_nm_ev(self, x, y):
        """ Bi-directional conversion between a spectrum per nanometer and a spectrum per electronvolt.

        Example:
        1) nm --> eV conversion
        wavelength_nm
        photon_flux_per_nm
        energy_ev, photon_flux_per_ev = spectral_conversion_nm_ev(wavelength_nm, photon_flux_per_nm)

        2) eV --> nm conversion
        energy_ev
        photon_flux_per_ev
        wavelength_nm, photon_flux_per_nm = spectral_conversion_nm_ev(energy_ev, photon_flux_per_ev)

        Discussion:
        A physical quantities such as total number of photon in a spectrum or
        total energy of a spectrum should remain invariant after a transformation
        to different units. This is called a spectral conversion. This function
        is bi-directional because the mathematics of the conversion processes
        is symmetrical.

        >>> import numpy as np
        >>> x = np.array([1,2,3])
        >>> y = np.array([1,1,1])
        >>> area_before = np.trapz(y, x=x)
        >>> x_new, y_new = spectral_conversion_nm_ev(x, y)
        >>> area_after = np.trapz(y_new, x=x_new)
        >>> compare_floats(area_before, area_after, relative_precision=0.2)
        True

        :param x: abscissa of the spectrum in units of [nm] or [eV]
        :param y: ordinate of the spectrum in units of [something/nm] or [something/eV]
        :return: A tuple (x, y) which has units either [eV, something/eV] or [nm. something/nm].
        """
        x_prime = self.eVnm(x)
        conversion_constant = self.asUnit(h, "eV s") * self.asUnit(c, "nm s-1")
        y_prime = y * conversion_constant / x_prime ** 2
        y_prime = reverse(y_prime)  # Wavelength ascends as electronvolts decends therefore reverse arrays
        x_prime = reverse(x_prime)
        return (x_prime, y_prime)

    def spectral_conversion_nm_hz(self, x, y):
        """ Bi-directional conversion between a spectrum per nanometer and a spectrum per Hertz.

        Example:
        1) nm --> Hz conversion
        wavelength_nm
        photon_flux_per_nm
        energy_hz, photon_flux_per_hz = spectral_conversion_nm_hz(wavelength_nm, photon_flux_per_nm)

        2) Hz --> nm conversion
        energy_hz
        photon_flux_per_hz
        wavelength_nm, photon_flux_per_nm = spectral_conversion_nm_ev(energy_hz, photon_flux_per_hz)

        Discussion:
        A physical quantities such as total number of photon in a spectrum or
        total energy of a spectrum should remain invariant after a transformation
        to different units. This is called a spectral conversion. This function
        is bi-directional because the mathematics of the conversion processes
        is symmetrical.

        >>> import numpy as np
        >>> x = np.array([1,2,3])
        >>> y = np.array([1,1,1])
        >>> area_before = np.trapz(y, x=x)
        >>> x_new, y_new = spectral_conversion_nm_hz(x, y)
        >>> area_after = np.trapz(y_new, x=x_new)
        >>> compare_floats(area_before, area_after, relative_precision=0.2)
        True

        :param x: abscissa of the spectrum in units of [nm] or [Hz]
        :param y: ordinate of the spectrum in units of [something/nm] or [something/Hz]
        :return: A tuple (x, y) which has units either [eV, something/nm] or [nm. something/Hz].
        """
        x_prime = self.nmHz(x)
        conversion_constant = self.asUnit(c, "nm s-1")
        y_prime = y * conversion_constant / x_prime ** 2
        y_prime = reverse(y_prime)  # Wavelength ascends as frequency decends therefore reverse arrays
        x_prime = reverse(x_prime)
        return (x_prime, y_prime)

    def sensibleUnits(self, value, dimension, precision=2):
        """ Attempt to convert a physical quantity of a particular dimension to the most sensible units

        >>> print(sensibleUnits(0.001,"length",0))
        1 mm
        >>> print(sensibleUnits(1000,"length",0))
        1 km
        >>> print(sensibleUnits(si("0.141 days"),"time", 5))
        3.38400 h

        :param value: The value to re-calculate in SI units
        :param dimension: The dimension of the value. Possible values are: 'luminous intensity', 'pressure', 'time',
         'angle', 'temperature', 'current', 'force', 'charge', 'power', 'voltage', 'resistance', 'mass', 'length', 'energy'
        :param precision: Precission of the converted value.
        :return: A string with the value in the more 'sensible' units and the units.
        """
        negative = ""
        if value < 0:
            value *= -1
            negative = "-"
        formatting = "%s%%.%if %%s" % (negative, precision)
        d = self.dimensions[dimension]
        possibleUnits = d.keys()
        if value == 0:
            return formatting % (0, "")
        allValues = [abs(np.log10(self.asUnit(value, unit))) for unit in possibleUnits]
        bestUnit = possibleUnits[allValues.index(min(allValues))]
        return formatting % (self.asUnit(value, bestUnit), bestUnit)

    def eV(self, e):
        """ Transform an energy value in SI units in a string expresing the value in eV.

        >>> print(eV(1e-19))
        0.624 eV

        :param e: Energy value in SI units
        :return: A string with the energy converted in eV and its units.
        """
        return "%.3f eV" % self.asUnit(e, "eV")

    def guess_dimension(self, unit):
        """ Guess the dimension of a unit.

        >>> print guess_dimension("nm")
        length

        :param unit: the unit.
        :return: None
        """
        possibilities = [key for key in self.dimensions.keys() if unit in self.dimensions[key]]

        assert len(possibilities) != 0, "Guessing dimension of '%s': No candidates found" % unit
        assert len(
            possibilities) == 1, "Guessing dimension of '%s': Multiple candidates found, please convert manually. (%s)" % (
            unit, ", ".join(possibilities))

        return possibilities[0]

    def list_dimensions(self):
        for dim in self.dimensions.keys():
            print(
                "%s: %s" % (dim, ", ".join([k for k in self.dimensions[dim].keys() if k is not None and k is not ""])))


def compare_floats(a, b, absoulte_precision=1e-12, relative_precision=None):
    """ Returns true if the absolute difference between the numbers a and b is less than the precision.

    Arguments:
    a -- a float
    b -- a float

    Keyword Arguments (optional):
    absolute_precision -- the absolute precision, abs(a-b) of the comparison.
    relative_precision -- the relative precision, max(a,b)/min(a,b) - 1. of the comparison.

    Returns:
    True if the numbers are the same within the limits of the precision.
    False if the number are not the same within the limits of the precision.
    """

    if relative_precision is None:
        absolute = abs(a - b)
        if absolute < absoulte_precision:
            return True
        else:
            return False
    else:
        relative = max(a, b) / min(a, b) - 1.
        if relative < relative_precision:
            return True
        else:
            return False


def independent_nm_ev(x, y):
    return UnitsSystem().eVnm(x)[::-1], y[::-1]


def independent_nm_J(x, y):
    return UnitsSystem().nmJ(x)[::-1], y[::-1]


def independent_m_J(x, y):
    return reverse(UnitsSystem().mJ(x)), reverse(y)


def reverse(x):
    return x[::-1]


if __name__ == '__main__':
    UnitsSystem()

    print(solcore.eVnm(1240))
