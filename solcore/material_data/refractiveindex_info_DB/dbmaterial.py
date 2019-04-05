# from: https://github.com/HugoGuillen/refractiveindex.info-sqlite

import os
import yaml
import sys
import argparse
import numpy
import scipy.interpolate

class DBMaterial:
    """ Material class"""
    def __init__(self, filename,interpolation_points=100,empty=False):
        """

        :param filename:
        """
        self.refractiveIndex = None
        self.extinctionCoefficient = None
        self.points = interpolation_points
        if empty:
            return

        f = open(filename)
        try:
            material = yaml.safe_load(f)
        except:
            raise Exception('Bad Material YAML File.')
        finally:
            f.close()

        previous_formula = False
        for data in material['DATA']:
            if (data['type'].split())[0] == 'tabulated':
                rows = data['data'].split('\n')
                splitrows = [c.split() for c in rows]
                wavelengths = []
                n = []
                k = []
                for s in splitrows:
                    if len(s) > 0:
                        wavelengths.append(float(s[0]))
                        n.append(float(s[1]))
                        if len(s) > 2:
                            k.append(float(s[2]))
                self.points = len(wavelengths)

                if (data['type'].split())[1] == 'n':

                    if self.refractiveIndex is not None:
                        Exception('Bad Material YAML File')

                    self.refractiveIndex = RefractiveIndexData.SetupRefractiveIndex(formula=-1,
                                                                                    wavelengths=wavelengths,
                                                                                    values=n)
                elif (data['type'].split())[1] == 'k':

                    self.extinctionCoefficient = ExtinctionCoefficientData.SetupExtinctionCoefficient(wavelengths, n)

                    if previous_formula:
                        self.refractiveIndex = RefractiveIndexData.SetupRefractiveIndex(formula=formula,
                                                                                        rangeMin=rangeMin,
                                                                                        rangeMax=rangeMax,
                                                                                        coefficients=coefficients,
                                                                                        interpolation_points=self.points)
                elif (data['type'].split())[1] == 'nk':

                    if self.refractiveIndex is not None:
                        Exception('Bad Material YAML File')

                    self.refractiveIndex = RefractiveIndexData.SetupRefractiveIndex(formula=-1,
                                                                                    wavelengths=wavelengths,
                                                                                    values=n)
                    self.extinctionCoefficient = ExtinctionCoefficientData.SetupExtinctionCoefficient(wavelengths, k)
            elif (data['type'].split())[0] == 'formula':

                if self.refractiveIndex is not None:
                    Exception('Bad Material YAML File')

                formula = int((data['type'].split())[1])
                coefficients = [float(s) for s in data['coefficients'].split()]

                rangeMin = float(data['wavelength_range'].split()[0])
                rangeMax = float(data['wavelength_range'].split()[1])

                previous_formula = True
                self.refractiveIndex = RefractiveIndexData.SetupRefractiveIndex(formula=formula,
                                                                                rangeMin=rangeMin,
                                                                                rangeMax=rangeMax,
                                                                                coefficients=coefficients,
                                                                                interpolation_points=self.points)
        if self.refractiveIndex is not None:
            self.rangeMin = self.refractiveIndex.rangeMin
            self.rangeMax = self.refractiveIndex.rangeMax
        else:
            self.rangeMin = self.extinctionCoefficient.rangeMin
            self.rangeMax = self.extinctionCoefficient.rangeMax

    def get_refractiveindex(self, wavelength):
        """

        :param wavelength:
        :return: :raise Exception:
        """
        if self.refractiveIndex is None:
            raise Exception('No refractive index specified for this material')
        else:
            return self.refractiveIndex.get_refractiveindex(wavelength)

    def get_extinctioncoefficient(self, wavelength):
        """

        :param wavelength:
        :return: :raise NoExtinctionCoefficient:
        """
        if self.extinctionCoefficient is None:
            raise NoExtinctionCoefficient('No extinction coefficient specified for this material')
        else:
            return self.extinctionCoefficient.get_extinction_coefficient(wavelength)

    def get_complete_extinction(self):
        if self.has_extinction():
            return self.extinctionCoefficient.get_complete_extinction()
        else:
            return None

    def get_complete_refractive(self):
        if self.has_refractive():
            return self.refractiveIndex.get_complete_refractive()
        else:
            return None

    def has_refractive(self):
        return self.refractiveIndex is not None

    def has_extinction(self):
        return self.extinctionCoefficient is not None

    def get_page_info(self):
        return self.pageinfo

    def to_csv(self, output):
        refr = self.get_complete_refractive()
        ext = self.get_complete_extinction()
        #FizzFuzz
        if self.has_refractive() and self.has_extinction() and len(refr) == len(ext):
            header = "wl,n,k\n"
            output_f = open(output.replace(".csv","(nk).csv"),'w')
            output_f.write(header)
            for i in range(len(refr)):
                output_f.write(",".join(list(map(str,[refr[i][0],refr[i][1],ext[i][1]])))+"\n")
            output_f.close()
            print("Wrote",output.replace(".csv","(nk).csv"))
        else:
            if self.has_refractive():
                output_f = open(output.replace(".csv","(n).csv"),'w')
                header = "wl,n\n"
                output_f.write(header)
                for i in range(len(refr)):
                    output_f.write(",".join(list(map(str,[refr[i][0],refr[i][1]])))+"\n")
                output_f.close()
                print("Wrote",output.replace(".csv","(n).csv"))
            if self.has_extinction():
                output_f = open(output.replace(".csv","(k).csv"),'w')
                header = "wl,k\n"
                output_f.write(header)
                for i in range(len(ext)):
                    output_f.write(",".join(list(map(str,[ext[i][0],ext[i][1]])))+"\n")
                output_f.close()
                print("Wrote",output.replace(".csv","(k).csv"))

    def to_txt(self, output):
        refr = self.get_complete_refractive()
        ext = self.get_complete_extinction()
        #FizzFuzz

        if self.has_refractive():
            output_f = open(output.replace(".txt","_n.txt"),'w')
            for i in range(len(refr)):
                output_f.write(" ".join(list(map(str,[refr[i][0]*1e-6,refr[i][1]])))+"\n")
            output_f.close()
            print("Wrote",output.replace(".txt","_n.txt"))
        else:
            print("No n data available, file not written")

        if self.has_extinction():
            output_f = open(output.replace(".txt","_k.txt"),'w')
            for i in range(len(ext)):
                output_f.write(" ".join(list(map(str,[ext[i][0]*1e-6,ext[i][1]])))+"\n")
            output_f.close()
            print("Wrote",output.replace(".txt","_k.txt"))
        else:
            print("No k data available, file not written")

    @staticmethod
    def FromLists(pageinfo,wavelengths_r=None,refractive=None,wavelengths_e=None,extinction=None):
        mat = DBMaterial("",empty=True)
        mat.pageinfo = pageinfo
        if refractive is not None:
            mat.refractiveIndex = TabulatedRefractiveIndexData.FromLists(wavelengths_r,refractive)
            mat.rangeMin = mat.refractiveIndex.rangeMin
            mat.rangeMax = mat.refractiveIndex.rangeMax
        if extinction is not None:
            mat.extinctionCoefficient = ExtinctionCoefficientData.FromLists(wavelengths_e,extinction)
            mat.rangeMin = mat.extinctionCoefficient.rangeMin
            mat.rangeMax = mat.extinctionCoefficient.rangeMax
        return mat
#
# Refractive Index
#
class RefractiveIndexData:
    """Abstract RefractiveIndex class"""

    @staticmethod
    def SetupRefractiveIndex(formula, **kwargs):
        """

        :param formula:
        :param kwargs:
        :return: :raise Exception:
        """
        if formula >= 0:
            return FormulaRefractiveIndexData(formula, **kwargs)
        elif formula == -1:
            return TabulatedRefractiveIndexData(**kwargs)
        else:
            raise Exception('Bad RefractiveIndex data type')

    def get_refractiveindex(self, wavelength):
        """

        :param wavelength:
        :raise NotImplementedError:
        """
        raise NotImplementedError('Different for functionally and experimentally defined materials')


class FormulaRefractiveIndexData:
    """Formula RefractiveIndex class"""

    def __init__(self, formula, rangeMin, rangeMax, coefficients,interpolation_points):
        """

        :param formula:
        :param rangeMin:
        :param rangeMax:
        :param coefficients:
        """
        self.formula = formula
        self.rangeMin = rangeMin
        self.rangeMax = rangeMax
        self.coefficients = coefficients
        self.interpolation_points = interpolation_points


    def get_complete_refractive(self):
        #print(self.rangeMin, self.rangeMax)
        wavelength = numpy.linspace(self.rangeMin, self.rangeMax,num=self.interpolation_points)
        extlist = [[wavelength[i], self.get_refractiveindex(wavelength[i] * 1000)] for i in range(len(wavelength))]
        #return numpy.array(extlist)
        return extlist

    def get_refractiveindex(self, wavelength):
        """

        :param wavelength:
        :return: :raise Exception:
        """
        wavelength /= 1000.0
        if self.rangeMin <= wavelength <= self.rangeMax:
            formula_type = self.formula
            coefficients = self.coefficients
            n = 0
            if formula_type == 1:  # Sellmeier
                nsq = 1 + coefficients[0]
                g = lambda c1, c2, w: c1 * (w ** 2) / (w ** 2 - c2 ** 2)
                for i in range(1, len(coefficients), 2):
                    nsq += g(coefficients[i], coefficients[i + 1], wavelength)
                n = numpy.sqrt(nsq)
            elif formula_type == 2:  # Sellmeier-2
                nsq = 1 + coefficients[0]
                g = lambda c1, c2, w: c1 * (w ** 2) / (w ** 2 - c2)
                for i in range(1, len(coefficients), 2):
                    nsq += g(coefficients[i], coefficients[i + 1], wavelength)
                n = numpy.sqrt(nsq)
            elif formula_type == 3:  # Polynomal
                g = lambda c1, c2, w: c1 * w ** c2
                nsq = coefficients[0]
                for i in range(1, len(coefficients), 2):
                    nsq += g(coefficients[i], coefficients[i + 1], wavelength)
                n = numpy.sqrt(nsq)
            elif formula_type == 4:  # RefractiveIndex.INFO
                g = lambda wl, ci, cj, ck, cl: ci * wl**cj / (wl**2 - ck**cl)
                n = coefficients[0]
                n += g(wavelength, *coefficients[1:5])
                n += g(wavelength, *coefficients[5:9])
                for kk in range(len(coefficients[9:]) // 2):
                    n += coefficients[9+kk] * wavelength**coefficients[9+kk+1]

                n = numpy.sqrt(n)
            elif formula_type == 5:  # Cauchy
                g = lambda c1, c2, w: c1 * w ** c2
                n = coefficients[0]
                for i in range(1, len(coefficients), 2):
                    n += g(coefficients[i], coefficients[i + 1], wavelength)
            elif formula_type == 6:  # Gasses
                n = 1 + coefficients[0]
                g = lambda c1, c2, w: c1 / (c2 - w ** (-2))
                for i in range(1, len(coefficients), 2):
                    n += g(coefficients[i], coefficients[i + 1], wavelength)
            elif formula_type == 7:  # Herzberger
                n = coefficients[0]
                n += coefficients[1] / (wavelength**2 - 0.028)
                n += coefficients[2] * (1 / (wavelength**2 - 0.028))**2
                for i, cc in enumerate(coefficients[3:]):
                    n += cc * wavelength**(2*(i+1))
            elif formula_type == 8:  # Retro
                n = coefficients[0]
                n += coefficients[1] * wavelength**2 / (wavelength**2 - coefficients[2])
                n += coefficients[3] * wavelength**2
                n = numpy.sqrt(-(2 * n + 1) / (n - 1))
            elif formula_type == 9:  # Exotic
                n = coefficients[0]
                n += coefficients[1] / (wavelength**2 - coefficients[2])
                n += coefficients[3] * (wavelength - coefficients[4]) / \
                     ((wavelength - coefficients[4])**2 + coefficients[5])
                n = numpy.sqrt(n)
            else:
                raise Exception('Bad formula type')

            return n
        else:
            raise Exception(
                'Wavelength {} is out of bounds. Correct range(um): ({}, {})'.format(wavelength, self.rangeMin,
                                                                                     self.rangeMax))


class TabulatedRefractiveIndexData:
    """Tabulated RefractiveIndex class"""

    def __init__(self, wavelengths, values):
        """

        :param wavelengths:
        :param values:
        """
        self.rangeMin = numpy.min(wavelengths)
        self.rangeMax = numpy.max(wavelengths)

        if self.rangeMin == self.rangeMax:
            self.refractiveFunction = values[0]
        else:
            self.refractiveFunction = scipy.interpolate.interp1d(wavelengths, values)

        self.wavelengths = wavelengths
        self.coefficients = values

    @staticmethod
    def FromLists(wavelengths,values):
        return TabulatedRefractiveIndexData(wavelengths,values)

    def get_refractiveindex(self, wavelength):
        """

        :param wavelength:
        :return: :raise Exception:
        """
        wavelength /= 1000.0
        if self.rangeMin == self.rangeMax and self.rangeMin == wavelength:
            return self.refractiveFunction
        elif self.rangeMin <= wavelength <= self.rangeMax and self.rangeMin != self.rangeMax:
            return self.refractiveFunction(wavelength)
        else:
            raise Exception(
                'Wavelength {} is out of bounds. Correct range(um): ({}, {})'.format(wavelength, self.rangeMin,
                                                                                     self.rangeMax))
    def get_complete_refractive(self):
        extlist =  [[self.wavelengths[i],self.coefficients[i]] for i in range(len(self.wavelengths))]
        #return numpy.array(extlist)
        return extlist


#
# Extinction Coefficient
#
class ExtinctionCoefficientData:
    """ExtinctionCofficient class"""

    @staticmethod
    def SetupExtinctionCoefficient(wavelengths, values):
        """

        :param wavelengths:
        :param values:
        :return:
        """
        return ExtinctionCoefficientData(wavelengths, values)

    @staticmethod
    def FromLists(wavelengths,values):
        return ExtinctionCoefficientData(wavelengths, values)

    def __init__(self, wavelengths, coefficients):
        """

        :param wavelengths:
        :param coefficients:
        """
        self.extCoeffFunction = scipy.interpolate.interp1d(wavelengths, coefficients)
        self.rangeMin = numpy.min(wavelengths)
        self.rangeMax = numpy.max(wavelengths)
        self.wavelengths = wavelengths
        self.coefficients = coefficients

    def get_extinction_coefficient(self, wavelength):
        """

        :param wavelength:
        :return: :raise Exception:
        """
        wavelength /= 1000.0
        if self.rangeMin <= wavelength <= self.rangeMax:
            return self.extCoeffFunction(wavelength)
        else:
            raise Exception(
                'Wavelength {} is out of bounds. Correct range(um): ({}, {})'.format(wavelength, self.rangeMin,
                                                                                     self.rangeMax))
    def get_complete_extinction(self):
        extlist =  [[self.wavelengths[i],self.coefficients[i]] for i in range(len(self.wavelengths))]
        #return numpy.array(extlist)
        return extlist

#
# Custom Exceptions
#
class FormulaNotImplemented(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NoExtinctionCoefficient(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
