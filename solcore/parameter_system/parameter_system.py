import ast  # completely safe eval, but can't call functions
import itertools  # iterate over element/fractions in material string parser
import math  # hyperbolic functions etc in parameterisation
import re  # parsing material string
from functools import lru_cache  # cache function calls to stop things taking forever / recalculating smae things
import numpy
from typing import Optional, Callable

import solcore
from solcore import siUnits
from solcore.source_managed_class import SourceManagedClass


def safe_cacher(maxsize):
    def safewrap(uncached):
        cached = lru_cache(maxsize=maxsize)(uncached)

        def mucked_up_func(*arg, **kwarg):
            try:
                return cached(*arg, **kwarg)
            except:
                return uncached(*arg, **kwarg)

        mucked_up_func.cache_clear = cached.cache_clear
        return mucked_up_func

    return safewrap


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)  # python2: izip_longest


def bow(parent_0_value, parent_1_value, bowing_parameter, x):
    return parent_0_value * (1 - x) + parent_1_value * x - bowing_parameter * (1 - x) * x


class ParameterSystem(SourceManagedClass):
    """Parameter database/bowing system for compound materials.
    
    Once instantiated, this plugin loads the materials parameterisations Parameters for
    compound materials can be retrieved with the get_parameter function.
    """

    def __init__(self, sources: Optional[Callable] = None):
        super().__init__({k: sources(k) for k in sources()})
        # create a dictionary that's safe for the eval function to use, so that config
        # files don't have access to all of python
        self.__assemble_builtins()
        # Matches capital letter + n * small letter, e.g.: In, Ga, As
        self.element_RE = re.compile(
            "([A-Z][a-z]*)")

    def get_parameter(self, material, parameter, verbose=False, **others):
        """Calculate/look up parameters for materials, returns in SI units
        
        Usage: .get_parameter(material_name, parameter_name, **kwargs)
        - material_name is a string of element symbols/fractions, e.g.: In0.2GaAsP0.1
        - parameter_name is a string of 
        - **kwargs captures parameters that may be necessary for some calculations, eg. Temperature
            material fractions may also be specified here, e.g.: .get_parameter("InGaAs", "band_gap", In=0.2)
        
        If a compound material is bowed between two parent materials, the parent materials' parameters are calculated
        recursively with this function. The final parameter is calculated as:
            result=parent_0_value * (1-x) + parent_1_value*x - bowing_parameter * (1-x) * x
        
        The function is cached, so that multiple calls with the same parameters do not incur additional overhead.
        
        """
        material, relevant_parameters = self.__parse_material_string(material, others)

        def tryget(p, alternative):
            try:
                result = self.get_parameter(material, p, verbose=verbose, **relevant_parameters)
                return result
            except ValueError:
                return alternative

        calculation_environment = {
            "get": lambda p: self.get_parameter(material, p, verbose=verbose, **relevant_parameters),
            "tget": tryget

        }

        if verbose:
            print(material, parameter, relevant_parameters)

        assert material in self.database.sections(), "Material {} not in database".format(material)

        if "Final Calculables" in self.database.sections() and \
                        parameter in self.database.options(
                    "Final Calculables") and not parameter in self.database.options(
            material):  # since this is recursive, "Final" gets done first!
            expression = self.database.get("Final Calculables", parameter)
            calculation_environment.update(relevant_parameters)
            result = self.__eval_string_expression(expression, **calculation_environment)
            return result

        if "x" in self.database.options(material):  # material is bowed
            bowed_element = self.database.get(material, "x")
            x = relevant_parameters[bowed_element]
            # del relevant_parameters[bowed_element] # not propagating the element reduces cache misses
            parent0 = self.database.get(material, "parent0")
            parent1 = self.database.get(material, "parent1")

            if parameter in self.database.options(material):
                bowing_parameter = self.__eval_string_expression(
                    self.database.get(material, parameter),
                    **relevant_parameters)
            else:
                bowing_parameter = 0
            del relevant_parameters[bowed_element]  # not propagating the element reduces cache misses

            parent0_value = self.get_parameter(parent0, parameter, verbose=verbose, **relevant_parameters)
            parent1_value = self.get_parameter(parent1, parameter, verbose=verbose, **relevant_parameters)
            return bow(parent0_value, parent1_value, bowing_parameter, x)

        if parameter in self.database.options(material):
            return self.__eval_string_expression(self.database.get(material, parameter), **relevant_parameters)

        if "Immediate Calculables" in self.database.sections() and \
                        parameter in self.database.options("Immediate Calculables"):
            expression = self.database.get("Immediate Calculables", parameter)
            calculation_environment.update(relevant_parameters)
            result = self.__eval_string_expression(expression, **calculation_environment)
            return result

        raise ValueError(
            "Parameter '{}' not in material '{}', nor in calculable parameters.".format(parameter, material))

    def __parse_material_string(self, material_string, other_parameters):
        """parses the material identifier strings of these types:
        
            - In0.2GaAsP0.01
            - InGaAsP {"In":0.2, "P":0.01}
            
            into:
                tuple("InGaAsP", {"In":0.2, "P":0.01})
            
            other parameters are passed into the fractions dictionary. Chemical element Symbols are permitted as 
            sub-material strings, as well as longer words as long as they begin with a capital letter. 
        """
        if "{" in material_string:  # fractions given as a dictionary: InGaAsP {'In':0.2, 'P':0.01}
            identifier, arguments = material_string.split(" ")
            arguments = ast.literal_eval(arguments)
            assert type(arguments) == dict, "{} is not a dict".format(arguments)
            arguments.update(other_parameters)
            return identifier, arguments

        # else: fractions given in parameter or in string: In0.2GaAsP0.01
        elements_and_fractions = self.element_RE.split(material_string)[1:]
        arguments = {}
        for element, fraction in grouper(elements_and_fractions, 2):
            try:
                arguments[element] = float(fraction)
            except:
                pass
        arguments.update(other_parameters)
        # print ("".join(elements_and_fractions[::2]),arguments )
        return "".join(elements_and_fractions[::2]), arguments

    def __eval_string_expression(self, string_expression, **others):
        if " " in string_expression:  # treat second part as unit!
            string_expression, units = string_expression.split(" ", 1)
            use_units = True
        else:
            use_units = False

        if 'T' in string_expression and 'T' not in others.keys():
            raise KeyError('The temperature is needed to calculate this parameter. '
                             'Include keyword argument "T" when calling "get_parameter"')

        non_converted_unit = eval(string_expression, {"__builtins__": self.builtins_replacement}, others)
        in_si_units = siUnits(non_converted_unit, units) if use_units else non_converted_unit
        return in_si_units

    def __assemble_builtins(self):
        self.builtins_replacement = {"max": numpy.max, "min": numpy.min}
        self.builtins_replacement.update(math.__dict__)


if __name__ == "__main__":
    import os

    v = ParameterSystem()
    v.add_source("v", os.path.split(__file__)[0] + "/plugins/vurgaftman/builtins/endpoints.txt")
    v.add_source("v2", os.path.split(__file__)[0] + "/plugins/vurgaftman/builtins/bowing_tree.txt")
    print(solcore.asUnit(v.get_parameter("GaAsSb.5", "band_gap", verbose=True, T=300), "eV"))
    print((v.get_parameter("GaAsSb.75", "band_gap", verbose=True, T=300), "eV"))
    print((v.get_parameter("GaAsSb.5", "band_gap", verbose=True, T=300), "eV"))
