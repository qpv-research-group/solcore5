# A number of custom functions/ classes than are designed to aid with plotting of graphs in Python.
# Good excuse to test out some OOP in python, T. Wilson 01-06-2016

# Custom Plot Colours :: This uses the colours and hex codes given online at;
# http://www.rapidtables.com/web/color/RGB_Color.htm
# to give access to a whole range of different colours.
# Will write as a class for practice.

import os
import re
from cycler import cycler, Cycler
from typing import List, Tuple, Dict, Union

# define the file path to the Colours.csv file, which sits within the Graphing Package
this_dir, this_fname = os.path.split(__file__)
file_path = os.path.join(this_dir, 'Colours.txt')

# Need to read the file in as text, declare
# opens .txt as file object
raw_txt: List[str] = []

with open(file_path) as file:

    raw_txt = file.read().splitlines()


def colours(colour: str, type: str = 'hex') -> Union[str, Tuple[float, ...]]:
    """colours(colour, type) :: Function identifies the name of the colour and returns the correct HEX or RGB key.
        Default return is in HEX.

    :param colour: string variable, name of colour as listed in the table of colours at;
                    http://www.rapidtables.com/web/color/RGB_Color.htm
    :param type: string variable, choices of "hex or "rgb". If colour is not found an error is returned.
    :return: HEX string or RGB tuple
    """

    # Initialise variables
    Cols: Dict[str, Union[str, Tuple[float, ...]]] = {}
    Cols.setdefault("Colour", colour)
    HEX2DEC: List[str] = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']

    for Lines in raw_txt:
        # Using regular expressions package, matching the material name with the correct section of the file.
        # If the name is found, ID is switched to True.
        split_line: List[str] = Lines.split(" = ")

        if re.fullmatch(colour.lower(), split_line[0]) != None:

            Cols["HEX"] = split_line[1]
            HEX: str = split_line[1]

            R: int = HEX2DEC.index(HEX[1]) * 16 ** 1 + HEX2DEC.index(HEX[2]) * 16 ** 0
            G: int = HEX2DEC.index(HEX[3]) * 16 ** 1 + HEX2DEC.index(HEX[4]) * 16 ** 0
            B: int = HEX2DEC.index(HEX[5]) * 16 ** 1 + HEX2DEC.index(HEX[6]) * 16 ** 0

            Cols["RGB"] = (R/255, G/255, B/255)

            break

        if "$END$" in Lines:
            raise ValueError("Colour name not found... Check http://www.rapidtables.com/web/color/RGB_Color.htm")

    if type.upper() == "HEX":
        return Cols["HEX"]

    elif type.upper() == "RGB":
        return Cols["RGB"]
    else:
        raise ValueError("Invalid 'type' selection... Choose either 'hex' or 'rgb'.")


def colour_cycle(name: str) -> Cycler:
    """colour_cycle(name, type) :: Function returns a cycler instance of a list of colours (in HEX) that can be used
        when plotting multi-plots.

    :param name: string variable matching the list of available colour cycles. Current list includes:
                    'rainbow', 'origin_basic', 'origin_system_colour_list'.
    :return: colour cycle instance
    """

    Arg_Out : List[Union[str, Tuple[float, ...]]] = []

    if name.lower() == 'rainbow': # 29 colours in shades of the rainbow...

        Arg_Out = [colours('Dark Red'), colours('Crimson'), colours('Salmon'),
                   colours('Orange Red'), colours('orange'), colours('Gold'),
                   colours('yellow'), colours('Lawn Green'), colours('Olive Drab'),
                   colours('Green', ), colours('Lime Green'), colours('Sea Green'),
                   colours('Teal'), colours('Cyan'), colours('Turquoise'),
                   colours('Steel Blue'), colours('Deep Sky Blue'), colours('Dodger Blue'),
                   colours('Blue'), colours('Navy'), colours('Indigo'),
                   colours('Slate Blue'), colours('Dark Magenta'), colours('Dark Orchid'),
                   colours('Purple'), colours('Magenta'), colours('Medium Violet Red'),
                   colours('Deep Pink'), colours('Black')]

    elif name.lower() == 'origin_basic': # 23 basic colours as defined in Origin...

        Arg_Out = ["#000000", "#FF0000", "#0000FF", "#FF00FF", "#008000", "#000080", "#8000FF", "#800080",
                   "#800000", "#808000", "#2B63A2", "#1E9696", "#9B641A", "#10C73E", "#B9247A", "#2DC5CC",
                   "#3F4198", "#93AC2B", "#808080", "#966464", "#649664", "#2BA3CA", "#326496"]

    elif name.lower() == "origin_system_color_list": # 13 colours defined by Origin's system colour list...

        Arg_Out = ["#F00820", "#FA3C3C", "#F08228", "#E6AF2D", "#E6DC32", "#A0E632", "#00DC00", "#00D28C",
                   "#00C8C8", "#00A0FF", "#1E3CFF", "#6E00DC", "#A000C8"]

    return cycler('color', Arg_Out)
