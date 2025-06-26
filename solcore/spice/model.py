"""The classes help to create a distributed SPICE model for a solar cell. 
Think of these building blocks as sub-circuits that you can combine into a
3D structure to build the solar cell.
"""


class Header:
    """A class representing header information for the SPICE netlist."""

    def __init__(
            self,
            temperature : float = 300.0
        ):
        """
        Parameters
        ----------
        temperature : float (Optional)
            The temperature in Kelvin
        """
        self.temperature = temperature

    def netlist(self):
        return f"""
        * HEADER
        .options TNOM={self.temperature - 273.15} TEMP={self.temperature - 273.15}
        vin in 0 DC 0
        """


class Diodes:
    """A class representing header information for the SPICE file that contains diode models"""

    def __init__(
            self,
            Eg: float,
            j01: float,
            j02: float,
            area: float,
            n1: float = 1.0,
            n2: float = 2.0,
            junction_idx: int = 0
        ):
        """
        Parameters
        ----------
        Eg : float
            Band gap in Joules.
        j01 : float
            Saturation current density in neutral region in amps per metre sq.
        j02 : float
            Saturation current density in the bulk region in amps per metre sq.
        area : float
            The top surface area of this segment
        n1 : float (Optional)
            Ideality factor, default n1 = 1
        n2 : float (Optional)
            Ideality factor, default n2 = 2
        """
        self.Eg = Eg
        self.j01 = j01
        self.j02 = j02
        self.area = area
        self.n1 = n1
        self.n2 = n2
        self.junction_idx = int(junction_idx)

    def netlist(self):
        return f"""
        * HEADER
        .model __D1_{self.junction_idx} D(is={self.area * self.j01},n={self.n1},eg={self.Eg})
        .model __D2_{self.junction_idx} D(is={self.area * self.j02},n={self.n2},eg={self.Eg})
        """


class Metal:
    """A unit cell representing a SPICE model of metal grid finger.

    The resistance of the cell is calculated using the geometric
    properties of the contact and the resistivity of the metal.
    """

    def __init__(
            self,
            idx,
            metal_height: float,
            x_length: float,
            y_length: float,
            resistivity_metal: float,
            resistivity_contact: float
        ):
        """
        Parameters
        ----------
        idx : Tuple
            The 3D index of this element in the grid e.g. (1, 2, 3)
        metal_height : float
            The metalisation height of the grid finger in meters for this cell
        x_length  : float
            The metalisation width of the grid finger in meters for this cell. This
            depends on the resolution of the discretiation in the X direction.
        y_length : float
            The metalisation width of the grid finger in meters for this cell. This
            depends on the resolution of the discretiation in the Y direction.
        resistivity_metal : float
            The resistivity of the metal in Ohm meter for this cell.
        resistivity_contact : float
            The resistivity of the metal-semiconductor contact in Ohm meter^2 for
            this cell.
        """
        self.idx = idx
        self.metal_height = metal_height
        self.x_length = x_length
        self.y_length = y_length
        self.resistivity_metal = resistivity_metal
        self.resistivity_contact = resistivity_contact
        self.left = f"NX_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.near = f"NY_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.top = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.right = f"NX_{self.idx[0]+1}_{self.idx[1]}_{self.idx[2]}"
        self.far = f"NY_{self.idx[0]}_{self.idx[1]+1}_{self.idx[2]}"
        self.bottom = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]+1}"
        self.centre = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.element = f"{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"

    def netlist(self):
        """Return the SPICE netlist for this cell.

        Parameters
        ----------
        idx : Tuple
            The (i, j, k) indexes of this cell's location in the grid. 
        """
        r_contact = self.resistivity_contact / (self.x_length * self.y_length) # usually ~ 5 mOhms
        r_metal = self.resistivity_metal * self.y_length / (self.metal_height * self.x_length) # usually ~ 1 Ohm

        return f"""
        * METAL
        R_METAL_X1_{self.element} {self.left} {self.centre} {r_metal/2}
        R_METAL_X2_{self.element} {self.centre} {self.right} {r_metal/2}
        R_METAL_Y1_{self.element} {self.near} {self.centre} {r_metal/2}
        R_METAL_Y2_{self.element} {self.centre} {self.far} {r_metal/2}
        R_METAL_SEMI_Z_{self.element} {self.centre} {self.bottom} {r_contact}
        """


class Bus:
    """A unit cell representing a SPICE model of metal bus bar segment.

    The Bus bar class is very similar to the Metal class with one important
    exception: the voltage of this node connects to the a voltage source
    for sweeping the bias. Thus these nodes are all at the same voltage.
    """

    def __init__(
            self,
            idx,
            metal_height: float,
            x_length: float,
            y_length: float,
            resistivity_metal: float,
            resistivity_contact: float
        ):
        """
        Parameters
        ----------
        idx : Tuple
            The 3D index of this element in the grid e.g. (1, 2, 3)
        metal_height : float
            The metalisation height of the grid finger in meters for this cell
        x_length  : float
            The metalisation width of the grid finger in meters for this cell. This
            depends on the resolution of the discretiation in the X direction.
        y_length : float
            The metalisation width of the grid finger in meters for this cell. This
            depends on the resolution of the discretiation in the Y direction.
        resistivity_metal : float
            The resistivity of the metal in Ohm meter for this cell.
        resistivity_contact : float
            The resistivity of the metal-semiconductor contact in Ohm meter for
            this cell.
        """
        self.idx = idx
        self.height = metal_height
        self.width = x_length
        self.length = y_length
        self.resistivity_metal = resistivity_metal
        self.resistivity_contact = resistivity_contact
        self.left = f"NX_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.near = f"NY_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.top = "in"
        self.right = f"NX_{self.idx[0]+1}_{self.idx[1]}_{self.idx[2]}"
        self.far = f"NY_{self.idx[0]}_{self.idx[1]+1}_{self.idx[2]}"
        self.bottom = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]+1}"
        self.centre = "in"
        self.element = f"{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"

    def netlist(self):
        """Return the SPICE netlist for this cell.
        """
        r_metal = self.resistivity_metal * self.length / (self.height * self.width)
        r_contact = self.resistivity_contact / (self.width * self.length)

        return f"""
        * BUS
        R_METAL_X1_{self.element} {self.left} {self.centre} {r_metal/2}
        R_METAL_X2_{self.element} {self.centre} {self.right} {r_metal/2}
        R_METAL_Y1_{self.element} {self.near} {self.centre} {r_metal/2}
        R_METAL_Y2_{self.element} {self.centre} {self.far} {r_metal/2}
        R_METAL_SEMI_Z_{self.element} {self.centre} {self.bottom} {r_contact}
        """


class Device:
    """Two-diode model of a solar cell.

    The two diode model:
        - Current source connected from the bottom to the centre node
        - A j01 diode connected from the centre node (anode) to the 
        bottom node (cathode)
        - A j02 diode connected from the centre node (anode) to the 
        bottom node (cathode)
    """

    def __init__(
            self,
            idx,
            x_length: float,
            y_length: float,
            jsc: float,
            sheet_resistance: float,
            relative_illumination_intensity: float = 0.0,
            junction_idx : int = 0
        ):
        """
        Parameters
        ----------
        idx : Tuple
            The 3D index of this element in the grid e.g. (1, 2, 3)
        x_length  : float
            The width of this cell. This depends on the resolution of the 
            discretiation in the X direction.
        y_length : float
            The height of this cell. This depends on the resolution of the 
            discretiation in the Y direction.
        jsc : float
            The short-circuit current density (amps per metre sq.) generated by the solar cell.
        sheet_resistance : float
            The sheet resistance of the emitter in Ohms / square.
        relative_illumination_intensity : float (Default: 0.0)
            A number from 0 to 1 describing the intensity of light this node will generate. This
            number scales the current generated by the current source.
        junction_idx : int (Default: 0)
            Needed when constructing a model of a multijunction solar cell.
        """
        self.idx = idx
        self.x_length = x_length
        self.y_length = y_length
        self.jsc = jsc
        self.sheet_resistance = sheet_resistance
        self.relative_illumination_intensity = relative_illumination_intensity
        self.junction_idx = int(junction_idx)

        self.left = f"NX_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.near = f"NY_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.top = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.right = f"NX_{self.idx[0]+1}_{self.idx[1]}_{self.idx[2]}"
        self.far = f"NY_{self.idx[0]}_{self.idx[1]+1}_{self.idx[2]}"
        self.bottom = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]+1}"
        self.centre = f"N_{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"
        self.element = f"{self.idx[0]}_{self.idx[1]}_{self.idx[2]}"


    def netlist(self):
        # illumination_factor is a number between 0 and 1 to account for the
        # intensity of light on this cell.
        area = self.x_length * self.y_length
        isc = self.jsc * area * self.relative_illumination_intensity
        
        # The resistances are calculated from the resistivity which 
        # can be different in the X and Y direction if the cell's
        # dimensions are different in those directions
        r_sheet_x = self.sheet_resistance * self.y_length / self.x_length
        r_sheet_y = self.sheet_resistance * self.x_length / self.y_length

        return f"""
        * DEVICE
        R_SHEET_X1_{self.element} {self.left} {self.centre} {r_sheet_x/2}
        R_SHEET_X2_{self.element} {self.centre} {self.right} {r_sheet_x/2}
        R_SHEET_Y1_{self.element} {self.near} {self.centre} {r_sheet_y/2}
        R_SHEET_Y2_{self.element} {self.centre} {self.far} {r_sheet_y/2}
        D1_{self.element} {self.centre} {self.bottom} __D1_{self.junction_idx}
        D2_{self.element} {self.centre} {self.bottom} __D2_{self.junction_idx}
        I1_{self.element} {self.bottom} {self.centre} DC {isc}
        """


class Base:
    """Representation of a base layer. This layer is not distributed. It is 
    simply a single resistor, only one of these models is needed per solar 
    cell junction.
    """
    def __init__(
            self,
            idx,
            area: float,
            base_buffer_specific_contact_resistivity: float
    ):
        """
        Parameters
        ----------
        idx : Tuple
            The 3D index of this element in the grid e.g. (1, 2, 3)
        area  : float
            The full surface area of the solar cell.
        base_buffer_specific_contact_resistivity : float
            The specific contact resistivity (Ohm m2) of the base + buffer layers i.e. 
            any layer before the rear contact layer.
        """
        self.idx = idx
        self.area = area
        self.base_buffer_specific_contact_resistivity = base_buffer_specific_contact_resistivity

    def netlist(self):
        # FIXME: need to include element name so this can stack in MJ cells.
        r_base = self.base_buffer_specific_contact_resistivity / self.area  # Approx 1 mOhms
        k = self.idx[2]
        return f"""
        * BASE
        R_BASE_Z N_0_0_{k} N_0_0_{k+1} {r_base}
        """


class RearContact:
    """Representation of a rear contact. This layer is not distributed. It is 
    simply a single resistor, only one of these models is needed per solar 
    cell junction.
    """
    def __init__(
            self,
            idx,
            area: float,
            rear_contact_specific_contact_resistivity: float
    ):
        """
        Parameters
        ----------
        idx : Tuple
            The 3D index of this element in the grid e.g. (1, 2, 3)
        area  : float
            The full surface area of the solar cell.
        rear_contact_specific_contact_resistivity : float
            The resistivity (Ohm m) of the rear contact layer.
        """
        self.idx = idx
        self.area = area
        self.rear_contact_specific_contact_resistivity = rear_contact_specific_contact_resistivity

    def netlist(self):
        # FIXME: need to include element name so this can stack in MJ cells.
        r_rear_contact = self.rear_contact_specific_contact_resistivity / self.area  # Approx 1 mOhms
        k = self.idx[2]
        return f"""
        * REAR CONTACT
        R_REAR_CONTACT_Z N_0_0_{k} 0 {r_rear_contact}
        """


