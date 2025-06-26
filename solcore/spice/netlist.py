import numpy as np 
import matplotlib.pyplot as plt
import itertools
from PySpice.Spice.Netlist import Circuit
from solcore.spice.grid import GridPattern
from solcore.spice.model import (
    Header,
    Diodes,
    Metal,
    Bus,
    Device,
    Base,
    RearContact
)
from typing import TYPE_CHECKING, Optional
if TYPE_CHECKING:
    pass


def generate_netlist(
    cell_metalisation_pattern: np.ndarray | GridPattern,
    cell_illumination_map: np.ndarray,
    cell_size: tuple[float, float],
    junctions: list[dict],
    metal_height: float = 3e-6,
    metal_resistivity: float = 3.5e-8,
    metal_semiconductor_specific_contact_resistivity: float = 6.34e-10,
    base_buffer_specific_contact_resistivity: float = 1.2e-8,
    rear_contact_specific_contact_resistivity: float = 3.5e-6,
    temperature: float = 300.0,
    show_plots=False
) -> str:
    """
    Returns a string that is a SPICE netlist.

    Parameters
    ----------
    cell_metalisation_pattern : np.ndarray | GridPattern
        A 2D array showing how metalisation is applied to the solar cell, a GridPattern object, 
        or a Path to an image specified as a string or a pathlib.Path

        The 2D image will be interpretted as followed where "px" is the gray scale value of the pixel:
            - Bus bar: px > 0.8
            - Grid finger: 0.2 < px < 0.8
            - No Metal: px < 0.2

        If the 2D array is read from an image file it will be normalised and values scaled between 0 and 1.

    cell_illumination_map: np.ndarray
        A 2D array providing the illumination distribution over the solar cell's surface.
    
    cell_size : Tuple[float, float]
        The tuple gives the edge length in the x and y direction of solar cell. This is used to 
        calculate the length and width of each pixel in the input image. Units: m
    
    junction : List[Dict]
        A list of dictionaries containing solar cell junction parameters. All keys are required, for example,

            junctions = [
                {
                    "jsc": 30000, # A / m2
                    "emitter_sheet_resistance": 100.0, # Ohm / sq
                    "j01": 4e-16, # A / m2
                    "j02": 2e-7, # A / m2
                    "Eg": 1.41, # eV
                    "n1": 1.0, # dimensionless
                    "n2": 2.0  # dimensionless
                }
            ]
        
        This has been implemented as a list to support multi-junction devices, however, that is not yet
        working.

    metal_height : float
        The height of the bus bar and the grid fingers. Units: m
    
    metal_resistivity : float
        The resisitivity of the metal used to form the bus bar and grid fingers. Units: Ohm m
    
    metal_semiconductor_specific_contact_resistivity : float
        The specific contact resisitivty of the metal-semiconductor layer. Units: Ohm m2
    
    base_buffer_specific_contact_resistivity : float
        The specific contact resisitivty of the base and buffer layers. Units: Ohm m2 
    
    rear_contact_specific_contact_resistivity : float
        The specific contact resisitivty of the rear contact layer. Units: Ohm m2 
    
    temperature : float
        The temperature of the solar cell. Unit: K
    
    show_plots : bool
        Plot some of the input data, this is useful for debugging the metalisation and illumination images.
    
    Discussion
    ----------

    The net list is intended to be run using PySpice Circuit objects; for this reason, the `.DC` command 
    is not included, nor is the `.end` statement. When PySpice runs, it will append these to the net
    list. To run the net list in an external simulator, simply append the two lines.

        .DC vin -0.1 1.4 0.1
        .end

    The first line is the `.DC` command that tells SPICE to sweep the voltage; in this example, the
    voltage starts at -0.1 and ends at 1.4 in steps of 0.1 volts. The second line is a command that 
    tells SPICE it has reached the end of the netlist.
    """

    # Check we have all the information we need to continue
    must_haves = [
        "jsc",
        "emitter_sheet_resistance",
        "j01",
        "j02",
        "Eg",
        "n1",
        "n2"
    ]
    try:
        for junction in junctions:
            for key in must_haves:
                junction[key]
    except KeyError:
        raise KeyError(f"Do all your junction dictionaries have the required keys? The key '{key}' is missing.")

    if isinstance(cell_metalisation_pattern, GridPattern):
        cell_metalisation_pattern = cell_metalisation_pattern.as_array()
    
    if not (cell_metalisation_pattern.shape == cell_illumination_map.shape):
        raise ValueError(f"The metalisation mask {cell_metalisation_pattern.shape} and cell illumination map {cell_illumination_map.shape} need to have the same size.")
    
    # The first layer of the solar cell contain SPICE model for the metalisation
    layers, generation_map = _create_grid_layer(
        cell_metalisation_pattern=cell_metalisation_pattern,
        cell_illumination_map=cell_illumination_map,
        cell_size=cell_size,
        metal_height=metal_height,
        metal_resistivity=metal_resistivity,
        metal_semiconductor_specific_contact_resistivity=metal_semiconductor_specific_contact_resistivity,
        show_plots=show_plots
    )
    X, Y = layers.shape[:2]
    unit_x = cell_size[0] / X
    unit_y = cell_size[1] / Y
    
    # For each junction create a layer containing a SPICE model of the PV cell.
    # This model defines diodes we also need to save the associated
    # header information.
    headers = [Header(temperature=temperature)]
    for junction in junctions:
        # Twiddle the dict around a bit so we can call the function below easily
        junction.update({"generation_map": generation_map})
        jsc = junction.pop("jsc")
        layers, header_info = _update_layers_with_junction_model(layers, cell_size, jsc, **junction)
        headers.append(header_info)

        # TODO: is this part sensible? Maybe we should treat the Base model as a sheet?
        # Add a buffer layer SPICE model. Note that this is not a distributed model,
        # it is just a single resistor that connects from N_0_0_z to N_0_0_z+1, where z 
        # is the index of base layer.
        layers = np.dstack((layers, np.array([None] * X * Y).reshape((X, Y, 1))))
        idx = (0, 0, layers.shape[2] - 1)
        layers[idx] = Base(
            idx,
            area=cell_size[0] * cell_size[1],
            base_buffer_specific_contact_resistivity=base_buffer_specific_contact_resistivity
        )

    # Add a rear contact layer SPICE model. Note that this is not a distributed model,
    # it is just a single resistor that connects to the N_0_0_z node, where z 
    # is the index of the rear contact layer.
    layers = np.dstack((layers, np.array([None] * X * Y).reshape((X, Y, 1))))
    idx = (0, 0, layers.shape[2] - 1)
    layers[0, 0, -1] = RearContact(idx, cell_size[0] * cell_size[1], rear_contact_specific_contact_resistivity=rear_contact_specific_contact_resistivity)

    # Convert the model objects to a net list string
    netlist = [info.netlist() for info in headers]
    X, Y, Z = layers.shape
    for (k, i, j) in itertools.product(range(Z), range(X), range(Y)):
        idx = (i, j, k)
        if layers[idx]:
            netlist.append(layers[idx].netlist())

    # netlist is list of strings, let's convert that into a single string so we can 
    # throw it into SPICE
    netlist = "".join(netlist)
    return netlist


def solve_netlist(net: str, temperature: float, Vstart: float, Vstop: float, Vstep : float):
    celsius = temperature - 273.15
    from solcore.spice.spice import spice as SpiceConfig

    with open("grid.net", "w") as f:
        f.write(net)

    cir = Circuit(net)
    simulator = cir.simulator(temperature=celsius, nominal_temperature=celsius, spice_command=SpiceConfig.engine)
    return simulator.dc(vin=slice(Vstart, Vstop, Vstep))


#
# Helper function private to this module
#


def _create_grid_layer(
    cell_metalisation_pattern: np.ndarray | GridPattern, # the cells that are grid fingers
    cell_illumination_map: np.ndarray,
    cell_size: tuple[float, float], # cell_size = (x_distance, y_distance), the edge lengths of the solar cell, assumes rectangular shape.
    metal_height: float = 3e-6, # Height: m, of grid fingers
    metal_resistivity: float = 3.5e-6, # Resistivity: Ohm m, of the metal used for front contacts
    metal_semiconductor_specific_contact_resistivity: float = 6.34e-6, # Specific contact resistivity: Ohm m2, of metal-semiconductor layer
    show_plots=False
):
    """Internal function that returns SPICE objects from cell and grid information.

    Parameters
    ----------
    cell_metalisation_pattern : np.ndarray | GridPattern
        A 2D array or a GridPattern object showing how metalisation is applied to the solar cell.
        
        The 2D image will be interpretted as followed where "px" is the gray scale value of the pixel:
            - Bus bar: px > 0.8
            - Grid finger: 0.2 < px < 0.8
            - No Metal: px < 0.2

    cell_illumination_map: np.ndarray
        A 2D array showing the illumination distribution over the solar cell's surface.

    cell_size : Tuple[float, float]
        The tuple gives the edge length in the x and y direction of solar cell. This is used to 
        calculate the length and width of each pixel in the input image. Units: m
    
    metal_height : float
        The height of the bus bar and the grid fingers. Units: m
    
    metal_resistivity : float
        The resisitivity of the metal used to form the bus bar and grid fingers. Units: Ohm m
    
    metal_semiconductor_specific_contact_resistivity : float
        The specific contact resisitivty of the metal-semiconductor layer. Units: Ohm m2

    show_plots : bool
        Plot some of the input data, this is useful for debugging the metalisation and illumination images.
    """
    if isinstance(cell_metalisation_pattern, GridPattern):
        cell_metalisation_pattern = cell_metalisation_pattern.as_array()
    
    if not (cell_metalisation_pattern.shape == cell_illumination_map.shape):
        raise ValueError(f"The metalisation mask {cell_metalisation_pattern.shape} and cell illumination map {cell_illumination_map.shape} need to have the same size.")
    
    # The image of the metalisation can be broken down into three distinct parts:
    # 1. The bus
    # 2. The grid fingers
    # 3. The device (i.e. the area without any metal)
    # The metalisation map is a gray scale image, any region that is above 90%
    # white is mapped to the bus, any region less that 10% white is mapped 
    # as the device and region between 40% and 60% is mapped to the grid fingers.
    norm_metalisation_map = cell_metalisation_pattern / cell_metalisation_pattern.max()
    is_bus = np.where(norm_metalisation_map > 0.8, 1, 0)
    is_finger = np.where((norm_metalisation_map > 0.20)&(norm_metalisation_map < 0.80), 1, 0)
    is_not_metal = np.where(norm_metalisation_map < 0.2, 1, 0)

    if show_plots:

        # These plots show how the image has been processed. If the metalisation
        # is not behaving then this should show you where things are going wrong.
        fig, ax = plt.subplots(2, 2)

        ax[0, 0].matshow(norm_metalisation_map, cmap="gray")
        ax[0, 0].set_title("Grid Image")

        ax[0, 1].matshow(is_bus, cmap="gray")
        ax[0, 1].set_title("Bus bar only")

        ax[1, 0].matshow(is_finger, cmap="gray")
        ax[1, 0].set_title("Grid fingers only")

        ax[1, 1].matshow(is_not_metal, cmap="gray")
        ax[1, 1].set_title("Solar cell only")
        
        fig.tight_layout()
        plt.show()
    
    # Create a 3D matrix to hold each cell:
    #  - Cells with indices [:,:,0] are Metal or Bus objects
    #  - Cells with indices [:,:,1] are semiconductor device cells
    X, Y = cell_illumination_map.shape
    xscale, yscale = cell_size
    unit_width = xscale / X  # the width of each spice cell
    unit_length = yscale / Y # the length of each spice cell
    grid = np.array([None] * X * Y).reshape(X, Y, 1)

    for idx in itertools.product(range(X), range(Y), [0]):
        i, j, k = idx
        # Layer #0
        if is_bus[i, j]:
            grid[i, j, k] = Bus(idx, metal_height, unit_width, unit_length, metal_resistivity, metal_semiconductor_specific_contact_resistivity)
        elif is_finger[i, j]:
            grid[i, j, k] = Metal(idx, metal_height, unit_width, unit_length, metal_resistivity, metal_semiconductor_specific_contact_resistivity)

    # Make a "generation map", this does not correspond to the illumiation map 
    # because of shading by the metalisation cells. Here we assume that any 
    # metalisation reduces the generation in that cell to zero.
    has_generation = np.where(~np.where((is_bus | is_finger), True, False), 1, 0)
    generation_map = has_generation * cell_illumination_map

    if show_plots:
        fig, (ax1, ax2) = plt.subplots(ncols=2)
        ax1.matshow(cell_illumination_map, vmin=0.0, vmax=1.0)
        ax1.set_title("Illumination map")

        p2 = ax2.matshow(generation_map, vmin=0.0, vmax=1.0)
        ax2.set_title("Generation map")

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(p2, cax=cbar_ax)
        plt.show()

    return grid, generation_map


def _update_layers_with_junction_model(
        layers: np.ndarray,
        cell_size: tuple[float, float],
        jsc: float, # Jsc: A / m2, short-circuit current generated by the solar cell (convention: positive)
        emitter_sheet_resistance: float = 100, # The emitter sheet resistance
        j01: float = 4.0e-20,
        j02: float = 2e-11,
        Eg: float = 1.41,
        n1: float = 1.0,
        n2: float = 2.0,
        generation_map: Optional[np.ndarray] = None
):
    """Internal function that appends solar cell SPICE objects to the layer structure.
    """
    X, Y, Z = layers.shape
    xscale, yscale = cell_size
    unit_width = xscale / X  # the width of each spice cell
    unit_length = yscale / Y # the length of each spice cell
    area = unit_length * unit_width

    def get_relative_illumination_intensity(i, j):
        if generation_map is not None:
            return generation_map[i, j]
        return 0.0

    # 0 if this is the first device layer,
    # 1 is this is the second device layer etc.
    junction_idx = len([isinstance(x, Device) for x in layers[0,0,:] if isinstance(x, Device)])
    layers = np.dstack((layers, np.array([None] * X * Y).reshape((X, Y, 1))))
    X, Y, Z = layers.shape
    k = Z - 1
    # NOTE: the device layer at the edge should be a different type, but it not yet implemented.

    for idx in itertools.product(range(X), range(Y), [k]):
        layers[idx] = Device(
            idx,
            unit_width,
            unit_length,
            jsc,
            emitter_sheet_resistance,
            relative_illumination_intensity=get_relative_illumination_intensity(*idx[:2]),
            junction_idx=junction_idx
        ) 
        layers[idx].bottom = f"N_0_0_{idx[2]+1}"

    # Create diode model that need to appear in the SPICE header
    header = Diodes(Eg, j01, j02, area, n1=n1, n2=n2, junction_idx=junction_idx)
    return layers, header

