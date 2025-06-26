import numpy as np
from solcore.spice.grid import HGridPattern
from solcore.spice.netlist import generate_netlist, solve_netlist
from solcore.spice.result import (
    get_characterisic_curve,
    get_electroluminescence,
    get_maximum_power_point,
    get_node_voltages,
    plot_characteristic_curve,
    plot_electroluminescence,
    plot_surface_voltages
)

if __name__ == "__main__":

    # Cell short-circuit current
    jsc = 3000.0

    # Temperature
    temperature = 300.0

    # Grid pattern
    nx, ny = 120, 120
    cell_grid = HGridPattern(10, [4, 4, 4, 4, 4, 4], 3, nx=nx, ny=ny)

    # Homogeneous illumination
    cell_illumination_map = np.ones(nx * ny).reshape((nx, ny))

    # The size of the solar is 3mm x 3mm
    cell_size = (0.003, 0.003) # meters
    
    # Define a list of properies that describe each junction in the solar cell.
    # NB: currently only one junction is working.
    junctions = [
        {
            "jsc": jsc,
            "emitter_sheet_resistance": 100.0,
            "j01": 4e-16,
            "j02": 2e-7,
            "Eg": 1.41,
            "n1": 1.0,
            "n2": 2.0
        }
    ]

    netlist = generate_netlist(
        cell_grid,
        cell_illumination_map,
        cell_size,
        junctions,
        temperature=300,
        show_plots=True
    )

    print("")
    print("This simulation will take a few minutes to run, please wait for results to appear...")
    result = solve_netlist(netlist, temperature, -0.1, 1.5, 0.01)

    V, I = get_characterisic_curve(result)

    plot_characteristic_curve(V, I)

    vmax, pmax, maxidx = get_maximum_power_point(result)

    voltages = get_node_voltages(result)

    plot_surface_voltages(voltages, bias_index=maxidx)

    pv_surface_voltages = voltages[:, :, 1, maxidx]
    el = get_electroluminescence(pv_surface_voltages, is_metal=cell_grid.is_metal)

    plot_electroluminescence(el)