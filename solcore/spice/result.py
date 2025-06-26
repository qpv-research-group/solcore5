"""Functions to extract useful results from the SPICE simulation data.
"""


import numpy as np 
from numpy import ndarray
import matplotlib.pyplot as plt
import itertools
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union
from PySpice.Probe.WaveForm import DcAnalysis


def get_characterisic_curve(result: DcAnalysis) -> tuple[ndarray, ndarray]:
    """Return the solar cells IV curve from the PySpice object.

    Parameters
    ----------
    result : DcAnalysis
        A PySpice analysis object containg solution of the net.

    Returns
    -------
    tuple : (ndarray, ndarray)
        A tuple like (voltage, current).
    """
    voltage = result.sweep.as_ndarray()
    current = result["vin"].as_ndarray()
    return (voltage, current)


def get_maximum_power_point(result: DcAnalysis) -> tuple[float, float, int]:
    """
    Parameters
    ----------
    result : DcAnalysis
        A PySpice analysis object containg solution of the net.
    
    Returns
    -------
    tuple : (float, float, int)
        A tuple like (vmax, pmax, maxidx), where:
            - vmax is the voltage at maximum power (units: V)
            - pmax is the maximum power (units: W)
            - maxidx is the index in the characterisic curve
              corresponding to the maximum power point.
    """
    v, i = get_characterisic_curve(result)
    p = v * i
    idx = np.argmax(p)
    vmax = v[idx]
    return vmax, p[idx], idx


def get_electroluminescence(voltage: ndarray, is_metal=Optional[ndarray], temperature=300.0) -> ndarray:
    """Returns voltage scaled by the Boltzmann approximation so they return values proportional to emission intensity.

    Parameters
    ----------
    voltage : ndarray
        Should be a 2D array containing layer voltages. For example,

            voltages = get_node_voltages(result)
            vmax, pmax, maxidx = get_maximum_power_point(result)
            layer_idx = 1
            voltage = voltages[:, :, layer_index, maxidx]
            get_electroluminescence(voltage)
        
        Here the second layer is selected which is the voltage between the surface
        of the solar cell to ground.
    
    is_metal : ndarray (Optional, Default = None)
        An optional array that is used to mask the prediced electroluminesence.

        Moreover, metalisation blocks the EL, use this array
        to tell the function where the metal is location so that it can 
        reduce the EL to zero in those regions.
    
    temperature : float (Optiona, Default = 300.0)
        The temperature of the solar cell in Kelvins.
    
    Returns
    -------
    el : ndarray
        An array the same shape as the input with values proportional to photon flux.
    """
    q = 1.6e-19
    k = 1.3e-23
    el = np.exp( q * voltage / (k * temperature) )
    if is_metal is not None:
        el[is_metal] = 0.0
    return el


def get_node_voltages(result: DcAnalysis, empty=float("nan")):
    """A 3+1 dimensional array of voltage.

    Parameters
    ----------
    result : DcAnalysis
        A PySpice analysis object containg solution of the net.
    
    empty : float
        The value to use for nodes that do not have a corresponding model. This
        occurs in the Base and Rear Contact layers because layers contain a 
        single object for the whole layer.
    
    Returns
    -------
    voltages : ndarray
        Dimensions are as follows:

            X, Y, Z, N = voltages.shape

        Where X, Y, Z are dimensions of the 3D discretisation and N
        is the number of voltage steps used in the DC sweep.

        For exmaple, say voltage index 10 corresponding the maximum
        power point index, then to get all voltages of the metal layer,

            voltages[:, :, 0, 10]
        
        to get all voltages of the emitter layer,

            voltages[:, :, 1, 10]
    """
    node_names = result.nodes.keys()
    node_names = [x for x in node_names if x.startswith("n")]  # node labels only, not voltage source etc.

    # Get the grid size from the result, we could pass this in, but it easy to determine
    xidxs, yidxs, zidxs = set(), set(), set()
    for key in node_names:
        x, y, z = key.split("_")[1:]
        xidxs.add(int(x))
        yidxs.add(int(y))
        zidxs.add(int(z))

    # The number of nodes in each dimension
    X, Y, Z = max(xidxs), max(yidxs), max(zidxs)

    # Empty array with the size we need, this will contain a voltage
    # at the node coordinate, we can return this so that user can 
    # plot voltage distribution and EL simuations
    bias_voltages = result.sweep.as_ndarray()
    N = bias_voltages.size
    voltages = np.array([empty] * X * Y * Z * N).reshape((X, Y, Z, N))

    # Use the cartesian product to flatten this nested loop - a bit nicer to read
    for i, j, k in itertools.product(range(X), range(Y), range(Z)):
        node_name = f"n_{i}_{j}_{k}"
        if node_name in result.nodes.keys():
            voltages[i, j, k, :] = result[node_name].as_ndarray()
        elif k==0 and "in" in result.nodes.keys():
            voltages[i, j, k, :] = result["in"].as_ndarray()
    
    return voltages


def plot_characteristic_curve(V, I, show=True, path : None | str | Path = None):
    """Plot the characteristic curve

    Parameters
    ----------
    V : ndarray
        Voltage values
    I : ndarray
        Current values
    show : bool (Default: True)
        Immediately render the plot using `plt.show()`.
    path : Optional, str, Path
        Saved plot image to location and file type given, will skip calling `plot.show()`.
    """
    plt.plot(V, I, label="IV curve")
    plt.xlabel("Bias (V)")
    plt.ylabel("Current (A)")
    plt.ylim(ymin=0, ymax=np.max(I)*1.1)
    plt.grid(ls="dotted")

    if show:
        plt.show()
        return

    if path is not None:
        plt.savefig(path)


def plot_surface_voltages(voltages, bias_index, show=True, path=None):
    """Plots the voltage distribution across the surface of the solar cell.

    Parameters
    ----------
    voltages : numpy.ndarray
        4D array of voltages returned from `get_node_voltages`. Dimensions are:
            1. X position index
            2. Y position index
            3. Z layer index
            4. Bias voltage index
    bias_index : int
        The index corresponding to the bias voltage to plot. Usually
        this is the voltage at the maximum power point returned from
        `get_maximum_power_point`.
    show : bool (Default: True)
        Immediately render the plot using `plt.show()`.
    path : Optional, str, Path
        Saved plot image to location and file type given, will skip calling `plot.show()`.
    """

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(11, 5))
    cx1 = ax1.matshow(voltages[:, :, 0, bias_index], cmap="inferno")
    ax1.set_title("Metal Layer Voltages")
    fig.colorbar(cx1)

    cx2 = ax2.matshow(voltages[:, :, 1, bias_index], cmap="inferno")
    ax2.set_title("Emitter Layer Voltages")
    fig.colorbar(cx2)

    if show:
        plt.show()
        return

    if path is not None:
        plt.savefig(path)


def plot_surface_voltages_shared_colormap(voltages, index, show=True, path=None):
    """Plots the voltage distribution across the surface of the solar cell with a single color bar.

    Parameters
    ----------
    voltages : numpy.ndarray
        4D array of voltages returned from `get_node_voltages`. Dimensions are:
            1. X position index
            2. Y position index
            3. Z layer index
            4. Bias voltage index
    bias_index : int
        The index corresponding to the bias voltage to plot. Usually
        this is the voltage at the maximum power point returned from
        `get_maximum_power_point`.
    show : bool (Default: True)
        Immediately render the plot using `plt.show()`.
    path : Optional, str, Path
        Saved plot image to location and file type given, will skip calling `plot.show()`.
    """

    vmin, vmax = np.nanmin(voltages[:, :, :2, index]), np.nanmax(voltages[:, :, :2, index])
    fig, axes = plt.subplots(ncols=2, figsize=(11, 5))
    for layer_idx, ax in enumerate(axes.flat):
        im = ax.matshow(voltages[:, :, layer_idx, index], vmin=vmin, vmax=vmax, cmap="inferno")
        if layer_idx == 0:
            ax.set_title("Metal Layer Voltages")
        else:
            ax.set_title("Emitter Layer Voltages")
    
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    if show:
        plt.show()
        return

    if path is not None:
        plt.savefig(path)


def plot_electroluminescence(el, show=True, path=None):
    """Plots the voltage distribution across the surface of the solar cell.

    Parameters
    ----------
    el : numpy.ndarray
        2D array of EL intensity returned from `get_electroluminescence`
    show : bool (Default: True)
        Display the plot immediately.
    path : bool (Optional)
        If `show=False` path can point to a save location for the plot.
    """

    cx1 = plt.matshow(el, cmap="inferno")
    plt.title("Predicted Electroluminescene")
    plt.colorbar(cx1)
    if show:
        plt.show()
        return
    
    if path is not None:
        plt.savefig(path)
    