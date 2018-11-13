from typing import Optional, Union
from collections import OrderedDict

import pandas as pd
import numpy as np
import pickle as pkl

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib import gridspec
from matplotlib.figure import Figure

from solcore.light_source import LightSource
from solcore.constants import q


class IVPlot(object):
    """ Class that creates a plot of IV data, using a sensible formatting for the axis.
    """

    def __init__(self, data: pd.DataFrame, light: bool = False, spectrum: Optional[LightSource] = None,
                 title: str = '') -> None:
        """ Creates the plot of IV data

        :param data: Data to plot as a Pandas dataframe
        """
        self.data = data
        self.fig, self.axes = plt.subplots(num=title, figsize=(8, 6))
        self.axes.margins(0)

        # This will store the solar cell parameters in the case of light IV
        self.outputs = OrderedDict()

        if light:
            self.plot_light_iv(data, spectrum)
        else:
            self.plot_dark_iv(data)

        self.axes.set_xlabel('Voltage (V)')
        self.axes.set_ylabel('Current (A/m$^2$)')
        self.axes.legend(frameon=False)
        plt.tight_layout()

        # A small pause is needed to actually draw the plot
        plt.draw()
        plt.pause(0.001)

    def plot_dark_iv(self, data: pd.DataFrame) -> None:
        """ Plots the data assuming it is dark IV data, ie. absolute values for the current and log scale.

        :param data: Data to plot as a Pandas dataframe
        :return: None
        """

        for name in data.columns.values:

            # Plotting tunnel junctions IV
            if 'J_TJ' in name:
                label = 'Isolated {}'.format(name.split('_')[1])
                self.axes.semilogy(data['V'], abs(data[name]), label=label)

            # Plotting junctions IV, when isolated and as part of the MJ structure
            elif 'V_MJ' in name:
                id = name.split('V_MJ')[-1]
                p = self.axes.semilogy(data[name], abs(data['J']), 'o--', fillstyle='none', label=name.split('_')[1])
                color = p[-1].get_color()

                new_name = 'J_' + id
                label = 'Isolated {}'.format(id)
                self.axes.semilogy(data['V'], abs(data[new_name]), color=color, label=label)

        # Finally the voltage drop due to series resistance and the total IV.
        self.axes.plot(data['V_Rs'], abs(data['J']), 'o--', fillstyle='none', label='Rs')
        self.axes.semilogy(data['V'], abs(data['J']), 'k', linewidth=4, label='Total IV')

    def plot_light_iv(self, data: pd.DataFrame, spectrum: LightSource) -> None:
        """ Plots the data assuming it is light IV data, ie. moving shifting the region of interest to the first quadrant. It also calculates the solar cell parameters and the cell efficiency if a spectrum is given.

        :param data: Data to plot as a Pandas dataframe
        :param spectrum: A LightSource object with the spectrum, needed to calculate the efficiency of the cell
        :return: None
        """

        # We move the IV in the first quadrant, considering where Isc and Voc are
        sj = 1 if np.interp(0, data['V'], data['J']) >= 0 else -1
        sv = 1 if np.interp(0, data['J'], data['V']) >= 0 else -1

        # We have to track the maximum values for the axes, too.
        max_j = max(-1, np.interp(0, data['V'], sj * data['J']))
        max_v = max(-1, np.interp(0, data['J'], sv * data['V']))

        # In the process, we get the solar cell parameters
        power = sv * data['V'] * sj * data['J']
        idx_max = pd.Series.idxmax(power)
        self.outputs['jsc'] = max_j
        self.outputs['voc'] = max_v
        self.outputs['jmpp'] = sj * data['J'][idx_max]
        self.outputs['vmpp'] = sv * data['V'][idx_max]
        self.outputs['mpp'] = self.outputs['jmpp'] * self.outputs['vmpp']
        self.outputs['ff'] = self.outputs['mpp'] / (self.outputs['jsc'] * self.outputs['voc'])

        if spectrum is not None:
            self.outputs['eta'] = self.outputs['mpp'] / spectrum.power_density

        for name in data.columns.values:

            # Plotting tunnel junctions IV. Actually, we need to review this as TJ have opposite polarity than
            # normal J and therefore should look more like Rs.
            if 'J_TJ' in name:
                sj = 1 if np.interp(0, data['V'], data[name]) >= 0 else -1
                sv = 1 if np.interp(0, data[name], data['V']) >= 0 else -1
                max_j = max(max_j, np.interp(0, data['V'], sj * data[name]))
                label = 'Isolated {}'.format(name.split('_')[1])
                self.axes.plot(sv * data['V'], sj * data[name], label=label)

            # Plotting junctions IV, when isolated and as part of the MJ structure
            elif 'V_MJ' in name:
                id = name.split('V_MJ')[-1]
                sj = 1 if np.interp(0, data[name], data['J']) >= 0 else -1
                sv = 1 if np.interp(0, data['J'], data[name]) >= 0 else -1
                p = self.axes.plot(sv * data[name], sj * data['J'], 'o--', fillstyle='none', label=name.split('_')[1])
                color = p[-1].get_color()

                new_name = 'J_' + id
                sj = 1 if np.interp(0, data['V'], data[new_name]) >= 0 else -1
                sv = 1 if np.interp(0, data[new_name], data['V']) >= 0 else -1
                max_j = max(max_j, np.interp(0, data['V'], sj * data[new_name]))
                label = 'Isolated {}'.format(id)
                self.axes.plot(sv * data['V'], sj * data[new_name], color=color, label=label)

        # Finally the voltage drop due to series resistance and the total IV.
        sj = 1 if np.interp(0, data['V'], data['J']) >= 0 else -1
        sv = 1 if np.interp(0, data['J'], data['V']) >= 0 else -1
        self.axes.plot(sv * data['V_Rs'], sj * data['J'], 'o--', fillstyle='none', label='Rs')
        self.axes.plot(sv * data['V'], sj * data['J'], 'k', linewidth=4, label='Total IV')

        # And adjust the limits of the axis to fit everything inside
        self.axes.set_xlim(right=1.1 * max_v)
        self.axes.set_ylim(bottom=0, top=1.1 * max_j)
        self.axes.axvline(x=0, ymin=0, ymax=1.1 * max_j, color='k', linestyle='--')


class QEPlot(object):
    """ Class that creates a plot of QE data, using a sensible formatting for the axis.
    """

    def __init__(self, data: pd.DataFrame, spectrum: Optional[LightSource] = None, title: str = '') -> None:
        """ Creates the plot of QE data

        :param data: Data to plot as a Pandas dataframe
        :param spectrum: Solar spectrum, only needed to calculate the currents
        """
        self.data = data
        self.fig, self.axes = plt.subplots(num=title, figsize=(8, 6))
        self.axes.margins(0)

        # This will store the Jsc calculated for each junction if a solar spectrum is given
        self.outputs = OrderedDict()

        self.plot_qe(data, spectrum)

        self.axes.set_xlabel('Wavelength (nm)')
        self.axes.set_ylabel('QE (-)')
        self.axes.legend(frameon=False)
        plt.tight_layout()

        # A small pause is needed to actually draw the plot
        plt.draw()
        plt.pause(0.001)

    def plot_qe(self, data: pd.DataFrame, spectrum: LightSource) -> None:
        """ Plots the QE data and the different components, if any. If a spectrum is provided, it also calculates the Jsc of each junction

        :param data: Data to plot as a Pandas dataframe
        :param spectrum: A LightSource object with the spectrum, needed to calculate the Jsc of each junction
        :return: None
        """
        wl = data['WL'] * 1e9  # Wavelength in nanometers
        color = 'red'  # Default color; not really use
        i = 0  # Counter
        markerstyle = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
        colors = ('r', 'g', 'b', 'm', 'c')
        losses = ('srh', 'rad', 'aug', 'surf', 'surb')
        losses_done = False

        # To accumulate the total light absorbed
        total_A = np.zeros_like(data['WL'])

        # And then the reflection
        y1 = 1 - data['R']
        y2 = np.ones_like(data['R'])
        self.axes.fill_between(wl, y1, y2, color='yellow', alpha=0.5, linewidth=3, label='R')

        # Now we plot what happens with light absorbed
        for name in data.columns.values:
            if 'A' in name:
                total_A += data[name]
                p = self.axes.plot(wl, data[name], '--', linewidth=3, label=name)
                color = p[-1].get_color()
                i = 0
                y1 = data[name]
            elif 'EQE' in name and len(name) <= 5:
                self.axes.plot(wl, data[name], linewidth=3, color=color, label=name)
            elif any(loss in name for loss in losses):
                # For the PDD model we have losses. We plot them as shaded areas
                label = losses[i] if (not losses_done and i < 5) else None
                y2 = y1
                y1 = y1 - data[name]
                self.axes.plot(wl, y1, ':', color=color)
                self.axes.fill_between(wl, y1, y2, facecolor=colors[i], alpha=0.5, label=label)
                i += 1
                losses_done = losses_done or i > 4
            elif 'EQE' in name:
                # And for the DA model we have contributions to the QE. We plot them with different markers.
                self.axes.plot(wl, data[name], color=color, marker=markerstyle[i], fillstyle='none', markevery=5,
                               label=name)
                i += 1

            # If a spectrum is given, we calculate the current related to each of the QE
            if spectrum is not None and name[:4] == 'EQE':
                current = data[name].values * spectrum.spectrum(x=wl, output_units='photon_flux_per_nm')[1]
                self.outputs[name] = q * np.trapz(current, wl)

        self.axes.plot(wl, total_A, 'k--', linewidth=2, label='Total A')

        self.axes.set_xlim(left=min(wl), right=max(wl))
        self.axes.set_ylim(bottom=0, top=1)


class BandstructurePlot(object):
    """ Class that creates a plot of bandstructure data and other position-dependent parameters data, using a sensible formatting for the axis.
    """

    def __init__(self, data: pd.DataFrame, title: str = '') -> None:
        """ Launches the creation of the plots, which actually consist on several separated subplots. We create:

        - Energy plot: Ec, Ev, Efe, Efh
        - Carrier density plot: Nc, Nv, Nd, Na, n, p, ni, rho  <-- Maybe a bit too much. Consider removing Nc, Nv
        - G-R plot: G and the recombination mechanisms SRH, Rad and Aug

        :param data: Data to plot as a Pandas dataframe
        """
        self.data = data
        self.fig, (self.energy, self.carriers, self.gr) = plt.subplots(num=title, figsize=(8, 8), nrows=3, ncols=1,
                                                                       sharex='all')

        self.plot_energy(data)
        self.plot_carriers(data)
        self.plot_gr(data)

        plt.tight_layout()

        # A small pause is needed to actually draw the plot
        plt.draw()
        plt.pause(0.001)

    def plot_energy(self, data: pd.DataFrame) -> None:
        """

        :param data:
        :return:
        """
        energy_data = ('Ec', 'Ev', 'Efe', 'Efh')
        colors_and_styles = {'Ec': 'b-', 'Ev': 'r-', 'Efe': 'b--', 'Efh': 'r--'}
        i = 0

        for name in data.columns.values:
            if any(ed in name for ed in energy_data):
                x = data['x_{}'.format(name.split('_')[-1])] * 1e9
                marker = colors_and_styles[name.split('_')[0]]
                label = name.split('_')[0] if i < 4 else ''
                self.energy.plot(x, data[name], marker, label=label)
                i += 1

        self.energy.set_ylabel('Energy (eV)')
        self.energy.tick_params(labelbottom=False)
        self.energy.legend(frameon=False)

    def plot_carriers(self, data: pd.DataFrame) -> None:
        """

        :param data:
        :return:
        """
        carrier_data = ('Nd', 'Na', 'n_', 'p_')
        colors_and_styles = {'Nd': 'b--', 'Na': 'r--', 'n': 'b-', 'p': 'r-'}
        i = 0

        for name in data.columns.values:
            if any(ed in name for ed in carrier_data):
                x = data['x_{}'.format(name.split('_')[-1])] * 1e9
                marker = colors_and_styles[name.split('_')[0]]
                label = name.split('_')[0] if i < 4 else ''
                self.carriers.semilogy(x, data[name] * 1e-6, marker, label=label)
                i += 1

        self.carriers.set_ylabel('Carrier density,\nDoping (cm$^{-3}$)')
        self.carriers.tick_params(labelbottom=False)
        self.carriers.legend(frameon=False)

    def plot_gr(self, data: pd.DataFrame) -> None:
        """

        :param data:
        :return:
        """
        self.gr.margins(0)
        gr_data = ('G', 'Rsrh', 'Rrad', 'Raug')
        colors_and_styles = {'G': 'r-', 'Rsrh': 'b-', 'Rrad': 'c-', 'Raug': 'g-'}
        i = 0

        for name in data.columns.values:
            if any(ed in name for ed in gr_data) and 'GR' not in name:
                x = data['x_{}'.format(name.split('_')[-1])] * 1e9
                marker = colors_and_styles[name.split('_')[0]]
                label = name.split('_')[0] if i < 4 else ''
                self.gr.semilogy(x, data[name] * 1e-6, marker, label=label)
                i += 1

        self.gr.set_ylabel('Generation,\nRecombination\n(cm$^{-3}$s$^{-1}$)')
        self.gr.set_xlabel('Position (nm)')
        self.gr.legend(frameon=False)


def save_plot(plot: Union[IVPlot, QEPlot, BandstructurePlot], filename: str) -> None:
    """ Saves a plot as a pickle object

    :param plot:
    :param filename:
    :return:
    """
    with open(filename, 'wb') as f:  # should be 'wb' rather than 'w'
        pkl.dump(plot, f)


def load_plot(filename: str, plot: bool = False) -> Union[IVPlot, QEPlot, BandstructurePlot]:
    """ Loads a plot from a pickle object

    :param filename: The name of the file to load
    :param plot: whether the figure should be plotted or not. Default=False.
    :return: An object containing the plot, of class IVPlot, QEPlot or BandstructurePlot
    """
    with open(filename, 'rb') as f:  # should be 'wb' rather than 'w'
        out = pkl.load(f)

    # A small pause is needed to actually draw the plot
    if plot:
        out.fig.canvas.draw_idle()
        plt.pause(0.001)

    return out


if __name__ == '__main__':
    import os

    # # df = pd.read_csv('DA_2J_IV_light.csv')
    # df = pd.read_csv('DA_2J_IV_dark.csv')
    #
    # # sp = LightSource(source_type='standard', version='AM1.5g', x=df['WL'])
    #
    df = pd.read_csv(os.path.join(os.path.expanduser('~'), 'Desktop', 'DA_2J_PDD_IV_light.csv'))
    my_fig = IVPlot(df, light=True, title='DA_2J_PDD_IV_light')
    plt.show()

    # df = pd.read_csv(os.path.join(os.path.expanduser('~'), 'Desktop','DA_2J_PDD_QE.csv'))
    #
    # sp = LightSource(source_type='standard', version='AM1.5g', x=df['WL'])
    #
    # my_fig = QEPlot(df, spectrum=sp, title='DA_2J_PDD_QE')
    # plt.show()

    # df = pd.read_csv(os.path.join(os.path.expanduser('~'), 'Desktop', 'DA_2J_PDD_equilibrium.csv'))
    #
    # my_fig = BandstructurePlot(df, title='DA_2J_PDD_equilibrium')
    # plt.show()

    # filename = 'DA_2J_PDD_short_circuit.csv'
    # df = pd.read_csv(os.path.join(os.path.expanduser('~'), 'Desktop', filename))
    #
    # my_fig = BandstructurePlot(df, title=filename)
    # plt.show()
