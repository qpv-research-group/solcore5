import tkinter as tk
from tkinter import ttk

from matplotlib import axes, figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

import numpy as np

from solcore.light_source import LightSource
from solcore import eVnm, nmHz, nmJ, mJ
from solcore.constants import q, h, c

from typing import Callable, Dict, Union


class SpectrumTab(ttk.Frame):
    """ This tab is used to create illumination spectra
    """

    def __init__(self, parent: ttk.Notebook) -> None:
        """ Constructor of the class

        :param parent: The notebook that serves as parent of this tab
        """
        super(SpectrumTab, self).__init__(parent)

        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        # We create the variables containing the details of the spectrum
        self.source_type = tk.StringVar(value='standard')
        self.concentration = tk.DoubleVar(value=1.0)
        self.units = tk.StringVar(value='power_density_per_nm')
        self.units_x = tk.StringVar(value='nm')
        self.x_min = tk.DoubleVar(value=300.0)
        self.x_max = tk.DoubleVar(value=1100.0)
        self.x_points = tk.IntVar(value=200)
        self.old_units = 'power_density_per_nm'
        self.old_units_x = 'nm'

        # We add some callbacks to update the plot when the range entries are changed
        self.x_min.trace_add("write", self.update_plot)
        self.x_max.trace_add("write", self.update_plot)
        self.x_points.trace_add("write", self.update_plot)

        # Next we create the variable that will hold the currently calculated LightSource object
        self.light_source = None

        # We add two frames, one for the controls on the left and another one for the plot on the right
        self.fig, self.axes = self.create_plot_area()
        self.control_frame = self.create_control_area()

        # We initialize the tab with the Standard AM1.5g solar spectrum
        self.spectrum = StandardSpectrum(parent=self.control_frame, row=4, update=self.update_plot)
        self.update_plot()

    def create_plot_area(self) -> (figure.Figure, axes.Axes):
        """ Creates the plot area in which to show the spectra

        :return: The figure and axes
        """
        plot_frame = ttk.Frame(self)
        plot_frame.grid(column=1, row=0, sticky=tk.NSEW)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)

        fig = figure.Figure(tight_layout=True)
        ax = fig.add_subplot(111)

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(column=0, row=0, sticky=tk.NSEW)

        toolbar_frame = ttk.Frame(master=plot_frame)
        toolbar_frame.grid(column=0, row=1, sticky=tk.NSEW)
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()

        return fig, ax

    def create_control_area(self):
        """

        :return:
        """
        control_frame = ttk.Frame(self)
        control_frame.grid(column=0, row=0, sticky=tk.NSEW)

        general_properties_frame = ttk.Frame(control_frame)
        general_properties_frame.grid(column=0, row=0, sticky=tk.NSEW)

        # Source type
        source_type_label = ttk.Label(general_properties_frame, text="Source type")
        spectrum_type_box = ttk.Combobox(general_properties_frame, state='readonly',
                                         values=LightSource.type_of_source, textvariable=self.source_type)
        source_type_label.grid(column=0, row=0, sticky=tk.NSEW)
        spectrum_type_box.grid(column=1, row=0, columnspan=2, sticky=tk.NSEW)
        spectrum_type_box.bind('<<ComboboxSelected>>', self.update_plot)

        # Units
        units_label = ttk.Label(general_properties_frame, text="Units")
        units_box = ttk.Combobox(general_properties_frame, state='readonly', values=LightSource.output_units,
                                 textvariable=self.units)
        units_label.grid(column=0, row=1, sticky=tk.NSEW)
        units_box.grid(column=1, row=1, columnspan=2, sticky=tk.NSEW)
        units_box.bind('<<ComboboxSelected>>', self.update_plot)

        # Range
        range_label = ttk.Label(general_properties_frame, text="Range (min, max, points) in: ")
        x_units_label = ttk.Label(general_properties_frame, textvar=self.units_x)
        x_max_label = ttk.Entry(general_properties_frame, textvar=self.x_min, width=10)
        x_min_label = ttk.Entry(general_properties_frame, textvar=self.x_max, width=10)
        x_points_label = ttk.Entry(general_properties_frame, textvar=self.x_points, width=10)
        range_label.grid(column=0, row=2, columnspan=2, sticky=tk.NSEW)
        x_units_label.grid(column=2, row=2, sticky=tk.NSEW)
        x_max_label.grid(column=0, row=3, sticky=tk.NSEW)
        x_min_label.grid(column=1, row=3, sticky=tk.NSEW)
        x_points_label.grid(column=2, row=3, sticky=tk.NSEW)

        return control_frame

    def update_plot(self, *args) -> None:
        """ Updates the plot whenever one of the settings changes and also update the LightSource variable

        :return: None
        """
        self.update_units_labels_and_values()
        self.light_source = self.calculate_light_source()
        x_label, y_label = self.get_axes_labels()

        self.axes.clear()
        self.axes.plot(*self.light_source.spectrum(), 'r')
        self.axes.set_xlabel(x_label)
        self.axes.set_ylabel(y_label)

        # recompute the ax.dataLim and update ax.viewLim using the new dataLim
        self.axes.relim()
        self.axes.autoscale_view()

        # re-draw the canvas
        self.fig.canvas.draw_idle()

    def update_units_labels_and_values(self) -> None:
        """ Updates the x units labels and fields

        :return: None
        """
        # If x units haven't changed, we do nothing
        new_x = self.units.get().split('_')[-1]
        old_x = self.units_x.get()
        if new_x == old_x:
            return

        self.units_x.set(new_x)

        old_min = self.x_min.get()
        old_max = self.x_max.get()

        if all(t in ('nm', 'eV') for t in (old_x, new_x)):
            new_min = eVnm(old_min)
            new_max = eVnm(old_max)
        elif all(t in ('nm', 'J') for t in (old_x, new_x)):
            new_min = nmJ(old_min)
            new_max = nmJ(old_max)
        elif all(t in ('nm', 'm') for t in (old_x, new_x)):
            factor = 1e-9 if old_x == 'nm' else 1e9
            new_min = factor * old_min
            new_max = factor * old_max
        elif all(t in ('nm', 'hz') for t in (old_x, new_x)):
            new_min = nmHz(old_min)
            new_max = nmHz(old_max)
        elif all(t in ('m', 'J') for t in (old_x, new_x)):
            new_min = mJ(old_min)
            new_max = mJ(old_max)
        elif all(t in ('m', 'eV') for t in (old_x, new_x)):
            factor = h * c / q
            new_min = factor / old_min
            new_max = factor / old_max
        elif all(t in ('m', 'hz') for t in (old_x, new_x)):
            factor = c
            new_min = factor / old_min
            new_max = factor / old_max
        elif all(t in ('J', 'eV') for t in (old_x, new_x)):
            factor = q if old_x == 'eV' else 1 / q
            new_min = factor * old_min
            new_max = factor * old_max
        elif all(t in ('J', 'hz') for t in (old_x, new_x)):
            factor = 1 / h if old_x == 'J' else h
            new_min = factor * old_min
            new_max = factor * old_max
        else:
            # eV <-> hz
            factor = q / h if old_x == 'eV' else h / q
            new_min = factor * old_min
            new_max = factor * old_max

        # Now we have to check if maximum and minimum are in the correct order, reversing them, otherwise
        if new_min > new_max:
            new_min, new_max = new_max, new_min

        self.x_min.set(format(new_min, '.4'))
        self.x_max.set(format(new_max, '.4'))

    def calculate_light_source(self) -> LightSource:
        """ Calculate the lightsource based on the chosen options

        :return: None
        """
        # First we create the x values
        try:
            x_min = self.x_min.get()
            x_max = self.x_max.get()
            x_points = self.x_points.get()
        except Exception:
            # If any of the parameters causes problems, we return the current light source
            return self.light_source

        # Now we retrieve the options
        options = self.spectrum.get_options()

        # And now we are ready to create the spectrum
        light_source = LightSource(source_type=self.source_type.get(),
                                   x=np.linspace(x_min, x_max, x_points),
                                   output_units=self.units.get(),
                                   concentration=self.concentration.get(),
                                   **options)

        return light_source

    def get_axes_labels(self) -> (str, str):
        """ Gets the correct labels for the axis

        :return: A tuple with the x and y labels
        """
        units = self.units.get()

        # First, the x axes
        x_units = units.split('_')[-1]
        if x_units in ('nm', 'm'):
            x_label = 'Wavelength ({})'.format(x_units)
        elif x_units == 'hz':
            x_label = 'Frequency (hz)'
        else:
            x_label = 'Energy ({})'.format(x_units)

        # And now the y axes
        if units.split('_')[0] == 'power':
            y_label = 'Power density (W m$^{{-2}}$ {}$^{{-1}}$)'.format(x_units)
        else:
            y_label = 'Photon flux (photons m$^{{-2}}$ {}$^{{-1}}$)'.format(x_units)

        return x_label, y_label


class StandardSpectrum(ttk.LabelFrame):
    """

    """

    def __init__(self, parent: ttk.Frame, row: int, update: Callable[[], None], col: int = 0) -> None:
        """

        :param parent:
        """
        super(StandardSpectrum, self).__init__(parent, text='Standard spectrum')
        self.grid(column=col, row=row, sticky=tk.NSEW)

        self.version = tk.StringVar(value='AM1.5g')

        source_version_label = ttk.Label(self, text="Version")
        source_version_box = ttk.Combobox(self, state='readonly', values=['AM0', 'AM1.5g', 'AM1.5d'],
                                          textvariable=self.version)
        source_version_label.grid(column=0, row=0, sticky=tk.NSEW)
        source_version_box.grid(column=1, row=0, columnspan=2, sticky=tk.NSEW)
        source_version_box.bind('<<ComboboxSelected>>', update)

    def get_options(self) -> Dict:
        """ Returns a dictionary with the options of this type of spectrum

        :return: A ditcionary with the options
        """
        out = {'version': self.version.get()}
        return out
