import tkinter as tk
from tkinter import ttk

from solcore import material
from solcore.material_data import calculate_mobility
from solcore.constants import q

from solcore.solar_cell_solver import rcwa_options, pdd_options, asc_options

# The default material is a GaAs at 298 K and provides a set of minimum properties that are used if they are not
# available for the material of interest and not alternative is given
# Default values MUST BE in SI units
DefaultMaterial = material("GaAs")(T=298)
default_layer_properties = {'band_gap': DefaultMaterial.band_gap,  # J
                            'electron_affinity': DefaultMaterial.electron_affinity,  # J
                            'eff_mass_electron': DefaultMaterial.eff_mass_electron,  # relative to m0
                            'eff_mass_hh_z': DefaultMaterial.eff_mass_hh_z,  # relative to m0
                            'eff_mass_lh_z': DefaultMaterial.eff_mass_lh_z,  # relative to m0
                            'electron_mobility': calculate_mobility("GaAs", 0, 1),  # m2 V-1 s-1
                            'hole_mobility': calculate_mobility("GaAs", 1, 1),  # m2 V-1 s-1
                            'ni': DefaultMaterial.ni,  # m-3
                            'Nc': DefaultMaterial.Nc,  # m-3
                            'Nv': DefaultMaterial.Nv,  # m-3
                            'electron_minority_lifetime': 3e-6,  # s
                            'hole_minority_lifetime': 2.5e-7,  # s
                            'relative_permittivity': DefaultMaterial.relative_permittivity,  # relative to epsilon0
                            'electron_auger_recombination': 1e-42,  # m6 s-1,
                            'hole_auger_recombination': 1e-42,  # m6 s-1
                            'radiative_recombination': DefaultMaterial.radiative_recombination,  # m3 s-1
                            'Nd': 1,  # m-3
                            'Na': 1,  # m3
                            'electron_diffusion_length': 10e-6,  # m
                            'hole_diffusion_length': 40e-6}  # m

default_junction_properties = {'Sn': 1e10,
                               'Sp': 1e10,
                               'R_shunt': 1e14,
                               'R_series': 0,
                               'A': 1,
                               'Eg': 1.42,
                               'n': 3.5,
                               'jsc': 0,
                               'j01': 1e-15,
                               'j02': 1e-8,
                               'n1': 1,
                               'n2': 2,
                               'pn': 1,
                               'R': 1e-16,
                               'j_peak': 1,
                               'v_peak': 0.1,
                               'j_valley': 0.01,
                               'v_valley': 0.5,
                               'prefactor': 10}

# Units in which the parameters are shown in the GUI
units = {'Sn': 'cm2 s-1',
         'Sp': 'cm2 s-1',
         'R_shunt': 'Ohm',
         'R_series': 'Ohm',
         'A': '-',
         'Eg': 'eV',
         'n': '-',
         'jsc': 'mA cm-2',
         'j01': 'mA cm-2',
         'j02': 'mA cm-2',
         'n1': '-',
         'n2': '-',
         'band_gap': 'eV',
         'electron_affinity': 'eV',
         'relative_permittivity': '-',
         'Na': 'cm-3',
         'Nd': 'cm-3',
         'electron_mobility': 'cm2 V-1 s-1',
         'hole_mobility': 'cm2 V-1 s-1',
         'eff_mass_hh_z': '-',
         'eff_mass_lh_z': '-',
         'eff_mass_electron': '-',
         'electron_minority_lifetime': 'ns',
         'hole_minority_lifetime': 'ns',
         'radiative_recombination': 'cm3 s-1',
         'electron_auger_recombination': 'cm6 s-1',
         'hole_auger_recombination': 'cm6 s-1',
         'electron_diffusion_length': 'nm',
         'hole_diffusion_length': 'nm',
         'pn': '-',
         'R': 'Ohm',
         'j_peak': 'mA cm-2',
         'v_peak': 'V',
         'j_valley': 'mA cm-2',
         'v_valley': 'V',
         'prefactor': 'V-1'}

# Conversion factor such that parameter_GUI(units) = parameter(SI) / conversion factor
conversion = {'Sn': 1e-4,
              'Sp': 1e-4,
              'R_shunt': 1,
              'R_series': 1,
              'A': 1,
              'Eg': 1,
              'n': 1,
              'jsc': 10,
              'j01': 10,
              'j02': 10,
              'n1': 1,
              'n2': 1,
              'band_gap': q,
              'electron_affinity': q,
              'relative_permittivity': 1,
              'Na': 1e6,
              'Nd': 1e6,
              'electron_mobility': 1e-4,
              'hole_mobility': 1e-4,
              'eff_mass_hh_z': 1,
              'eff_mass_lh_z': 1,
              'eff_mass_electron': 1,
              'electron_minority_lifetime': 1e-9,
              'hole_minority_lifetime': 1e-9,
              'radiative_recombination': 1e-6,
              'electron_auger_recombination': 1e-12,
              'hole_auger_recombination': 1e-12,
              'electron_diffusion_length': 1e-9,
              'hole_diffusion_length': 1e-9,
              'pn': 1,
              'R': 1,
              'j_peak': 0.1,
              'v_peak': 1,
              'j_valley': 0.1,
              'v_valley': 1,
              'prefactor': 1}

# Set of parameters required for each layer or junction
properties_junctions = {'Junction-PDD': ['Sn', 'Sp', 'R_shunt', 'R_series'],
                        'Junction-DA': ['Sn', 'Sp', 'R_shunt', 'R_series'],
                        'Junction-2D': ['j01', 'j02', 'n1', 'n2', 'jsc', 'R_shunt', 'R_series'],
                        'Junction-DB': ['A', 'Eg', 'n', 'R_shunt', 'R_series'],
                        'TJ-resistive': ['R', 'pn'],
                        'TJ-parametric': ['j_peak', 'v_peak', 'j_valley', 'v_valley', 'prefactor', 'j01', 'pn']}

properties_layers = {'Junction-PDD': ['band_gap', 'electron_affinity', 'relative_permittivity', 'Na', 'Nd',
                                      'electron_mobility', 'hole_mobility',
                                      'eff_mass_hh_z', 'eff_mass_lh_z', 'eff_mass_electron',
                                      'electron_minority_lifetime', 'hole_minority_lifetime',
                                      'radiative_recombination', 'electron_auger_recombination',
                                      'hole_auger_recombination'],
                     'Junction-DA': ['band_gap', 'relative_permittivity', 'Na', 'Nd',
                                     'electron_mobility', 'hole_mobility', 'eff_mass_hh_z', 'eff_mass_lh_z',
                                     'eff_mass_electron', 'electron_diffusion_length', 'hole_diffusion_length'],
                     'Junction-2D': [],
                     'Junction-DB': [],
                     'TJ-resistive': [],
                     'TJ-parametric': []}

global_settings = [('T', 298, 'K'),
                   ('T_ambient', 298, 'K')]

electrical_settings = [('voltages', '0, 1.2, 100', 'V, V, points'),
                       ('internal_voltages', '-6, 4, 1000', 'V, V, points'),
                       ('mpp', False, ''),
                       ('light_iv', False, ''),
                       ('radiative_coupling', False, ''),
                       ('position', None, 'nm, nm, points')]

optical_settings = [('optics_method', 'BL', ''),
                    ('wavelength', '300, 1800, 251', 'nm, nm, points'),
                    ('source_type', 'standard', ''),
                    ('version', 'AM1.5g', ''),
                    ('output_units', 'photon_flux_per_m', '')]


class popupWindow(tk.Toplevel):
    def __init__(self, master, text, value, units):
        super(popupWindow, self).__init__(master)

        self.new_value = tk.StringVar(value=value)

        frame = ttk.Frame(self, padding=(10, 10, 10, 10))
        frame.grid()
        ttk.Label(frame, text=text).grid(column=0, row=0, columnspan=2, sticky=tk.NSEW)
        ttk.Entry(frame, textvariable=self.new_value).grid(column=0, row=1, sticky=tk.NSEW)
        ttk.Label(frame, text=units).grid(column=1, row=1, sticky=tk.NSEW)
        ttk.Button(frame, text='Ok', command=self.cleanup).grid(column=0, row=2, sticky=tk.NSEW)
        ttk.Button(frame, text='Cancel', command=self.exit).grid(column=1, row=2, sticky=tk.NSEW)

    def cleanup(self):
        self.value = self.new_value.get()
        self.destroy()

    def exit(self):
        self.value = None
        self.destroy()


class SettingsWindow(tk.Toplevel):
    def __init__(self, master, settings, model):
        super(SettingsWindow, self).__init__(master)

        frame = ttk.Frame(self, padding=(10, 10, 10, 10))
        frame.grid()

        self.settings = settings
        self.model = model

        # Junction properties
        self.settings_list = ttk.Treeview(frame)
        self.settings_list.grid(column=0, row=0, columnspan=2, sticky=(tk.NSEW))

        self.settings_list['columns'] = ('value', 'units')
        self.settings_list.heading('#0', text='Variable', anchor='w')
        self.settings_list.column("#0", minwidth=0, width=180, stretch=tk.NO)
        self.settings_list.heading('value', text='Value', anchor='w')
        self.settings_list.column("value", minwidth=0, width=180, stretch=tk.NO)
        self.settings_list.heading('units', text='Units', anchor='w')
        self.settings_list.column("units", minwidth=0, width=180, stretch=tk.NO)
        self.settings_list.tag_configure('modified', foreground='red')

        # When selecting an item, we update the properties frame
        self.settings_list.bind('<Double-Button-1>', self.settings_popup)

        ttk.Button(frame, text='Ok', command=self.cleanup).grid(column=0, row=2, sticky=tk.NSEW)
        ttk.Button(frame, text='Cancel', command=self.exit).grid(column=1, row=2, sticky=tk.NSEW)

        if settings == 'Global':
            self.load_global_settings()

        elif settings == 'Electrical':
            self.load_electrical_settings()

        elif settings == 'Optical':
            self.load_optical_settings()
        else:
            self.exit()

    def load_global_settings(self):
        """

        :return:
        """
        for value in global_settings:
            if value[0] in self.model['Global'].keys():
                self.settings_list.insert('', 'end', text=value[0], values=(self.model['Global'][value[0]], value[2]),
                                          tags=('modified'))
            else:
                self.settings_list.insert('', 'end', text=value[0], values=(value[1], value[2]))

    def load_electrical_settings(self):
        """

        :return:
        """
        for value in electrical_settings:
            if value[0] in self.model['Electrical'].keys():
                self.settings_list.insert('', 'end', text=value[0],
                                          values=(self.model['Electrical'][value[0]], value[2]),
                                          tags=('modified'))
            else:
                self.settings_list.insert('', 'end', text=value[0], values=(value[1], value[2]))

        for key in pdd_options.keys():
            if key in self.model['Electrical'].keys():
                self.settings_list.insert('', 'end', text=key, values=(self.model['Electrical'][key], ''),
                                          tags=('modified'))
            else:
                self.settings_list.insert('', 'end', text=key, values=(pdd_options[key], ''))

        for key in asc_options.keys():
            if key in self.model['Electrical'].keys():
                self.settings_list.insert('', 'end', text=key, values=(self.model['Electrical'][key], ''),
                                          tags=('modified'))
            else:
                self.settings_list.insert('', 'end', text=key, values=(asc_options[key], ''))

    def load_optical_settings(self):
        """

        :return:
        """
        for value in optical_settings:
            if value[0] in self.model['Optical'].keys():
                self.settings_list.insert('', 'end', text=value[0], values=(self.model['Optical'][value[0]], value[2]),
                                          tags=('modified'))
            else:
                self.settings_list.insert('', 'end', text=value[0], values=(value[1], value[2]))

        for key in rcwa_options.keys():
            if key in self.model['Optical'].keys():
                self.settings_list.insert('', 'end', text=key, values=(self.model['Optical'][key], ''),
                                          tags=('modified'))
            else:
                self.settings_list.insert('', 'end', text=key, values=(rcwa_options[key], ''))

    def settings_popup(self, *args):
        """

        :param args:
        :return:
        """
        item_id = self.settings_list.focus()
        text = self.settings_list.item(item_id, "text")
        value = self.settings_list.set(item_id, "value")
        units = self.settings_list.set(item_id, "units")

        w = popupWindow(self.master, 'New value for: ' + text, value, units)
        w.grab_set()
        self.master.wait_window(w)

        if w.value is not None and w.value != value:
            self.settings_list.set(item_id, "value", w.value)
            self.settings_list.item(item_id, tags=('modified'))

    def cleanup(self):
        """

        :return:
        """

        tagged = self.settings_list.tag_has('modified')

        for item_id in tagged:
            text = self.settings_list.item(item_id, "text")
            value = self.settings_list.set(item_id, "value")

            self.model[self.settings][text] = value

        self.destroy()

    def exit(self):
        """

        :return:
        """
        self.destroy()
