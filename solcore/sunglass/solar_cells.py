import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox

import yaml, os, copy, glob, time

import numpy as np
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt

import solcore
from solcore import material
from solcore.solar_cell import SolarCell
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell_solver import solar_cell_solver

from .properties_and_settings import *
from .run_calculation import run_solar_cell_model
from .save_data import save_bandstructure_data, save_iv_data, save_qe_data
from .plot import BandstructurePlot, IVPlot, QEPlot, save_plot, load_plot


class SolarCellsTab(ttk.Frame):
    """ Class that contains all widgets related to the creation and simulation of solar cells

    """

    def __init__(self, parent, master):
        """ Constructor of the class

        :param parent: The notebook that serves as parent of this tab
        """
        super(SolarCellsTab, self).__init__(parent)

        self.master = master

        self.create_buttons()
        self.create_solar_cell_list()
        self.create_layer_properties_frame()
        self.create_junction_properties_frame()
        self.create_solar_cell_properties_frame()
        self.create_output_frame()

        # The following dictionary contains all the information related to the simulation, except material properties
        # taken as the default values.
        self.model = {'Global': {'T': 298},
                      'Electrical': {},
                      'Optical': {},
                      'Solar cell': {}}

        self.current_filename = 'Untitled.sol'
        self.current_location = os.path.expanduser('~')
        self.wd_label_var = tk.StringVar(value=self.current_location)
        wd_label = ttk.Label(self, textvar=self.wd_label_var, relief='sunken')
        wd_label.grid(column=0, row=500, columnspan=500, sticky=tk.NSEW, pady=5, padx=5)

        self.unselect_item()

    def create_solar_cell_list(self):
        """

        :return: None
        """

        self.solar_cell_list = ttk.Treeview(self)
        solar_list_scroll = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.solar_cell_list.yview)
        self.solar_cell_list.configure(yscrollcommand=solar_list_scroll.set)
        self.solar_cell_list.grid(column=2, row=0, columnspan=3, rowspan=99, sticky=tk.NSEW)
        solar_list_scroll.grid(column=1, row=0, sticky=tk.NS)

        # To un-select items an item that was selected
        self.solar_cell_list.bind("<Button-1>", self.unselect_item)
        self.solar_cell_list.bind("<Escape>", self.unselect_item)

        # When selecting an item, we update the properties frame
        self.solar_cell_list.bind("<<TreeviewSelect>>", self.show_properties_frame)

        self.solar_cell_list['columns'] = ('type', 'width', 'material', 'options')
        self.solar_cell_list.heading('#0', text='Name/Role', anchor='w')
        self.solar_cell_list.column("#0", minwidth=0, width=150, stretch=tk.NO)
        self.solar_cell_list.heading('type', text='Type', anchor='w')
        self.solar_cell_list.column("type", minwidth=0, width=120, stretch=tk.NO)
        self.solar_cell_list.heading('width', text='width (nm)', anchor='w')
        self.solar_cell_list.column("width", minwidth=0, width=80, stretch=tk.NO)
        self.solar_cell_list.heading('material', text='material', anchor='w')
        self.solar_cell_list.column("material", minwidth=0, width=80, stretch=tk.NO)
        self.solar_cell_list.heading('options', text='options', anchor='w')
        self.solar_cell_list.column("options", minwidth=0, width=200, stretch=tk.YES)

    def unselect_item(self, *args):
        """ Un-select the item from the list

        :param args:
        :return:
        """
        self.solar_cell_list.selection_remove(self.solar_cell_list.focus())
        self.solar_cell_list.focus('')

        # If no item is selected, we have to disable the layer properties frame
        self.layer_properties_frame.grid_remove()
        self.sc_properties_frame.grid()

    def create_buttons(self):
        """

        :return: None
        """
        # Buttons for loading and saving data
        loadsave_frame = ttk.LabelFrame(self, text='Load/Save model')
        loadsave_frame.grid(column=0, row=0, sticky=tk.NSEW)
        loadsave_frame.columnconfigure(0, weight=1)

        save = ttk.Button(loadsave_frame, text='Save', command=self.save_model)
        load = ttk.Button(loadsave_frame, text='Load', command=self.load_model)
        clear = ttk.Button(loadsave_frame, text='Clear', command=self.clear_solar_cell)
        save.grid(column=0, row=0, sticky=tk.NSEW)
        load.grid(column=0, row=1, sticky=tk.NSEW)
        clear.grid(column=0, row=2, sticky=tk.NSEW)

        # Buttons for adding the structure
        structure_frame = ttk.LabelFrame(self, text='Structure')
        structure_frame.grid(column=0, row=1, sticky=tk.NSEW)

        layer = ttk.Button(structure_frame, text='Add layer', width=15, command=self.add_layer)
        junction = ttk.Button(structure_frame, text='Add junction', width=15, command=self.add_junction)
        up = ttk.Button(structure_frame, text='▲', width=2, command=self.move_up)
        out = ttk.Button(structure_frame, text='×', width=2, command=self.delete_line)
        dup = ttk.Button(structure_frame, text='⥀', width=2, command=self.duplicate)
        down = ttk.Button(structure_frame, text='▼', width=2, command=self.move_down)
        layer.grid(column=0, row=0, columnspan=4, sticky=tk.NSEW)
        junction.grid(column=0, row=1, columnspan=4, sticky=tk.NSEW)
        up.grid(column=0, row=3, sticky=tk.NSEW)
        out.grid(column=1, row=3, sticky=tk.NSEW)
        dup.grid(column=2, row=3, sticky=tk.NSEW)
        down.grid(column=3, row=3, sticky=tk.NSEW)

        # Buttons for setting the conditions
        settings_frame = ttk.LabelFrame(self, text='Settings')
        settings_frame.grid(column=0, row=2, sticky=(tk.NSEW))
        settings_frame.columnconfigure(0, weight=1)

        glob = ttk.Button(settings_frame, text='Global', command=lambda: self.settings_popup(settings='Global'))
        elect = ttk.Button(settings_frame, text='Electrical',
                           command=lambda: self.settings_popup(settings='Electrical'))
        opt = ttk.Button(settings_frame, text='Optical', command=lambda: self.settings_popup(settings='Optical'))
        glob.grid(column=0, row=0, columnspan=1, sticky=tk.NSEW)
        elect.grid(column=0, row=1, columnspan=1, sticky=tk.NSEW)
        opt.grid(column=0, row=2, columnspan=1, sticky=tk.NSEW)

        # Buttons for running the solver
        running_frame = ttk.LabelFrame(self, text='Run')
        running_frame.grid(column=0, row=3, sticky=(tk.NSEW))
        running_frame.columnconfigure(0, weight=1)

        self.run_var = tk.IntVar(value=2)
        iv = ttk.Radiobutton(running_frame, text="IV", variable=self.run_var, value=2)
        qe = ttk.Radiobutton(running_frame, text="QE", variable=self.run_var, value=3)
        equi = ttk.Radiobutton(running_frame, text="Equilibrium", variable=self.run_var, value=4)
        sc = ttk.Radiobutton(running_frame, text="Short Circuit", variable=self.run_var, value=5)
        iv.grid(column=0, row=1, sticky=tk.EW)
        qe.grid(column=0, row=2, sticky=tk.EW)
        equi.grid(column=0, row=3, sticky=tk.EW)
        sc.grid(column=0, row=4, sticky=tk.EW)

        run = ttk.Button(running_frame, text='Run', command=self.run)
        run.grid(column=0, row=5, sticky=tk.NSEW)

    def move_up(self):
        """ Moves up an item

        :return: None
        """
        item_id = self.solar_cell_list.focus()
        self.solar_cell_list.move(item_id, self.solar_cell_list.parent(item_id),
                                  self.solar_cell_list.index(item_id) - 1)

    def move_down(self):
        """ Moves down an item

        :return: None
        """

        item_id = self.solar_cell_list.focus()
        self.solar_cell_list.move(item_id, self.solar_cell_list.parent(item_id),
                                  self.solar_cell_list.index(item_id) + 1)

    def delete_line(self):
        """ Deletes an item and moves the focus to the previous one within the same parent. When there are no more previous items, the command is ignored.

        :return: None
        """
        try:
            item_id = self.solar_cell_list.focus()
            other_item = self.solar_cell_list.prev(item_id)
            self.solar_cell_list.delete(item_id)
            self.model['Solar cell'].pop(item_id)
            self.solar_cell_list.focus(other_item)
        except:
            return

    def duplicate(self):
        """ Duplicates an item. Currently, it only works well for single layers as it does not duplicate the contents of a junction or tunnel junction, if present.

        :return: None
        """
        item_id = self.solar_cell_list.focus()
        parent_id = self.solar_cell_list.parent(item_id)
        item = self.solar_cell_list.item(item_id)
        new = self.solar_cell_list.insert(parent_id, tk.END, text=item['text'], values=item['values'], open=True)
        self.model['Solar cell'][new] = copy.copy(self.model['Solar cell'][item_id])

        # And we duplicate the children, if any
        for c in self.solar_cell_list.get_children(item_id):
            child_item = self.solar_cell_list.item(c)
            new_child = self.solar_cell_list.insert(new, tk.END, text=child_item['text'],
                                                    values=copy.copy(child_item['values']))
            self.model['Solar cell'][new_child] = copy.copy(self.model['Solar cell'][c])

    def add_layer(self):
        """ Adds a layer to the structure, either in the root or within a junction. The layer is also added to the model dictionary

        :return: None
        """
        item_id = self.solar_cell_list.focus()
        if 'Layer' not in self.solar_cell_list.item(item_id)['values']:
            new = self.solar_cell_list.insert(item_id, tk.END, text='New layer',
                                              values=('Layer', 100.0, 'GaAs', ''))
            self.model['Solar cell'][new] = {'name': 'New layer', 'type': 'Layer', 'width': 100.0,
                                             'material': 'GaAs',
                                             'options': {}}

    def add_junction(self):
        """ Adds a junction to the structure

        :return: None
        """
        item_id = self.solar_cell_list.insert('', tk.END, text='New junction', values=('Junction-PDD', '', '', ''),
                                              open=True)
        self.model['Solar cell'][item_id] = {'name': 'New junction', 'type': 'Junction-PDD', 'width': '',
                                             'material': '',
                                             'options': {}}

    def clear_solar_cell(self, confirmed=False):
        """ Delets all elements in the slar cell

        :return:
        """

        clear = True
        if not confirmed:
            clear = messagebox.askyesno("Are you sure you want to delete current model?",
                                        "This action cannot be undone!",
                                        icon='warning')

        if clear:
            old_list = self.solar_cell_list.get_children()
            if len(old_list) > 0:
                self.solar_cell_list.delete(*old_list)

            self.model = {'Global': {'T': 298},
                          'Electrical': {},
                          'Optical': {},
                          'Solar cell': {}}

            self.current_filename = 'Untitled.sol'
            self.current_location = os.path.expanduser('~')

    def settings_popup(self, settings):
        """

        :return:
        """
        w = SettingsWindow(self.master, settings, self.model)
        w.grab_set()
        self.master.wait_window(w)

    def properties_popup(self, *args):
        """

        :return:
        """
        element_id = self.solar_cell_list.focus()

        if self.model['Solar cell'][element_id]['type'] == 'Layer':
            # We are editing the properties of a layer
            widget = self.properties_list
        else:
            # We are editing the properties of a junction
            widget = self.junction_properties_list

        item_id = widget.focus()
        text = widget.item(item_id, "text")
        value = widget.set(item_id, "value")
        units = widget.set(item_id, "units")

        w = popupWindow(self.master, 'New value for: ' + text, value, units)
        w.grab_set()
        self.master.wait_window(w)

        if w.value is not None and float(w.value) != value:
            widget.set(item_id, "value", float(w.value))
            widget.item(item_id, tags=('modified'))
            self.model['Solar cell'][element_id]['options'][text] = float(w.value)

    def create_solar_cell_properties_frame(self, *args):
        """ Creates the frame that contain the general solar cell properties, like name, substrate, front relfection...

        :param args:
        :return:
        """
        self.sc_properties_frame = ttk.LabelFrame(self, text='Solar cell properties', width=380)
        self.sc_properties_frame.grid(column=10, row=0, rowspan=30, sticky=tk.NSEW)
        self.sc_properties_frame.columnconfigure(1, weight=1)
        self.sc_properties_frame.grid_propagate(0)

        # Name of the cell
        ttk.Label(self.sc_properties_frame, text="Name:").grid(column=0, row=0, sticky=tk.NSEW)
        self.sc_name_var = tk.StringVar(value='New solar cell')
        sc_name = ttk.Entry(self.sc_properties_frame, textvariable=self.sc_name_var)
        sc_name.grid(column=1, row=0, sticky=tk.NSEW)

        # Substrate material
        ttk.Label(self.sc_properties_frame, text="Substrate:").grid(column=0, row=1, sticky=(tk.NSEW))
        self.substrate_var = tk.StringVar(value='GaAs')
        substrate_box = ttk.Combobox(self.sc_properties_frame, state='readonly',
                                     values=self.master.materials.valid_materials,
                                     textvariable=self.substrate_var)
        substrate_box.grid(column=1, row=1, sticky=tk.NSEW)
        substrate_box.bind('<<ComboboxSelected>>', self.update_substrate_composition_label)

        # Composition
        self.subs_composition_label_var = tk.StringVar(value='-')
        lab = ttk.Label(self.sc_properties_frame, textvariable=self.subs_composition_label_var)
        lab.grid(column=0, row=2, sticky=tk.NSEW)
        self.subs_composition_var = tk.DoubleVar(value=0)
        subs_composition = tk.Spinbox(master=self.sc_properties_frame, from_=0, to=1,
                                      textvariable=self.subs_composition_var)
        subs_composition.grid(column=1, row=2, sticky=tk.NSEW)

        # Reflectivity
        self.reflectivity_check = tk.IntVar(value=0)
        reflect = ttk.Checkbutton(self.sc_properties_frame, text='Reflectivity:', var=self.reflectivity_check,
                                  command=self.choose_reflectivity)
        reflect.grid(column=0, row=3, sticky=tk.NSEW)
        self.reflectivity_var = tk.StringVar(value='')
        reflectivity = ttk.Entry(self.sc_properties_frame, textvariable=self.reflectivity_var)
        reflectivity.grid(column=1, row=3, sticky=tk.NSEW)

        # Shading
        ttk.Label(self.sc_properties_frame, text="Shading").grid(column=0, row=4, sticky=tk.NSEW)
        self.shading_var = tk.DoubleVar(value=0)
        shading = tk.Spinbox(master=self.sc_properties_frame, from_=0, to=1, textvariable=self.shading_var)
        shading.grid(column=1, row=4, sticky=tk.NSEW)

        # Size
        ttk.Label(self.sc_properties_frame, text="Size (cm2):").grid(column=0, row=5, sticky=tk.NSEW)
        self.size_var = tk.DoubleVar(value=1)
        size = tk.Spinbox(self.sc_properties_frame, textvariable=self.size_var, from_=0.0001, to=1e6)
        size.grid(column=1, row=5, sticky=tk.NSEW)

        # Series resistance
        ttk.Label(self.sc_properties_frame, text="R series (Ohm):").grid(column=0, row=6, sticky=tk.NSEW)
        self.rs_var = tk.DoubleVar(value=0)
        rs = tk.Spinbox(self.sc_properties_frame, textvariable=self.rs_var, from_=0)
        rs.grid(column=1, row=6, sticky=tk.NSEW)

    def choose_reflectivity(self, *args):
        """ Allows to choose the reflectivity of the solar cell

        :param args:
        :return:
        """
        if self.reflectivity_check.get() == 0:
            return
        else:
            # A pop up dialog allows to choose the relfectivity
            filename = filedialog.askopenfilename(initialdir=self.current_location)

            if filename is not '':
                self.reflectivity_var.set(filename)
            else:
                self.reflectivity_check.set(0)

    def update_substrate_composition_label(self, *args):
        """ Updates the composition label of the substrate with the relevant material name

        :param args:
        :return:
        """
        # Checking if the material is ternary (or more) by counting capital letters.
        # If it is not, we ignore the compositionon
        mat = self.substrate_var.get()
        if sum(1 for c in mat if c.isupper()) > 2:
            x_label = solcore.ParameterSystem().database.get(mat, 'x')
            self.subs_composition_label_var.set(x_label)
        else:
            self.subs_composition_label_var.set('-')

    def confirm_changes_to_solar_cell(self, *args):
        """ Confirm all changes done to the solar cell properties by storing them in the model

        :param args:
        :return:
        """
        self.model['Solar cell']['name'] = self.sc_name_var.get()
        self.model['Solar cell']['substrate'] = self.substrate_var.get()
        self.model['Solar cell']['r_series'] = self.rs_var.get()
        self.model['Solar cell']['reflectivity'] = self.reflectivity_var.get()
        self.model['Solar cell']['shading'] = self.shading_var.get()
        self.model['Solar cell']['size'] = self.size_var.get()

        if self.subs_composition_label_var.get() != '-':
            self.model['Solar cell']['element'] = self.subs_composition_label_var.get()
            self.model['Solar cell']['composition'] = self.subs_composition_var.get()

    def create_junction_properties_frame(self, *args):
        """ Creates the frame that contain the junction properties

        :param args:
        :return:
        """
        self.junction_properties_frame = ttk.LabelFrame(self, text='Junction properties', width=380)
        self.junction_properties_frame.grid(column=10, row=0, rowspan=30, sticky=tk.NSEW)
        self.junction_properties_frame.columnconfigure(1, weight=1)
        self.junction_properties_frame.grid_propagate(0)

        # Name of the cell
        ttk.Label(self.junction_properties_frame, text="Name:").grid(column=0, row=0, sticky=tk.NSEW)
        self.junction_name_var = tk.StringVar(value='New junction')
        junction_name = ttk.Entry(self.junction_properties_frame, textvariable=self.junction_name_var)
        junction_name.grid(column=1, row=0, sticky=tk.NSEW)

        # Type of junction
        keys = [k for k in properties_junctions.keys()]
        ttk.Label(self.junction_properties_frame, text="Type:").grid(column=0, row=1, sticky=(tk.NSEW))
        self.junction_type_var = tk.StringVar(value='Junction-PDD')
        junction_type_box = ttk.Combobox(self.junction_properties_frame, state='readonly',
                                         values=keys,
                                         textvariable=self.junction_type_var)
        junction_type_box.grid(column=1, row=1, sticky=tk.NSEW)
        junction_type_box.bind('<<ComboboxSelected>>', self.update_junction_properties_list)

        # Junction properties
        self.junction_properties_list = ttk.Treeview(self.junction_properties_frame)
        self.junction_properties_list.grid(column=0, row=3, columnspan=2, sticky=(tk.NSEW))

        self.junction_properties_list['columns'] = ('value', 'units')
        self.junction_properties_list.heading('#0', text='Variable', anchor='w')
        self.junction_properties_list.column("#0", minwidth=0, width=180, stretch=tk.NO)
        self.junction_properties_list.heading('value', text='Value', anchor='w')
        self.junction_properties_list.column("value", minwidth=0, width=100, stretch=tk.NO)
        self.junction_properties_list.heading('units', text='Units', anchor='w')
        self.junction_properties_list.column("units", minwidth=0, width=100, stretch=tk.NO)
        self.junction_properties_list.tag_configure('modified', foreground='red')

        # When selecting an item, we update the properties frame
        self.junction_properties_list.bind('<Double-Button-1>', self.properties_popup)

        confirm = ttk.Button(self.junction_properties_frame, text='Confirm changes to junction',
                             command=self.confirm_changes_to_junction)
        confirm.grid(column=0, row=99, columnspan=2, sticky=tk.NSEW)

    def update_junction_properties_list(self, *args):
        """

        :return:
        """

        item_id = self.solar_cell_list.focus()

        self.model['Solar cell'][item_id]['type'] = self.junction_type_var.get()

        # We populate the properties TreeView with the values for this junction,
        # starting by eliminating the current ones
        old_list = self.junction_properties_list.get_children()
        if len(old_list) > 0:
            self.junction_properties_list.delete(*old_list)

        # Get the new set of parameters
        prop_list = self.get_properties_list(item_id)

        # We scan the list of required properties and load the overwritten ones of the default ones
        for key in prop_list:
            try:
                prop = self.model['Solar cell'][item_id]['options'][key]
                self.junction_properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]),
                                                     tags=('modified'))
            except:
                prop = default_junction_properties[key] / conversion[key]

                self.junction_properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]))

    def confirm_changes_to_junction(self, *args):
        """

        :param args:
        :return:
        """
        item_id = self.solar_cell_list.focus()

        junction_type = self.model['Solar cell'][item_id]['name'] = self.junction_type_var.get()
        name = self.model['Solar cell'][item_id]['name'] = self.junction_name_var.get()

        options = ''
        for key in self.model['Solar cell'][item_id]['options']:
            options += '{0}={1} '.format(key, self.model['Solar cell'][item_id]['options'][key])

        # We update the layers list
        self.solar_cell_list.item(item_id, text=name, values=[junction_type, '', '', options])

    def create_layer_properties_frame(self, *args):
        """

        :return: None
        """

        name = 'New layer'
        mat = 'GaAs'
        width = 100

        self.layer_properties_frame = ttk.LabelFrame(self, text='Layer properties', width=380)
        self.layer_properties_frame.grid(column=10, row=0, rowspan=30, sticky=tk.NSEW)
        self.layer_properties_frame.columnconfigure(1, weight=1)
        self.layer_properties_frame.rowconfigure(4, weight=1)
        self.layer_properties_frame.grid_propagate(0)

        ttk.Label(self.layer_properties_frame, text="Name/Role:").grid(column=0, row=0, sticky=(tk.NSEW))
        self.name_var = tk.StringVar(value=name)
        self.name = ttk.Entry(self.layer_properties_frame, textvariable=self.name_var)
        self.name.grid(column=1, row=0, sticky=(tk.NSEW))

        ttk.Label(self.layer_properties_frame, text="Width (nm):").grid(column=0, row=1, sticky=(tk.NSEW))
        self.width_var = tk.DoubleVar(value=width)
        self.width = tk.Spinbox(master=self.layer_properties_frame, from_=0, textvariable=self.width_var)
        self.width.grid(column=1, row=1, sticky=(tk.NSEW))

        ttk.Label(self.layer_properties_frame, text="Material:").grid(column=0, row=2, sticky=(tk.NSEW))
        self.selected_material = tk.StringVar(value=mat)
        self.materials_box = ttk.Combobox(self.layer_properties_frame, state='readonly',
                                          values=self.master.materials.valid_materials,
                                          textvariable=self.selected_material)
        self.materials_box.grid(column=1, row=2, sticky=(tk.NSEW))
        self.materials_box.bind('<<ComboboxSelected>>', self.update_composition_label)

        self.composition_label_var = tk.StringVar(value='-')
        ttk.Label(self.layer_properties_frame, textvariable=self.composition_label_var).grid(column=0, row=3,
                                                                                             sticky=(tk.NSEW))
        self.composition_var = tk.DoubleVar(value=0)
        self.composition = tk.Spinbox(master=self.layer_properties_frame, from_=0, to=1,
                                      textvariable=self.composition_var)
        self.composition.grid(column=1, row=3, sticky=(tk.NSEW))

        self.properties_list = ttk.Treeview(self.layer_properties_frame)
        self.properties_list.grid(column=0, row=4, columnspan=2, sticky=(tk.NSEW))

        self.properties_list['columns'] = ('value', 'units')
        self.properties_list.heading('#0', text='Variable', anchor='w')
        self.properties_list.column("#0", minwidth=0, width=180, stretch=tk.NO)
        self.properties_list.heading('value', text='Value', anchor='w')
        self.properties_list.column("value", minwidth=0, width=100, stretch=tk.NO)
        self.properties_list.heading('units', text='Units', anchor='w')
        self.properties_list.column("units", minwidth=0, width=100, stretch=tk.NO)
        self.properties_list.tag_configure('modified', foreground='red')

        # When selecting an item, we update the properties frame
        self.properties_list.bind('<Double-Button-1>', self.properties_popup)

        confirm = ttk.Button(self.layer_properties_frame, text='Confirm changes to layer',
                             command=self.confirm_changes_to_layer)
        confirm.grid(column=0, row=99, columnspan=2, sticky=tk.NSEW)

    def show_properties_frame(self, *args):
        """

        :return: None
        """

        item_id = self.solar_cell_list.focus()

        if item_id == '':
            return
        elif 'Layer' in self.model['Solar cell'][item_id]['type']:
            self.sc_properties_frame.grid_remove()
            self.junction_properties_frame.grid_remove()
            self.show_layer_properties_frame(item_id)
        elif 'Junction' in self.model['Solar cell'][item_id]['type']:
            self.sc_properties_frame.grid_remove()
            self.layer_properties_frame.grid_remove()
            self.show_junction_properties_frame(item_id)

    def show_junction_properties_frame(self, item_id):
        """

        :return: None
        """
        current_item = self.model['Solar cell'][item_id]

        # We bring the properties frame back, in case it was hidden
        self.junction_properties_frame.grid()

        self.junction_name_var.set(current_item['name'])
        self.junction_type_var.set(current_item['type'])

        # Finally, we populate the properties TreeView with the values for this junction,
        # starting by eliminating the current ones
        old_list = self.junction_properties_list.get_children()
        if len(old_list) > 0:
            self.junction_properties_list.delete(*old_list)

        # Get the new set of parameters
        prop_list = self.get_properties_list(item_id)

        # We scan the list of required properties and load the overwritten ones of the default ones
        for key in prop_list:
            try:
                prop = current_item['options'][key]
                self.junction_properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]),
                                                     tags=('modified'))
            except:
                prop = default_junction_properties[key] / conversion[key]

                self.junction_properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]))

    def show_layer_properties_frame(self, item_id):
        """

        :return: None
        """

        current_item = self.model['Solar cell'][item_id]

        # We bring the properties frame back, in case it was hidden
        self.layer_properties_frame.grid()

        # # And enable the input
        # self.enable_disable_layer_properties_frame('enable')

        self.name_var.set(current_item['name'])
        self.width_var.set(current_item['width'])
        self.selected_material.set(current_item['material'])

        # Checking if the material is ternary (or more) by counting capital letters.
        # If it is not, we ignore the composition
        if sum(1 for c in current_item['material'] if c.isupper()) > 2:
            try:
                x_label = current_item['options']['element']
            except:
                x_label = solcore.ParameterSystem().database.get(current_item['material'], 'x')
                current_item['options']['element'] = x_label

            self.composition_label_var.set(x_label)
            x = current_item['options']['x'] if 'x' in current_item['options'].keys() else 0
            self.composition_var.set(x)
        else:
            x_label = '-'
            x = 0
            self.composition_label_var.set('-')
            self.composition_var.set(0)

        # Finally, we populate the properties TreeView with the values for this layer,
        # starting by eliminating the current ones
        old_list = self.properties_list.get_children()
        if len(old_list) > 0:
            self.properties_list.delete(*old_list)

        # Get the new set of parameters
        prop_list = self.get_properties_list(item_id)

        Na = current_item['options']['Na'] if 'Na' in current_item['options'].keys() else 0.0
        Nd = current_item['options']['Nd'] if 'Nd' in current_item['options'].keys() else 0.0

        # Most properties are the default ones defined for that material, so we need to open that material
        solcore_mat = solcore.material(current_item['material'])(
            **{'T': self.model['Global']['T'], x_label: x, 'Na': Na, 'Nd': Nd})

        # We scan the list of required properties and load the overwritten ones of the default ones
        for key in prop_list:
            try:
                prop = current_item['options'][key]
                self.properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]),
                                            tags=('modified'))
            except:
                try:
                    # We try the default property fo that material
                    prop = solcore_mat.__getattr__(key) / conversion[key]
                    self.properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]))
                except:
                    # Or a default value.
                    prop = default_layer_properties[key] / conversion[key]
                    current_item['options'][key] = prop
                    self.properties_list.insert('', 'end', text=key, values=(format(prop, '.4'), units[key]),
                                                tags=('modified'))

    def get_properties_list(self, item_id):
        """ Return the properties that a layer needs to have. Which properties are relevant depend on where the layer is located, if its outside any juction, if the junction is a PDD one, a DA, etc.

        :param item_id: The id of the item we want the properties
        :return: The list of relevant properties
        """
        item = self.model['Solar cell'][item_id]
        parent_id = self.solar_cell_list.parent(item_id)
        try:
            parent = self.model['Solar cell'][parent_id]
            return properties_layers[parent['type']]
        except KeyError:
            if 'Layer' in item['type']:
                return []
            else:
                return properties_junctions[item['type']]

    def confirm_changes_to_layer(self, *args):
        """ This command confirm the changes made to the layer, validating the values and updating the TreeView and the Model

        :param args: Dummy
        :return: None
        """

        item_id = self.solar_cell_list.focus()

        if item_id == '':
            return

        name = self.name_var.get()
        mat = self.selected_material.get()
        width = self.width_var.get()
        width = max(width, 0)  # <-- In case width is negative
        self.width_var.set(width)

        self.model['Solar cell'][item_id]['name'] = name
        self.model['Solar cell'][item_id]['width'] = width
        self.model['Solar cell'][item_id]['material'] = mat

        # Checking if the material is ternary (or more) by counting capital letters.
        # If it is not, we ignore the composition
        if sum(1 for c in mat if c.isupper()) == 3:
            x = self.composition_var.get()
            x = min(max(x, 0), 1)  # <-- In case x is outside the 0-1 range
            self.composition_var.set(x)
            self.model['Solar cell'][item_id]['options']['x'] = x
            self.model['Solar cell'][item_id]['options']['element'] = self.composition_label_var.get()
        else:
            try:
                self.model['Solar cell'][item_id]['options'].pop('x')
                self.model['Solar cell'][item_id]['options'].pop('element')
            except:
                pass

        options = ''
        for key in self.model['Solar cell'][item_id]['options']:
            options += '{0}={1} '.format(key, self.model['Solar cell'][item_id]['options'][key])

        # We update the layers list
        self.solar_cell_list.item(item_id, text=name, values=['Layer', width, mat, options])

        # And the whole properties frame. Part of what is done here is redundant and can be improved a lot
        self.show_properties_frame()

    def update_composition_label(self, *args):
        """ Updates the composition label with the relevant material name

        :param args:
        :return:
        """
        # Checking if the material is ternary (or more) by counting capital letters.
        # If it is not, we ignore the compositionon
        mat = self.materials_box.get()
        if sum(1 for c in mat if c.isupper()) > 2:
            x_label = solcore.ParameterSystem().database.get(mat, 'x')
            self.composition_label_var.set(x_label)
        else:
            self.composition_label_var.set('-')

    def get_solar_cell_structure(self):
        """ Gets the structure of the solar cell from the tree view (order of the layers, junctions and layers within each junction).

        :return:
        """
        # First we get the structure of the solar cell, just using the labels. The first label in each list is the
        # parent
        structure = []
        children = self.solar_cell_list.get_children()
        for c in children:
            structure.append([c])
            grandchildren = self.solar_cell_list.get_children(c)
            for gc in grandchildren:
                structure[-1].append(gc)

        self.model['Solar cell']['structure'] = structure
        self.confirm_changes_to_solar_cell()

    def save_model(self):
        """ Saves the solar cell model to a Yaml-formatted file

        :return: None
        """

        self.get_solar_cell_structure()

        try:
            filename = filedialog.asksaveasfilename(defaultextension='.sol', initialdir=self.current_location,
                                                    initialfile=self.current_filename)
            with open(filename, 'w') as f:
                yaml.dump(self.model, f)

            self.current_location, self.current_filename = os.path.split(filename)
            self.wd_label_var.set(self.current_location)
        except:
            return 1

    def load_model(self):
        """ Loads a solar cell model previously saved in a Yalm-formatted file

        :return: None
        """

        # First we eliminate the existing elements of the solar cell list
        self.clear_solar_cell(confirmed=True)

        # We load the data from a Solcore file
        try:
            filename = filedialog.askopenfilename(initialdir=self.current_location,
                                                  defaultextension='.sol',
                                                  filetypes=[('All files', '*.*'), ('Solcore', '*.sol')],
                                                  initialfile=self.current_filename)
            with open(filename, 'r') as f:
                self.model = yaml.load(f)

            self.current_location, self.current_filename = os.path.split(filename)
            self.wd_label_var.set(self.current_location)
            self.refresh_plots_list()
        except:
            return

        try:
            # And now we populate the list with the loaded elements: junctions, tunnel junctions or layers
            for c in self.model['Solar cell']['structure']:
                item_id = c[0]
                current_item = self.model['Solar cell'][item_id]

                options = ''
                for key in current_item['options']:
                    options += '{0}={1} '.format(key, current_item['options'][key])

                self.solar_cell_list.insert('', tk.END, iid=item_id, text=current_item['name'], open=True,
                                            values=(
                                                current_item['type'], current_item['width'], current_item['material'],
                                                options))

                # Now we populate the layers within the junctions, if any
                for id in c[1:]:
                    layer = self.model['Solar cell'][id]

                    options = ''
                    for key in layer['options']:
                        options += '{0}={1} '.format(key, layer['options'][key])

                    self.solar_cell_list.insert(item_id, tk.END, iid=id, text=layer['name'], open=True,
                                                values=(layer['type'], layer['width'], layer['material'], options))
        except KeyError:
            print("The model cannot be loaded. There is missing information or the format is not correct.")

    def run(self, *args):
        """ Runs the selected calculation, starting by saving the solar cell model

        :param args:
        :return:
        """

        # We get the solar cell structure and if there is no structure, we do nothing else.
        self.get_solar_cell_structure()
        if len(self.model['Solar cell']['structure']) == 0:
            return

        # Get problem to solve
        list_of_tasks = ['optics', 'iv', 'qe', 'equilibrium', 'short_circuit']
        task = list_of_tasks[self.run_var.get() - 1]

        # Run model
        output = run_solar_cell_model(task, self.model)

        # Get output
        # And now we save the data to individual files for IV, QE and bandstructure (only for PDD junctions)
        # And plot the results
        # And show the calculated outputs, if any
        root_filename = os.path.join(self.current_location, self.current_filename.split('.')[0])
        if task == 'iv':
            light_iv = self.model['Electrical']['light_iv'] if 'light_iv' in self.model['Electrical'].keys() else False
            df, filename = save_iv_data(root_filename, output, light_iv)
            plot_name = '{} {}'.format(os.path.basename(filename.split('.')[0]), time.asctime())
            plot = IVPlot(df, light=light_iv, title=plot_name)
            save_plot(plot, filename.replace('csv', 'pkl'))
            self.show_plot_outputs(plot)

        elif task == 'qe':
            df, filename = save_qe_data(root_filename, output)
            plot_name = '{} {}'.format(os.path.basename(filename.split('.')[0]), time.asctime())
            plot = QEPlot(df, title=plot_name)
            save_plot(plot, filename.replace('csv', 'pkl'))
            self.show_plot_outputs(plot)

        elif task in ['equilibrium', 'short_circuit']:
            print('WARNING: Only data for PDD junctions will be saved.')
            df, filename = save_bandstructure_data(root_filename, output, task)
            plot_name = '{} {}'.format(os.path.basename(filename.split('.')[0]), time.asctime())
            plot = BandstructurePlot(df, title=plot_name)
            save_plot(plot, filename.replace('csv', 'pkl'))
            
        elif task == 'optics':
            pass

        self.refresh_plots_list()

    def create_output_frame(self):
        """

        :return:
        """

        # Frame containing buttons controlling the output and a few buttons
        out_frame = ttk.LabelFrame(self, text='Plot actions')
        out_frame.grid(column=0, row=100, rowspan=99, sticky=tk.NSEW)
        out_frame.columnconfigure(0, weight=1)

        choose_directory_button = ttk.Button(out_frame, text='Choose working directory', command=self.choose_directory)
        choose_directory_button.grid(column=0, row=0, sticky=tk.NSEW)
        refresh_plots_button = ttk.Button(out_frame, text='Refresh plots list', command=self.refresh_plots_list)
        refresh_plots_button.grid(column=0, row=1, sticky=tk.NSEW)

        # We create the list of figures available in the current directory
        self.figures_list = ttk.Treeview(self)
        figures_list_scroll = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.figures_list.yview)
        self.figures_list.configure(yscrollcommand=figures_list_scroll.set)
        self.figures_list.grid(column=2, row=100, columnspan=3, rowspan=99, sticky=tk.NSEW)
        figures_list_scroll.grid(column=1, row=100, sticky=tk.NS)

        # When selecting an item, we load the corresponding figure
        self.figures_list.bind("<Double-Button-1>", self.open_plot)
        self.figures_list.bind("<<TreeviewSelect>>", self.show_plot)

        self.figures_list['columns'] = ('type')
        self.figures_list.heading('#0', text='Filename', anchor='w')
        self.figures_list.column("#0", stretch=tk.YES)
        self.figures_list.heading('type', text='Plot type', anchor='w')
        self.figures_list.column("type", stretch=tk.YES)

        # And a list containing the calculated parameters, if any
        self.plot_outputs_list = ttk.Treeview(self)
        self.plot_outputs_list.grid(column=10, row=100, columnspan=2, sticky=(tk.NSEW))

        self.plot_outputs_list['columns'] = ('value', 'units')
        self.plot_outputs_list.heading('#0', text='Variable', anchor='w')
        self.plot_outputs_list.column("#0", minwidth=0, width=180, stretch=tk.NO)
        self.plot_outputs_list.heading('value', text='Value', anchor='w')
        self.plot_outputs_list.column("value", minwidth=0, width=100, stretch=tk.NO)
        self.plot_outputs_list.heading('units', text='Units', anchor='w')
        self.plot_outputs_list.column("units", minwidth=0, width=100, stretch=tk.NO)

    def show_plot(self, *args) -> None:
        """ Brings to the foregrund a figure already plotted

        :param args: dummy
        :return: None
        """
        id = self.figures_list.focus()

        if self.plot_outputs[id]['object'] is not None:
            raise_window(self.plot_outputs[id]['object'].fig.number)
            self.show_plot_outputs(self.plot_outputs[id]['object'])

    def open_plot(self, *args) -> None:
        """ Opens the selected plot. If it was already oppened, closes it - actually just destroys the figure - and opens it again.

        :return: None
        """
        id = self.figures_list.focus()

        if self.plot_outputs[id]['object'] is not None:
            plt.close(self.plot_outputs[id]['object'].fig.number)
            self.plot_outputs[id]['object'] = None

        self.plot_outputs[id]['object'] = load_plot(self.plot_outputs[id]['path'], plot=True)
        self.show_plot_outputs(self.plot_outputs[id]['object'])

    def show_plot_outputs(self, plot):
        """ Gets the output calculations of the plot, if any, and shows them

        :return:
        """
        if not hasattr(plot, 'outputs'):
            return

        # First we clear the existing list
        old_list = self.plot_outputs_list.get_children()
        if len(old_list) > 0:
            self.plot_outputs_list.delete(*old_list)

        for key in plot.outputs.keys():
            self.plot_outputs_list.insert('', tk.END, text=key, values=(format(plot.outputs[key], '.4'), ''))

    def refresh_plots_list(self) -> None:
        """

        :return:
        """
        # First we clear the existing list
        old_list = self.figures_list.get_children()
        if len(old_list) > 0:
            self.figures_list.delete(*old_list)

        # The we populate it with new data
        available_plots = glob.glob(os.path.join(self.current_location, '*.pkl'))
        self.plot_outputs = {}

        for f in available_plots:
            name = os.path.basename(f).split('.')[0]

            if 'QE' in name:
                fig_type = 'QE'
            elif 'IV' in name:
                fig_type = 'IV'
            elif any(a in name for a in ['equilibrium', 'short_circuit']):
                fig_type = 'Bandstructure'
            else:
                fig_type = 'Unknown'

            id = self.figures_list.insert('', tk.END, text=name, values=[fig_type])
            self.plot_outputs[id] = {'path': f, 'object': None}

    def choose_directory(self) -> None:
        """ Opens a pop-up window to choose the working directory

        :return: None
        """
        out = filedialog.askdirectory(initialdir=self.current_location)

        if out != '':
            self.current_location = out
            self.wd_label_var.set(self.current_location)
            self.refresh_plots_list()


def raise_window(figname=None):
    """
    Raise the plot window for Figure figname to the foreground.  If no argument
    is given, raise the current figure.

    This function will only work with a Tk graphics backend.  It assumes you
    have already executed the command 'import matplotlib.pyplot as plt'.

    Function taken from:
    http://physicalmodelingwithpython.blogspot.com/2015/07/raising-figure-window-to-foreground.html

    """

    if figname:
        plt.figure(figname)
    cfm = plt.get_current_fig_manager()
    cfm.window.attributes('-topmost', True)
    cfm.window.attributes('-topmost', False)
    return cfm
