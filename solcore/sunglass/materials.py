import tkinter as tk
from tkinter import ttk

import solcore


class MaterialsTab(ttk.Frame):
    """ Class that contains all widgets related to the managing of Solcore materials

    """
    def __init__(self, parent, master):
        """ Constructor of the class

        :param parent: The notebook that serves as parent of this tab
        """

        super(MaterialsTab, self).__init__(parent)

        self.master = master

        self.populate_materials_list()

    def populate_materials_list(self):
        """

        :return:
        """

        self.mat_list = tk.Listbox(self)
        self.mat_list.bind('<<ListboxSelect>>', self.load_material_properties)
        mat_list_scroll = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.mat_list.yview)
        self.mat_list.configure(yscrollcommand=mat_list_scroll.set)
        self.mat_list.grid(column=1, row=0, columnspan=3, sticky=(tk.NSEW))
        mat_list_scroll.grid(column=0, row=0, sticky=(tk.NS))

        self.mat_properties = ttk.Treeview(self)
        mat_properties_scroll = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.mat_properties.yview)
        self.mat_properties.configure(yscrollcommand=mat_properties_scroll.set)
        self.mat_properties.grid(column=5, row=0, columnspan=3, sticky=(tk.NSEW))
        mat_properties_scroll.grid(column=4, row=0, sticky=(tk.NS))

        self.mat_properties['columns'] = ('value')
        self.mat_properties.heading('#0', text='Property', anchor='w')
        self.mat_properties.heading('value', text='Value', anchor='w')
        self.mat_properties.bind('<ButtonRelease-1>', self.print_properties)

        self.mat_database = solcore.ParameterSystem().database

        self.valid_materials = sorted(self.mat_database.sections())
        self.valid_materials.remove('Immediate Calculables')
        self.valid_materials.remove('Final Calculables')

        for material in self.valid_materials:
            self.mat_list.insert(tk.END, material)

    def load_material_properties(self, *args):
        """

        :param args:
        :return:
        """

        if args[0] in self.valid_materials:
            # This option gets all the properties of a given material and returns it
            material = args[0]
            T = args[1]
            x = args[2] if len(args) == 3 else 1

            parameters = {}

            for option in self.mat_database.options(material):
                parameters[option] = solcore.get_parameter(material, option, T=T, x=x)

            print(parameters)
            return parameters

        elif len(self.mat_list.curselection()) != 1:
            return

        else:
            self.mat_properties.delete(*self.mat_properties.get_children())
            self.mat_selected = self.mat_list.curselection()[0]
            material = self.valid_materials[self.mat_selected]

            for option in self.mat_database.options(material):
                x = solcore.get_parameter(material, option, T=300)
                self.mat_properties.insert('', 'end', text=option, values=(x))

    def print_properties(self, *args):
        curItem = self.mat_properties.focus()
        print(self.mat_properties.item(curItem))

if __name__ == '__main__':

    valid_materials = sorted(solcore.ParameterSystem().database.sections())
    valid_materials.remove('Immediate Calculables')
    valid_materials.remove('Final Calculables')

    for mat in valid_materials:
        print(mat)
        for option in solcore.ParameterSystem().database.options(mat):
            x = solcore.ParameterSystem().database.get(mat, option)
            print('\t {}  \t{}'.format(option, x))