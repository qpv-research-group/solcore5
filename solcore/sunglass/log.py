import tkinter as tk
from tkinter import ttk

class LogTab(ttk.Frame):
    """ This tab is used to print the output log of any running process
    """
    def __init__(self, parent):
        """ Constructor of the class

        :param parent: The notebook that serves as parent of this tab
        """

        super(LogTab, self).__init__(parent)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.text = tk.Text(self, wrap="word")
        scroll = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.text.yview)
        self.text.grid(column=0, row=0, sticky=(tk.NSEW))
        self.text.configure(yscrollcommand=scroll.set)
        scroll.grid(column=99, row=0, sticky=(tk.NS))

        self.text.tag_configure('red', foreground='red')
