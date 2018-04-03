import re
import numpy
from configparser import ConfigParser
import os
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt




interactive_terminal = False

float_re = re.compile(r"[-+]?[0-9]*\.?[0-9]+")

class CPPicker:
    """ This class can be used to pick the critical points out of a collection of n and k curves for a range of
    compositions in ternary alloys. Later, this critical points can be used later to calculate the n and k data at any
    other intermediate composition.

    >>> directory = '/InGaSb-Material/k.txt'
    >>> pickerApp = CPPicker(directory)
    >>> pickerApp.loop()

    """
    def __init__(self, directory):
        self.current_critical_point = 0
        self.current_curve = 0
        self.directory = directory
        self.gather_files()
        self.fraction_order = list(sorted(self.fractions_dict.keys()))
        self.n_curves = len(self.fraction_order)
        self.critical_points = [[]]
        self.finished = False
        self.log = True
        self.snap_to_data = True

    def guess_fraction(self, file_name):
        re_matches = float_re.findall(file_name)
        return re_matches[0] if len(re_matches) != 0 else None

    def gather_files(self):
        full_directory_path = self.directory
        candidate_files = [
            os.path.join(directory, f)
            for f in os.listdir(full_directory_path)
            if ".txt" in f and "critical" not in f
            ]
        results = {}
        for file_path in sorted(candidate_files):
            _, file_name = os.path.split(file_path)
            fraction_guess = self.guess_fraction(file_name)
            print("File '{}': material fraction guess = {} {}".format(
                file_name,
                fraction_guess,
                ("OK?\n> " if interactive_terminal else "\n")
            ), end="")

            answer = input() if interactive_terminal else ""
            fraction = fraction_guess if answer == "" else answer

            if interactive_terminal:
                print("using: {}".format(fraction))

            file_info = {
                "file_path": file_path,
                "file_name": file_name,
                "fraction": fraction,
                "data": numpy.loadtxt(file_path, unpack=True)
            }
            results[fraction] = file_info
        self.fractions_dict = results

    def plot_a_curve(self, curve_id, background=False):

        data = self.fractions_dict[curve_id]["data"]

        if (curve_id == self.fraction_order[self.current_curve] and not background):
            colour = "red"
            picker = 1000
            width = 3
        elif (self.current_curve != 0 and curve_id == self.fraction_order[self.current_curve - 1]):
            colour = "black"
            picker = None
            width = 1
        else:
            colour = "lightgrey"
            picker = None
            width = 1

        plt.plot(data[0] * 1e9, data[1], "-", color=colour, picker=picker, linewidth=width)

    def setup_plot(self, keep_limits=False):
        self.fig = plt.figure()
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        plt.ion()
        self.plot_all(keep_limits=False)

    def plot_all(self, keep_limits=True):
        if keep_limits:
            ranges = plt.xlim(), plt.ylim()
            self.fig.clear()
            plt.xlim(*ranges[0])
            plt.ylim(*ranges[1])
        else:
            self.fig.clear()
        if self.log:
            plt.yscale("log")
        else:
            plt.yscale("linear")

        for fraction in self.fractions_dict:
            self.plot_a_curve(fraction, background=True)
        self.plot_a_curve(self.fraction_order[self.current_curve], background=False)

        x_old = []
        y_old = []

        for critical_point_sublist in self.critical_points:
            x, y = [point[0] for point in critical_point_sublist], [point[1] for point in critical_point_sublist]
            plt.plot(x, y, ":+", color="blue", markersize=10, )

        plt.draw()
        plt.show()

    def loop(self):
        current_critical_point = 0
        self.setup_plot()

        while self.finished == False:
            self.fig.canvas.get_tk_widget().update()

        self.save()

    def onpick(self, event):

        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        x = event.mouseevent.xdata
        y = numpy.interp(x, xdata, ydata) if self.snap_to_data else event.mouseevent.ydata

        if self.current_curve < self.n_curves:
            if len(self.critical_points[-1]) < self.n_curves:
                print("critical for fraction {} at: x,y=({},{})".format(self.fraction_order[self.current_curve], x, y))
                self.critical_points[-1].append((x, y))

        if self.current_curve < self.n_curves - 1:
            self.current_curve += 1

        self.plot_all()

    def onkey(self, event):
        print("event.key", event.key)
        if event.key == "backspace":
            if len(self.critical_points[-1]) != 0:
                self.critical_points[-1].pop(-1)
                self.current_curve = len(self.critical_points[-1])
        if event.key == "enter":
            if len(self.critical_points[-1]) != self.n_curves:
                print("Insufficient critical points, continuing...")
            else:
                self.critical_points.append(list())
                self.current_curve = 0
        if event.key == "escape":
            self.finished = True
        if event.key == 'l':
            self.log = not self.log
        if event.key == 'n':
            self.snap_to_data = not self.snap_to_data

        self.plot_all()

    def save(self):

        out_path = os.path.join(self.directory, "temp_critical_points.txt")
        print("output path: {}".format(out_path))
        if interactive_terminal:
            print("chante_to? >")

        out_path = input() if interactive_terminal else out_path
        if interactive_terminal:
            print("using: {}".format(out_path))

        self.result_config = ConfigParser()

        for section in self.fraction_order:
            self.result_config.add_section(section)
            self.result_config.set(section, "file", self.fractions_dict[section]["file_name"])

        for critical_point_number, critical_points_subset in enumerate(self.critical_points[:-1]):
            for fraction_number, fraction in enumerate(self.fraction_order):
                c = critical_points_subset[fraction_number]
                self.result_config.set(fraction, str(critical_point_number + 1), "{} {}".format(c[0] * 1e-9, c[1]))

        with open(out_path, "w") as f:
            self.result_config.write(f)


if __name__ == '__main__':
    directory = '/Users/diego/OneDrive - Imperial College London/Developement/solcore5/solcore/material_data/AlInP-Material/k'
    pickerApp = CPPicker(directory)
    pickerApp.loop()
