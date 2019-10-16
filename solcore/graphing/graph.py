import tempfile
from copy import copy

from matplotlib import font_manager
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.colors import Normalize  # ColorConverter
import matplotlib.cm as cmx
from typing import List, Tuple, Dict, Union, Any
from .graph_support import flatten, open_with_os

graph_defaults: Dict[str, Any] = {
    "format": "pdf",
    "linetype": "-",
    "dpi": 100,
    "legendframe": True,
    "font_family": "serif",
    "palette": "bright hue wheel",
    "square": True,
    "colormap": "nipy_spectral",
    "grid": {"dash_capstyle": "round", "zorder": 1, "alpha": 0.5}
}

mpl_simple_properties: Tuple[str, ...] = (
"title", "xlabel", "ylabel", "zlabel", "xlim", "ylim", "zlim3d", "xticks", "xticklabels", "yticks", "yticklabels",
"yscale", "xscale",)
directions: Tuple[str, ...] = ("left", "right", "top", "bottom")


def make_square(fig: Figure, ax: Axes) -> None:
    bb: Any = ax.get_position()
    fwidth: float = fig.get_figwidth()
    fheight: float = fig.get_figheight()
    axwidth: float = fwidth * (bb.x1 - bb.x0)
    axheight: float = fheight * (bb.y1 - bb.y0)

    # square_edge = min((axwidth, axheight))
    if axwidth > axheight:
        narrow_by: float = (axwidth - axheight) / fwidth
        bb.x0 += narrow_by / 2
        bb.x1 -= narrow_by / 2
    elif axheight > axwidth:
        shrink_by: float = (axheight - axwidth) / fheight
        bb.y0 += shrink_by / 2
        bb.y1 -= shrink_by / 2
    ax.set_position(bb)


class Graph:
    suppress_show: bool = False
    plotted = 0

    def __init__(self, *data: Any, **options: Any) -> None:
        self.axis: Any = None
        self.options = copy(graph_defaults)
        self.options.update(options)
        self.data = list(flatten(data))
        self.extra_artists: List = []
        self.subplots: List = []
        self.create_figure()

    def main_drawing_flow(self, figure: Figure) -> None:
        # print ("mdf")
        for subplot in self.subplots:
            subplot.main_drawing_flow(figure)
        self.create_axis(figure)
        self.setup_matplotlib()
        self.render_data()
        # if "yscale" in self.options and "xscale" in self.options and self.options["xscale"] == self.options["yscale"]:
        #     self.axis.set_aspect(1./self.axis.get_data_ratio())
        # else:
        #     print ("cannot guarantee square log-lin plots yet")

        if "legend" in self.options:
            self.draw_legend()
            self.extra_artists.append(self.legend)

        if "text" in self.options:
            self.extra_artists.append(self.text)

        if "square" in self.options and self.options["square"]:
            make_square(self.figure, self.axis)

    def add_subplot(self, *data: Any, **options: Any) -> None:
        self.subplots.append(Graph(*data, figure=self.figure, **options))

    def create_figure(self):
        if "figure" in self.options:
            self.figure = self.options["figure"]
        elif "figsize" in self.options:
            self.figure = Figure(figsize=self.options["figsize"])
        else:
            self.figure = Figure()
        self.figure.subplots_adjust(**{k: self.options[k] for k in self.options.keys() if k in directions})

    def create_axis(self, figure: Figure) -> None:
        if "axis" in self.options:
            self.axis = self.options["axis"]
        elif self.axis is not None:
            # print ("already has axis")
            return
        else:
            subplot_args = self.options["subplot_args"] if "subplot_args" in self.options else {}
            self.axis = figure.add_subplot(
                self.options["subplot"] if "subplot" in self.options else 111,
                **subplot_args
            )

    def setup_matplotlib(self) -> None:
        simple_setters = {k: self.options[k] for k in mpl_simple_properties if k in self.options}
        for key, value in simple_setters.items():
            getattr(self.axis, "set_{}".format(key))(value)

        if "grid" in self.options and self.options["grid"] is not None:
            self.axis.grid(**self.options["grid"])

        if "xhide" in self.options: self.axis.get_xaxis().set_visible(False)
        if "yhide" in self.options: self.axis.get_yaxis().set_visible(False)
        if "xlabelcolor" in self.options: self.axis.get_xaxis().label.set_color(self.options["xlabelcolor"])
        if "ylabelcolor" in self.options: self.axis.get_yaxis().label.set_color(self.options["ylabelcolor"])
        if "x_major_formatter" in self.options: self.axis.get_xaxis().set_major_formatter(
            self.options["x_major_formatter"])
        if "y_major_formatter" in self.options: self.axis.get_yaxis().set_major_formatter(
            self.options["y_major_formatter"])

        # if "xticks" in self.options:
        #     x_coordinate, x_label= list(zip(*self.options["xticks"]))
        #     self.axis.set_xticks(x_coordinate)
        #     self.axis.set_xticklabels(x_label)

        if "yticklabelcolor" in self.options and self.options["yticklabelcolor"] is not None:
            for tl in self.axis.get_yticklabels():
                tl.set_color(self.options["yticklabelcolor"])
        if "xticklabelcolor" in self.options and self.options["xticklabelcolor"] is not None:
            for tl in self.axis.get_xticklabels():
                tl.set_color(self.options["xticklabelcolor"])

    def render_data(self) -> None:

        all_color_indexes = [d.color_index if (
            hasattr(d, "color_index")
            and d.color_index is not None
            and not ("color" in d.line_format)
        ) else -i - 1 for i, d in enumerate(self.data) if not (
            "color" in d.line_format or
            "facecolor" in d.line_format or
            "edgecolor" in d.line_format)]

        unique_color_indexes = len(set(all_color_indexes))

        cNorm = Normalize(vmin=0, vmax=unique_color_indexes)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=self.options["colormap"])

        self.legend_additions: List = []

        for i, datum in enumerate(self.data):
            if "color" not in datum.line_format and "edgecolor" not in datum.line_format and "facecolor" not in datum.line_format:
                if hasattr(datum, "color_index") and datum.color_index is not None:
                    color = scalarMap.to_rgba(datum.color_index)
                else:
                    color = scalarMap.to_rgba(i)

                datum.line_format["color"] = color

            if "edit" in self.options:
                datum.edit = self.options["edit"]
            legend_entry = datum.draw(self.axis)
            if legend_entry is not None:
                self.legend_additions.append((i, legend_entry[0], legend_entry[1])
                                             )

    def draw_legend(self) -> None:
        font_args = {}

        if "legend_font_size" in self.options:
            font_args["size"] = self.options["legend_font_size"]
        if "legend_font_name" in self.options:
            font_args["fontname"] = self.options["legend_font_name"]
        if "legend_font_weight" in self.options:
            font_args["weight"] = self.options["legend_font_weight"]
        if "legend_font_color" in self.options:
            font_args["color"] = self.options["legend_font_color"]

        self.legend_font_properties = font_manager.FontProperties(**font_args)

        legend_location = self.options["legend"]
        handles, labels = self.axis.get_legend_handles_labels()
        self.legend_additions.sort()

        for i, (position, handle, label) in enumerate(self.legend_additions):
            handles.insert(position, handle)
            labels.insert(position, label)

        final_handles, final_labels = [], []
        for handle, label in zip(handles, labels):
            if label in ("None", None):
                continue
            final_labels.append(label)
            final_handles.append(handle)

        legend_args = (final_handles, final_labels)

        legend_kwargs = {
            "frameon": self.options["legendframe"] if "legendframe" in self.options else False,
            "title": self.options["legendtitle"] if "legendtitle" in self.options else False,
            "ncol": self.options["legendcolumns"] if "legendcolumns" in self.options else 1,
            "numpoints": 4,
            "prop": self.legend_font_properties,
            "handlelength": self.options["handlelength"] if "handlelength" in self.options else 2
        }
        if legend_location == "outside":
            legend_kwargs['bbox_to_anchor'] = self.options["legendbox"] if "legendbox" in self.options else (1.6, 1)
        else:
            legend_kwargs["loc"] = legend_location

        self.legend = self.axis.legend(*legend_args, **legend_kwargs)

    def draw(self, path: str = None, show: bool = True) -> None:
        self.main_drawing_flow(self.figure)

        if path is None:
            handle, path = tempfile.mkstemp(prefix="tmp_solcore_", suffix=".%s" % self.options["format"])

        if self.options["format"].upper() == "PDF":
            self.canvas = FigureCanvasPdf(self.figure)
        else:
            self.canvas = FigureCanvasAgg(self.figure)

        self.canvas.print_figure(path, dpi=self.options["dpi"], bbox_extra_artists=self.extra_artists,
                                 bbox_inches='tight')

        if show and not self.suppress_show:
            open_with_os(path)

    def show(self) -> None:
        import pylab
        self.main_drawing_flow(pylab.gcf())
        pylab.show()


class TwinnedAxisGraph(Graph):
    def __init__(self, *args: Any, twin:str = "x", **kwargs: Any) -> None:
        self.twin = twin
        self.line_format: Any = {}
        Graph.__init__(self, *args, **kwargs)

    def draw(self, axis: Axes) -> None:
        # print ("Td")
        if self.twin == "x":
            self.axis = axis.twinx()
        elif self.twin == "y":
            self.axis = axis.twiny()
        else:
            raise ValueError("Incorrect value for twinning_axis, must be x or y")
        self.main_drawing_flow(self.figure)
