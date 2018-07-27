from copy import copy
import numpy

from matplotlib.patches import Rectangle
from matplotlib import font_manager

line_defaults = {
    "linewidth": 0.5,
}


class GraphData_Base(object):
    def __init__(self, *data, color_index=None, edit=None, **options):
        if not (hasattr(self, "suppress_defaults") and self.suppress_defaults):
            self.line_format = copy(line_defaults)
        else:
            self.line_format = {}

        self.color_index = color_index
        if len(data) == 1:
            self.data = numpy.array(data[0])
        else:
            self.data = numpy.array(data)
        self.line_format.update(options)
        if edit is not None:
            self.edit = edit

    def draw(self, axis, figure=None):
        raise NotImplementedError("Subclasses must define the draw method")

    @property
    def data(self):
        if hasattr(self, "edit"):
            if len(self._data) == 2:
                return self.edit(*self._data)
            if len(self._data) == 3:
                x, y1 = self.edit(self._data[0], self._data[1])
                x, y2 = self.edit(self._data[0], self._data[2])
                return numpy.array((x, y1, y2))
        return self._data

    @data.setter
    def data(self, new_value):
        self._data = new_value
        # def legend_entry(self):
        #     return self.line_format["legend"]


class GraphDataHRect(GraphData_Base):
    def draw(self, axis, figure=None):
        x1, x2 = self.data
        axis.axhspan(x1, x2, **self.line_format)


class GraphDataVRect(GraphData_Base):
    def draw(self, axis, figure=None):
        x1, x2 = self.data
        axis.axvspan(x1, x2, **self.line_format)
        # if "label" in self.line_format:
        #     patch_color = {}
        #     patch_color.update(self.line_format)
        #     p = Rectangle((0, 0), 1, 1, **patch_color)
        #     return p, self.line_format["label"]


class GraphData(GraphData_Base):
    def draw(self, axis, figure=None):
        x, y = self.data
        axis.plot(x, y, **self.line_format)


class GraphDataErrorbars(GraphData_Base):
    def draw(self, axis, figure=None):
        assert (
        "xerr" in self.line_format or "yerr" in self.line_format), "GraphDataErrorbars requres xerr, yerr or both to be set"
        x, y = self.data
        axis.errorbar(x, y, **self.line_format)


class GraphDataBetween(GraphData_Base):
    def draw(self, axis, figure=None):
        x, y1, y2 = self.data
        axis.fill_between(x, y1, y2, **self.line_format)
        if "label" in self.line_format:
            patch_color = {}
            patch_color.update(self.line_format)
            p = Rectangle((0, 0), 1, 1, **patch_color)
            return p, self.line_format["label"]


class LegendEntry(GraphData_Base):
    def draw(self, axis, figure=None):
        p = Rectangle((0, 0), 1, 1, edgecolor='none', facecolor='none')

        return p, self.data


class GraphDataLabels(GraphData_Base):
    def __init__(self, *data, color_index=None, font_size=None, fontname=None, weight=None, font_color=None, **options):
        font_args = {}

        if font_size is not None:
            font_args["size"] = font_size
        if fontname is not None:
            font_args["font_name"] = fontname
        if weight is not None:
            font_args["weight"] = weight
        if font_color is not None:
            font_args["color"] = font_color

        self.font_properties = font_manager.FontProperties(**font_args)
        GraphData_Base.__init__(self, *data, color_index=color_index, **options)

    def draw(self, axis, figure=None):
        x, y, labels = self.data
        axis.scatter(x, y, **self.line_format)
        for X, Y, L in zip(x, y, labels):
            axis.annotate(" " + L, xy=(X, Y), fontproperties=self.font_properties)
            # if "label" in self.line_format:
            #     patch_color = {}
            #     patch_color.update(self.line_format)
            #     p = Rectangle((0, 0), 1, 1, **patch_color)
            #     return p, self.line_format["label"]


class Text(GraphData_Base):
    def draw(self, axis, figure=None, ):
        print(self._data[0], self._data[1], self._data[2])
        figure.text(self._data[0], self._data[1], self._data[2])


class GraphDataImage(GraphData_Base):
    suppress_defaults = True

    def draw(self, axis, figure=None):
        z = self.data
        image = axis.imshow(z, **{k: self.line_format[k] for k in self.line_format.keys() if k not in ["color", ]})
