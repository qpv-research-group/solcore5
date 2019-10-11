# from scipy import interpolate as I
import numpy as np
from scipy.interpolate import make_interp_spline, splev
from typing import Union


class interp1d(object):
    """ Interpolate a 1D function.

    See Also
    --------
    splrep, splev - spline interpolation based on FITPACK
    UnivariateSpline - a more recent wrapper of the FITPACK routines
    """

    def __init__(self, x: np.ndarray, y: np.ndarray, kind: str = 'linear',
                 axis: int = -1, copy: bool = True, bounds_error: bool = False,
                 fill_value: float = np.nan):
        """ Initialize a 1D linear interpolation class.

        Description
        -----------
        x and y are arrays of values used to approximate some function f:
            y = f(x)
        This class returns a function whose call method uses linear
        interpolation to find the value of new points.

        Parameters
        ----------
        x : array
            A 1D array of monotonically increasing real values.  x cannot
            include duplicate values (otherwise f is overspecified)
        y : array
            An N-D array of real values.  y's length along the interpolation
            axis must be equal to the length of x.
        kind : str or int
            Specifies the kind of interpolation as a string ('linear',
            'nearest', 'zero', 'slinear', 'quadratic, 'cubic') or as an integer
            specifying the order of the spline interpolator to use.
        axis : int
            Specifies the axis of y along which to interpolate. Interpolation
            defaults to the last axis of y.
        copy : bool
            If True, the class makes internal copies of x and y.
            If False, references to x and y are used.
            The default is to copy.
        bounds_error : bool
            If True, an error is thrown any time interpolation is attempted on
            a value outside of the range of x (where extrapolation is
            necessary).
            If False, out of bounds values are assigned fill_value.
            By default, an error is raised.
        fill_value : float
            If provided, then this value will be used to fill in for requested
            points outside of the data range.
            If not provided, then the default is NaN.
        """

        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value

        if kind in ['zero', 'slinear', 'quadratic', 'cubic']:
            order = {'nearest': 0, 'zero': 0, 'slinear': 1,
                     'quadratic': 2, 'cubic': 3}[kind]
            kind = 'spline'
        elif isinstance(kind, int):
            order = kind
            kind = 'spline'
        elif kind not in ('linear', 'nearest'):
            raise NotImplementedError("%s is unsupported: Use fitpack " \
                                      "routines for other types." % kind)
        x = np.array(x, copy=self.copy)
        y = np.array(y, copy=self.copy)

        if x.ndim != 1:
            raise ValueError("the x array must have exactly one dimension.")
        if y.ndim == 0:
            raise ValueError("the y array must have at least one dimension.")

        # Force-cast y to a floating-point type, if it's not yet one
        if not issubclass(y.dtype.type, np.inexact):
            y = y.astype(np.float_)

        # Normalize the axis to ensure that it is positive.
        self.axis = axis % len(y.shape)
        self._kind = kind

        if kind in ('linear', 'nearest'):
            # Make a "view" of the y array that is rotated to the interpolation
            # axis.
            axes = list(range(y.ndim))
            del axes[self.axis]
            axes.append(self.axis)
            oriented_y = y.transpose(axes)
            minval = 2
            len_y = oriented_y.shape[-1]
            if kind == 'linear':
                self._call = self._call_linear
            elif kind == 'nearest':
                self.x_bds = (x[1:] + x[:-1]) / 2.0
                self._call = self._call_nearest
        else:
            axes = list(range(y.ndim))
            del axes[self.axis]
            axes.insert(0, self.axis)
            oriented_y = y.transpose(axes)
            minval = order + 1
            len_y = oriented_y.shape[0]
            self._call = self._call_spline
            self._spline = make_interp_spline(x, oriented_y, k=order)

        len_x = len(x)
        if len_x != len_y:
            raise ValueError("x and y arrays must be equal in length along "
                             "interpolation axis.")
        if len_x < minval:
            raise ValueError("x and y arrays must have at "
                             "least %d entries" % minval)
        self.x = x
        self.y = oriented_y

    def _call_linear(self, x_new: np.ndarray) -> np.ndarray:

        # 2. Find where in the orignal data, the values to interpolate
        #    would be inserted.
        #    Note: If x_new[n] == x[m], then m is returned by searchsorted.
        x_new_indices = np.searchsorted(self.x, x_new)

        # 3. Clip x_new_indices so that they are within the range of
        #    self.x indices and at least 1.  Removes mis-interpolation
        #    of x_new[n] = x[0]
        x_new_indices = x_new_indices.clip(1, len(self.x) - 1).astype(int)

        # 4. Calculate the slope of regions that each x_new value falls in.
        lo = x_new_indices - 1
        hi = x_new_indices

        x_lo = self.x[lo]
        x_hi = self.x[hi]
        y_lo = self.y[..., lo]
        y_hi = self.y[..., hi]

        # Note that the following two expressions rely on the specifics of the
        # broadcasting semantics.
        slope = (y_hi - y_lo) / (x_hi - x_lo)

        # 5. Calculate the actual value for each entry in x_new.
        y_new = slope * (x_new - x_lo) + y_lo

        return y_new

    def _call_nearest(self, x_new: np.ndarray) -> np.ndarray:
        """ Find nearest neighbour interpolated y_new = f(x_new)."""

        # 2. Find where in the averaged data the values to interpolate
        #    would be inserted.
        #    Note: use side='left' (right) to searchsorted() to define the
        #    halfway point to be nearest to the left (right) neighbour
        x_new_indices = np.searchsorted(self.x_bds, x_new, side='left')

        # 3. Clip x_new_indices so that they are within the range of x indices.
        x_new_indices = x_new_indices.clip(0, len(self.x) - 1).astype(np.intp)

        # 4. Calculate the actual value for each entry in x_new.
        y_new = self.y[..., x_new_indices]

        return y_new

    def _call_spline(self, x_new: np.ndarray) -> np.ndarray:
        x_new = np.asarray(x_new)
        result = splev(x_new.ravel(), self._spline)
        return result.reshape(x_new.shape + result.shape[1:])

    def __call__(self, x_new: np.ndarray) -> np.ndarray:
        """Find interpolated y_new = f(x_new).

        Parameters
        ----------
        x_new : number or array
            New independent variable(s).

        Returns
        -------
        y_new : ndarray
            Interpolated value(s) corresponding to x_new.

        """

        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        x_new = np.asarray(x_new)
        out_of_bounds = self._check_bounds(x_new)

        y_new = self._call(x_new)

        # Rotate the values of y_new back so that they correspond to the
        # correct x_new values. For N-D x_new, take the last (for linear)
        # or first (for other splines) N axes
        # from y_new and insert them where self.axis was in the list of axes.
        nx = x_new.ndim
        ny = y_new.ndim

        # 6. Fill any values that were out of bounds with fill_value.
        # and
        # 7. Rotate the values back to their proper place.

        if nx == 0:
            # special case: x is a scalar
            if out_of_bounds:
                if ny == 0:
                    return np.asarray(self.fill_value)
                else:
                    y_new[...] = self.fill_value
            return np.asarray(y_new)
        elif self._kind in ('linear', 'nearest'):
            y_new[..., out_of_bounds] = self.fill_value
            axes = list(range(ny - nx))
            axes[self.axis:self.axis] = list(range(ny - nx, ny))
            return y_new.transpose(axes)
        else:
            y_new[out_of_bounds] = self.fill_value
            axes = list(range(nx, ny))
            axes[self.axis:self.axis] = list(range(nx))
            return y_new.transpose(axes)

    def _check_bounds(self, x_new: np.ndarray) -> np.ndarray:
        """Check the inputs for being in the bounds of the interpolated data.

        Parameters
        ----------
        x_new : array

        Returns
        -------
        out_of_bounds : bool array
            The mask on x_new of values that are out of the bounds.
        """

        # If self.bounds_error is True, we raise an error if any x_new values
        # fall outside the range of x.  Otherwise, we return an array indicating
        # which values are outside the boundary region.
        below_bounds = x_new < self.x[0]
        above_bounds = x_new > self.x[-1]

        # !! Could provide more information about which values are out of bounds
        if self.bounds_error and below_bounds.any():
            raise ValueError("A value in x_new is below the interpolation "
                             "range.")
        if self.bounds_error and above_bounds.any():
            raise ValueError("A value in x_new is above the interpolation "
                             "range.")

        # !! Should we emit a warning if some values are out of bounds?
        # !! matlab does not.
        out_of_bounds = np.logical_or(below_bounds, above_bounds)
        return out_of_bounds


class BilinearInterpolation(object):
    """docstring for BilinearInterpolation"""

    def __init__(self, x: Union[list, np.ndarray], y: Union[list, np.ndarray],
                 z: Union[list, np.ndarray], fill_value: float = 0.):
        super(BilinearInterpolation, self).__init__()
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        # Check the shapes of x,y and z
        x_s = self.x.shape
        y_s = self.y.shape
        z_s = self.z.shape
        assert x_s[0] >= 3, "Need at least 3 points in x."
        assert y_s[0] >= 3, "Need at least 3 points in y."
        assert x_s[0] == z_s[
            0], "The number of rows in z (which is %s), must be the same as number of elements in x (which is %s)" % (
        x_s[0], z_s[0])
        assert y_s[0] == z_s[
            1], "The number of columns in z (which is %s), must be the same as number of elements in y (which is %s)" % (
        y_s[0], z_s[1])

    def __call__(self, xvalue: float, yvalue: float) -> np.ndarray:
        """The intepolated value of the surface."""
        try:
            x_i = 0
            found = False
            for i in list(range(0, len(self.x) - 1)):
                # print self.x[i], "<=", xvalue, "<", self.x[i+1]
                if self.x[i] <= xvalue < self.x[i + 1]:
                    # print "Found, returning index", i
                    x_i = i
                    found = True
                    break

            if not found:
                # print "Not found, returning 0."
                return 0.
            x_v = (self.x[x_i], self.x[x_i + 1])
            x_i = (x_i, x_i + 1)

            y_i = 0
            found = False
            for i in list(range(0, len(self.y) - 1)):
                # print self.y[i], "<=", yvalue, "<", self.y[i+1]
                if self.y[i] <= yvalue < self.y[i + 1]:
                    # print "Found, returning index", i
                    y_i = i
                    found = True
                    break

            if not found:
                # print "Not found, returning 0."
                return 0.
            y_v = (self.y[y_i], self.y[y_i + 1])
            y_i = (y_i, y_i + 1)

            D = (x_v[1] - x_v[0]) * (y_v[1] - y_v[0])
            ta = self.z[x_i[0], y_i[0]]
            tb = (x_v[1] - xvalue)
            tc = (y_v[1] - yvalue)
            t1 = ta * tb * tc / D
            # print t1, "=", ta, "*", tb, "*", tc, "/", D

            ta = self.z[x_i[1], y_i[0]]
            tb = (xvalue - x_v[0])
            tc = (y_v[1] - yvalue)
            t2 = ta * tb * tc / D
            # print t1, "=", ta, "*", tb, "*", tc, "/", D

            ta = self.z[x_i[0], y_i[1]]
            tb = (x_v[1] - xvalue)
            tc = (yvalue - y_v[0])
            t3 = ta * tb * tc / D
            # print t1, "=", ta, "*", tb, "*", tc, "/", D

            ta = self.z[x_i[1], y_i[1]]
            tb = (xvalue - x_v[0])
            tc = (yvalue - y_v[0])
            t4 = ta * tb * tc / D
            # print t1, "=", ta, "*", tb, "*", tc, "/", D

            # print "Therefore,"
            # print t1, "+", t2, "+", t3, "+", t4

            # print self.z
            return t1 + t2 + t3 + t4

        except ValueError:
            #            import pdb; pdb.set_trace()
            # Is one of the inputs an array?
            if type(xvalue) == list or type(xvalue) == tuple or type(xvalue) == type(np.array([])):
                constructed_list = []
                for sub_xvalue in xvalue:
                    constructed_list.append(self(sub_xvalue, yvalue))
                return np.array(constructed_list)


if __name__ == "__main__":
    # Some wavelength points
    x = np.array([400., 600., 601, 1000.])

    # Some angle values
    y = np.array([0., np.pi / 4, np.pi / 2])

    # Some reflectivity values
    #             rads @ 400nm,    rads @ 600nm,        rads @ 601m,             rads @ 1000nm
    z = np.array([[0., 0., 0.], [1e-9, 1e-9, 1e-9], [1 - 1e-9, 1 - 1e-9, 1 - 1e-9], [1., 1., 1.]])
    # print z.shape

    z_i = BilinearInterpolation(x, y, z)
    print(z_i(610, 0.1))
