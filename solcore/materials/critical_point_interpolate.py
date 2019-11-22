"""
This module interpolates the n and k data of an alloy based on the known values at certain specific compositions.

The way it works is by interpolating in a smart way certain critical points (usually the Adachi critical points) and
then filling the gaps in between.
"""
import numpy as np
import os
from configparser import ConfigParser


def interpolate_critical_points(critical_points, target):
    fractions_available = np.array(sorted(critical_points.keys()))

    target_index = fractions_available.searchsorted(target)
    lower = fractions_available[target_index - 1]
    upper = fractions_available[target_index]

    fraction_along = (target - lower) / (upper - lower)
    lower_points = critical_points[lower]
    upper_points = critical_points[upper]

    result_critical_points = {k: lower_points[k] + fraction_along * (upper_points[k] - lower_points[k]) for k in
                              lower_points.keys()}
    return result_critical_points, (lower, upper), fraction_along


def load_data_from_directory(directory):
    result_config = ConfigParser()
    result_config.read(os.path.join(directory, "critical_points.txt"))
    sections = result_config.sections()

    fractions_config_file = [float(s) for s in sections]

    critical_points = {}
    fraction_dict = {}
    for f, s in zip(fractions_config_file, sections):
        options = sorted(result_config.options(s))

        critical_points[f] = {o: np.array([float(f) for f in result_config.get(s, o).split()]) for o in (options) if
                              o != "file"}
        datapath = os.path.join(directory, result_config.get(s, "file"))

        fraction_dict[f] = np.loadtxt(datapath, unpack=True)

    return fraction_dict, critical_points


def split_data_at_points(data, points):
    split_at = data[0].searchsorted(points[:, 0])

    split_data = (np.split(data, split_at, axis=1))
    return split_data


def transform(data, from_critical_points, to_critical_points):
    # sort the critical points according to "from" side order
    key_order, ordered_from_critical_points = zip(*sorted(from_critical_points.items(), key=lambda item: item[1][0]))
    moo, ordered_to_critical_points = zip(
        *sorted(to_critical_points.items(), key=lambda item: key_order.index(item[0])))
    ordered_to_critical_points = np.array(ordered_to_critical_points)
    ordered_from_critical_points = np.array(ordered_from_critical_points)
    split_data = split_data_at_points(data, ordered_from_critical_points)

    full_interval_bounds = np.concatenate(([0], ordered_from_critical_points[:, 0], [1e10]))
    full_interval_bounds_target = np.concatenate(([0], ordered_to_critical_points[:, 0], [1e10]))

    full_interval_bounds_y = np.concatenate(([0], ordered_from_critical_points[:, 1], [0]))
    full_interval_bounds_target_y = np.concatenate(([0], ordered_to_critical_points[:, 1], [0]))

    result_x, result_y = [], []

    for interval, data_subset in enumerate(split_data):
        interval_bounds = full_interval_bounds[interval], full_interval_bounds[interval + 1]
        interval_bounds_y = full_interval_bounds_y[interval], full_interval_bounds_y[interval + 1]

        interval_bounds_target = full_interval_bounds_target[interval], full_interval_bounds_target[interval + 1]
        interval_bounds_target_y = full_interval_bounds_target_y[interval], full_interval_bounds_target_y[interval + 1]

        if len(data_subset) == 0:
            continue
        x, y, *_ = data_subset  # heh heh heh, *_

        fraction_along = (x - interval_bounds[0]) / (interval_bounds[1] - interval_bounds[0])

        fraction_up = (y - interval_bounds_y[0]) / (interval_bounds_y[1] - interval_bounds_y[0])

        transformed_x = interval_bounds_target[0] + fraction_along * (
            interval_bounds_target[1] - interval_bounds_target[0])
        transformed_y = interval_bounds_target_y[0] + fraction_up * (
            interval_bounds_target_y[1] - interval_bounds_target_y[0])

        result_x.append(transformed_x)
        result_y.append(transformed_y)

    return np.array((np.concatenate(result_x), np.concatenate(result_y)))


def critical_point_interpolate(data, critical_points, target_fraction, grid):
    result_critical_points, (lower_fract, higher_fract), fraction_along = \
        interpolate_critical_points(critical_points, target_fraction)

    interpolate_from_lower = transform(
        data[lower_fract],
        from_critical_points=critical_points[lower_fract],
        to_critical_points=result_critical_points)
    interpolate_from_higher = transform(
        data[higher_fract],
        from_critical_points=critical_points[higher_fract],
        to_critical_points=result_critical_points)

    newgrid_lower = np.interp(grid, interpolate_from_lower[0], interpolate_from_lower[1])
    newgrid_higher = np.interp(grid, interpolate_from_higher[0], interpolate_from_higher[1])

    return grid, (newgrid_higher * fraction_along + newgrid_lower * (1 - fraction_along)), result_critical_points


if __name__ == "__main__":

    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import animation

    # Animate only works with Tk, so we tell matpotlib to use it
    matplotlib.use("TkAgg")

    data, critical_points = load_data_from_directory(os.path.join(os.path.split(__file__)[0], 'Data', 'AlGaAs nk'))

    fig = plt.figure()

    plt.yscale("log")
    nm = np.linspace(200, 1000, 1000)
    for k in data.keys():
        # print (data[k][0]*1e9, data[k][1])
        plt.plot(data[k][0] * 1e9, data[k][1], label=str(k), color="grey")

    resultx, resulty, crit = critical_point_interpolate(data, critical_points, 0, nm * 1e-9)

    lines, = plt.plot(resultx * 1e9, resulty, label="result", color="black", linewidth=2)

    plt.xlim(0, 1000)

    crit = np.array(list(crit.values())).transpose()
    print(crit, "lr")
    markers, = plt.plot(crit[0] * 1e9, crit[1], "+", label="result", color="red", linewidth=2, markersize=10)


    def animate(i):
        f = i / 100
        resultx, resulty, crit = critical_point_interpolate(data, critical_points, f, nm * 1e-9)
        # crit = crit.transpose()
        lines.set_ydata(resulty)
        # markers.set_ydata(crit[1])
        # markers.set_xdata(crit[0]*1e9)
        fig.canvas.get_tk_widget().update()  # process events
        return lines  # , markers


    anim = animation.FuncAnimation(fig, animate, init_func=None,
                                   frames=100, interval=2, blit=False)

    # anim.save('mpp.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    w = animation.FFMpegFileWriter(fps=10)
    w.fps = 10
    # anim.save('ellipse.mp4', writer=w, fps=10, bitrate=100, extra_args=['-vcodec', 'libx264'])

    plt.show()


    #
    # print (lines)
    # for f in np.linspace(0.01,1,100):
    #     resultx, resulty, crit=critical_point_interpolate(data, critical_points, f,nm)
    #     crit = crit.transpose()
    #     lines.set_ydata(resulty)
    #     markers.set_ydata(crit[1])
    #     markers.set_xdata(crit[0]*1e9)
    #     fig.canvas.get_tk_widget().update() # process events
    #     
    #     draw()
    # for f in np.linspace(1,0.01,100):
    #     resultx, resulty, crit=critical_point_interpolate(data, critical_points, f,nm)
    #     crit = crit.transpose()
    #     
    #     lines.set_ydata(resulty)
    #     markers.set_ydata(crit[1])
    #     markers.set_xdata(crit[0]*1e9)
    #     
    #     fig.canvas.get_tk_widget().update() # process events
    #     
    #     draw()
    #     
    # ioff()
    # show()
    # # legend()
    #
