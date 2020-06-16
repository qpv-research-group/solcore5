from pytest import approx
import numpy as np


def test_prepare_solar_cell_basics(prepare_test_cell):

    from solcore.solar_cell_solver import prepare_solar_cell
    from solcore.state import State

    solar_cell, widths = prepare_test_cell

    options = State()
    options.position = 1e-8

    prepare_solar_cell(solar_cell, options=options)

    # check total width
    assert solar_cell.width * 1e9 == approx(np.sum(widths))

    # check individual width calculations

    indiv_width_int = np.array([layer_obj.width for layer_obj in solar_cell]) * 1e9

    indiv_width = np.array(
        [
            widths[0],
            widths[1],
            np.sum(widths[2:5]),
            np.sum(widths[5:7]),
            np.sum(widths[7:9]),
        ]
    )

    assert indiv_width_int == approx(indiv_width)

    # check offset calculations

    offsets_int = np.array([layer_obj.offset for layer_obj in solar_cell]) * 1e9
    offsets = np.insert(
        np.cumsum(
            np.array(
                [
                    widths[0],
                    widths[1],
                    np.sum(widths[2:5]),
                    np.sum(widths[5:7]),
                    np.sum(widths[7:9]),
                ]
            )
        )[0:-1],
        0,
        0,
    )

    assert offsets_int == approx(offsets)

    # check junction kind
    assert solar_cell[2].kind == "DA"
    assert solar_cell[4].kind == "PDD"


def test_prepare_solar_cell_position(prepare_test_cell):
    from solcore.solar_cell_solver import process_position
    from solcore.state import State
    from solcore.structure import Layer, TunnelJunction, Junction

    solar_cell, widths = prepare_test_cell

    options = State()

    offset = 0
    layer_widths = []
    for j, layer_object in enumerate(solar_cell):

        # Independent layers, for example in a AR coating
        if type(layer_object) is Layer:
            layer_widths.append(layer_object.width)

        # Each Tunnel junctions can also have some layers with a given thickness.
        elif type(layer_object) is TunnelJunction:
            junction_width = 0
            for i, layer in enumerate(layer_object):
                junction_width += layer.width
                layer_widths.append(layer.width)
            solar_cell[j].width = junction_width

        # For each junction, and layer within the junction, we get the layer width.
        elif type(layer_object) is Junction:
            junction_width = 0
            for i, layer in enumerate(layer_object):
                layer_widths.append(layer.width)
                junction_width += layer.width
            solar_cell[j].width = junction_width

        solar_cell[j].offset = offset
        offset += solar_cell[j].width

    solar_cell.width = offset

    pos = 1e-8
    options.position = pos
    process_position(solar_cell, options, layer_widths)
    assert options.position == approx(np.arange(0, solar_cell.width, pos))

    pos = [1e-8]
    options.position = pos
    process_position(solar_cell, options, layer_widths)
    assert options.position == approx(np.arange(0, solar_cell.width, pos[0]))

    pos = np.random.rand(len(solar_cell)) * 1e-9
    options.position = pos
    process_position(solar_cell, options, layer_widths)
    assert options.position == approx(
        np.hstack(
            [
                np.arange(
                    layer_object.offset,
                    layer_object.offset + layer_object.width,
                    pos[j],
                )
                for j, layer_object in enumerate(solar_cell)
            ]
        )
    )

    pos = np.random.rand(len(widths)) * 1e-9
    options.position = pos
    process_position(solar_cell, options, layer_widths)
    layer_offsets = np.insert(np.cumsum(widths), 0, 0) * 1e-9
    assert options.position == approx(
        np.hstack(
            [
                np.arange(layer_offsets[j], layer_offsets[j] + layer_width, pos[j])
                for j, layer_width in enumerate(widths * 1e-9)
            ]
        )
    )

    options.position = None
    process_position(solar_cell, options, layer_widths)
    layer_offsets = np.insert(np.cumsum(layer_widths), 0, 0)
    pos = [max(1e-10, width / 5000) for width in layer_widths]
    expected = np.hstack(
        [
            np.arange(layer_offsets[j], layer_offsets[j] + layer_width, pos[j])
            for j, layer_width in enumerate(layer_widths)
        ]
    )
    assert options.position == approx(expected)

    pos = np.arange(0, solar_cell.width, 1e-8)
    options.position = pos
    process_position(solar_cell, options, layer_widths)
    assert np.all(options.position == pos)
