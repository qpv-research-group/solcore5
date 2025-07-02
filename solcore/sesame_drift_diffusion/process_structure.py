from __future__ import annotations

from solsesame.builder import Builder
import numpy as np
from scipy.optimize import root
from solcore.constants import q, kb
from solcore.state import State
from solcore.structure import Junction
from solcore import material

def process_structure(junction: Junction, options: State):
    """
    Process the layers in the junction and convert it to a format which can be used by Sesame. This includes unit
    conversions from Solcore base SI units to the units used by Sesame (eV, cm). This function also sets the mesh,
    either from a user-supplied junction.mesh (in m) attribute, or by trying to guess an appropriate mesh depending on the
    doping profile and layer thicknesses. It will also extract the doping profile, which can be supplied in
    the following formats:

        - a constant value for each layer, set by the Nd or Na attribute of the material
        - a doping profile for the whole junction, passed as a function which accepts a depth in the junction in m
          and returns doping in m-3
        - a doping profile for each layer, passed as a function which accepts a depth in the layer in m and returns
          doping in m-3
        - Combinations of the above

    Note that since we can specify a doping profile for the whole junction, the Sesame solver can support a junction
    which is made up of a single layer, unlike the depletion approximation and Fortran-based PDD solvers.

    This function will add or update the following attributes in the Junction object:

        - mesh: the mesh used by the Sesame solver, in m
        - mesh_cm: the mesh used by the Sesame solver, in cm
        - sesame_sys: a Sesame Builder object defining the junction in Sesame format, used to solve the IV and QE

    :param junction: a Junction object
    :param options: a State object containing options for the solver

    """

    materials = []  # list of dictionaries for Sesame
    layer_widths = []  # layer widths in cm
    doping_profile_functions = (
        []
    )  # list of functions which take an argument of depth in junction in m and
    # return doping in cm-3
    offset = 0  # offset relative to front of junction of the current layer in cm

    for layer in junction:
        layer_widths.append(layer.width * 1e2)  # m to cm

        # this will always be set, even if the layer has a doping profile
        # constant_doping = layer.material.Nd / 1e6 if layer.material.Nd > layer.material.Na else -layer.material.Na /1e6

        if hasattr(layer, "doping_profile"):
            doping_profile_functions.append(
                lambda i, x: junction[i].doping_profile(x) / 1e6
            )
            # function which takes an argument of depth in junction in m and returns doping in cm-3

        else:
            doping_profile_functions.append(
                lambda i, x: junction[i].material.Nd / 1e6
                if junction[i].material.Nd > junction[i].material.Na
                else -junction[i].material.Na / 1e6
            )
            # to be consistent with functions passed by user, this should be a function which takes an argument in m
            # and returns doping in cm-3

        offset += layer.width * 1e2
        # Note: in Sesame, all of these can be functions of position as well as constants - but this is untested!

        new_mat = get_material_parameters(layer.material)
        materials.append(new_mat)

    # see if we can determine the junction depth:

    edges = np.insert(
        np.cumsum(layer_widths), 0, 0
    )  # locations of interfaces between layers
    edges[-1] = (
        edges[-1] + 1e-10
    )  # otherwise final point will not be assigned any values

    if "sesame_periodic" in options:
        periodic_bcs = options.sesame_periodic

    else:
        periodic_bcs = True

    if hasattr(junction, "doping_profile"):
        junction_depth = (
            100 * root(junction.doping_profile, np.sum(layer_widths) / (2 * 100)).x
        )  # in cm
        # this is where doping crosses zero - not necessarily the junction depth, but not sure how best to determine
        # this. Maybe want to find slope of doping profile and identify locations where it is steepest?

    else:
        # identify where the doping profile changes sign for constant doping layers or individual doping profiles
        # per layer
        doping_sign = [
            np.sign(doping_profile_functions[i1](i1, edges[i1] / 100))
            for i1 in range(len(edges) - 1)
        ]
        junction_depth = edges[np.where(np.diff(doping_sign))[0] + 1]

    if not hasattr(junction, "mesh"):
        # calculate appropriate mesh (hopefully) if not supplied as input in junction.mesh
        make_mesh(junction, layer_widths, options, junction_depth)
        user_mesh = False

    else:
        # user-provided mesh
        junction.mesh_cm = junction.mesh * 1e2  # m to cm
        user_mesh = True

    if hasattr(junction, "doping_profile"):
        doping_profile_x = junction.doping_profile(junction.mesh) / 1e6  # in cm-3

    else:
        doping_profile_x = np.zeros(len(junction.mesh_cm))
        for i1, doping_profile_fn in enumerate(doping_profile_functions):
            doping_profile_x[
                np.all(
                    (junction.mesh_cm >= edges[i1], junction.mesh_cm < edges[i1 + 1]),
                    axis=0,
                )
            ] = doping_profile_fn(
                i1,
                junction.mesh[
                    np.all(
                        (
                            junction.mesh_cm >= edges[i1],
                            junction.mesh_cm < edges[i1 + 1],
                        ),
                        axis=0,
                    )
                ]
                - edges[i1] / 100,
            )

    # Make Sesame Builder object for simulations
    junction.sesame_sys = Builder(junction.mesh_cm, T=options.T, periodic=periodic_bcs)  # Sesame system

    junction.sesame_sys.rho = doping_profile_x / junction.sesame_sys.scaling.density

    if (
        not user_mesh
    ):  # do not update mesh after finding doping profile if a specific mesh was supplied by the user
        update_mesh(junction, layer_widths, options)

        if hasattr(junction, "doping_profile"):
            doping_profile_x = junction.doping_profile(junction.mesh) / 1e6  # in cm-3

        else:
            doping_profile_x = np.zeros(len(junction.mesh_cm))
            for i1, doping_profile_fn in enumerate(doping_profile_functions):
                doping_profile_x[
                    np.all(
                        (
                            junction.mesh_cm >= edges[i1],
                            junction.mesh_cm < edges[i1 + 1],
                        ),
                        axis=0,
                    )
                ] = doping_profile_fn(
                    i1,
                    junction.mesh[
                        np.all(
                            (
                                junction.mesh_cm >= edges[i1],
                                junction.mesh_cm < edges[i1 + 1],
                            ),
                            axis=0,
                        )
                    ]
                    - edges[i1] / 100,
                )

        junction.sesame_sys = Builder(junction.mesh_cm, T=options.T, periodic=periodic_bcs)  # Sesame system

        # set doping profile in Sesame Builder objects
        junction.sesame_sys.rho = doping_profile_x / junction.sesame_sys.scaling.density

    for i1 in range(len(junction)):
        junction.sesame_sys.add_material(
            materials[i1], lambda x: np.all((x >= edges[i1], x < edges[i1 + 1]), axis=0)
        )

    junction.sesame_sys.contact_type("Ohmic", "Ohmic")
    # should also be able to choose other options, since Sesame can handle them

    # get surface recombination velocities
    junction.sesame_sys.contact_S(*get_srv(junction))

    # make sure that the mesh points are all inside the range covered by
    # options.position, otherwise will get NaNs:



def get_material_parameters(mat: material):
    """
    Extract approriate material parameters for Sesame from a Solcore material object.

    :param mat: a Solcore material
    """
    try:
        bulk_recombination_energy = mat.bulk_recombination_energy / q
        # this is an option which can be used by Sesame, but is not currently used by Solcore
        # anywhere else, or defined in the database for any materials
    except ValueError:
        bulk_recombination_energy = 0  # default value

    try:
        auger_electron = (
            mat.electron_auger_recombination * 1e12
        )  # m6/s to cm6/s

    except ValueError:
        auger_electron = 0

    try:
        auger_hole = mat.hole_auger_recombination * 1e12  # m6/s to cm6/s

    except ValueError:
        auger_hole = 0

    # mobility must be provided, but if the user has provided a diffusion length, we can calculate
    # the minority carrier lifetime from the mobility and diffusion length
    try:
        electron_minority_lifetime = mat.electron_minority_lifetime
    except ValueError:
        electron_minority_lifetime = carrier_constants("electron_minority_lifetime", mat)

    try:
        hole_minority_lifetime = mat.hole_minority_lifetime
    except ValueError:
        hole_minority_lifetime = carrier_constants("hole_minority_lifetime", mat)

    new_mat = {
        "Nc": mat.Nc
        * 1e-6,  # effective density of states at CB edge (cm-3)
        "Nv": mat.Nv
        * 1e-6,  # effective density of states at VB edge (cm-3)
        "Eg": mat.band_gap / q,  # material bandgap (eV)
        "affinity": mat.electron_affinity / q,  # electron affinity (eV)
        "epsilon": mat.relative_permittivity,
        # relative permittivity (dimensionless)
        "mu_e": mat.electron_mobility
        * 1e4,  # electron mobility (cm2/(V s))
        "mu_h": mat.hole_mobility * 1e4,  # hole mobility (cm2/(V s))
        "tau_e": electron_minority_lifetime,  # electron bulk lifetime (s)
        "tau_h": hole_minority_lifetime,  # hole bulk lifetime (s)
        "Et": bulk_recombination_energy,  # energy level of bulk recombination centres (eV)
        "B": mat.radiative_recombination * 1e6,
        # radiative recombination constant (cm3/s)
        "Cn": auger_electron,  # Auger recombination constant for electrons (cm6/s)
        "Cp": auger_hole,  # Auger recombination constant for holes (cm6/s)
    }

    return new_mat

def carrier_constants(quantity, mat):
    if quantity == "electron_diffusion_length":
        mobility, lifetime = mat.electron_mobility, mat.electron_minority_lifetime
        return np.sqrt(kb * mat.T * mobility * lifetime / q)

    if quantity == "hole_diffusion_length":
        mobility, lifetime = mat.hole_mobility, mat.hole_minority_lifetime
        return np.sqrt(kb * mat.T * mobility * lifetime / q)

    if quantity == "electron_minority_lifetime":
        mobility, diffusion_length = mat.electron_mobility, mat.electron_diffusion_length
        return q * diffusion_length ** 2 / (kb * mat.T * mobility)

    if quantity == "hole_minority_lifetime":
        mobility, diffusion_length = mat.hole_mobility, mat.hole_diffusion_length
        return q * diffusion_length ** 2 / (kb * mat.T * mobility)


def make_mesh(
    junction: Junction, layer_width: float, options: State, junction_depth: list[float]
):
    """
    Make a mesh for Sesame, in case one was not supplied. This will try to make the mesh dense near
    the front surface, at the (guessed) junction depth, and near the rear surface. The user can specify
    the minimum and maximum spacing allowed (in m) through options.minimum_spacing or options.maximum_spacing.
    Otherwise these will be set so that there is a maximum of 10,000 and a minimum of 500 points in the
    mesh

    :param junction: a Junction object
    :param layer_width: a list of layer widths in m
    :param options: a State object containing options for the solver
    :param junction_depth: a list of junction depths in cm
    """

    # TODO: general improvement of the meshing; should generate a smooth mesh rather than either coarse or fine

    total_width = np.sum(layer_width)

    if "minimum_spacing" in options.keys():
        minimum_spacing = options.minimum_spacing * 1e2

    else:
        minimum_spacing = total_width / 1e4  # maximum of 10000 points

    if "maximum_spacing" in options.keys():
        maximum_spacing = options.maximum_spacing * 1e2

    else:
        maximum_spacing = total_width / 500  # minimum of 500 points

    if maximum_spacing == minimum_spacing:
        mesh = np.arange(0, total_width, maximum_spacing)

    # identify location of junction. doping_profile is a function

    else:
        # edges = np.insert(np.cumsum(layer_width), 0, 0)

        front_spacing = np.max([minimum_spacing, 0.5e-7])

        dense_front_width = 200e-7

        end_front = np.min([dense_front_width - 1e-10, total_width])

        front = np.arange(0, end_front, front_spacing)

        # 1 nm in cm = 1e-7
        if end_front != total_width:
            if junction_depth[0] < total_width / 2:
                end_dense = (
                    junction_depth[0] + 0.1 * total_width
                    if (junction_depth[0] + 0.1 * total_width < total_width)
                    else total_width
                )

                include_end = True if end_dense == total_width else False

                if end_dense > end_front:
                    n_points = int((end_dense - end_front) / minimum_spacing)
                    mesh = np.linspace(
                        end_front, end_dense, n_points, endpoint=include_end
                    )

                    mesh = np.concatenate((front, mesh))

                else:
                    mesh = front

                if not include_end:
                    n_points = int((total_width - end_dense) / maximum_spacing)
                    mesh_rest = np.linspace(
                        end_dense, total_width, n_points, endpoint=True
                    )

                    mesh = np.concatenate((mesh, mesh_rest))

            else:
                n_points = int((0.2 * total_width - end_front) / minimum_spacing)

                mesh_dense_front = np.linspace(
                    end_front, 0.2 * total_width, n_points, endpoint=False
                )

                n_points = int(
                    (junction_depth[0] - 0.3 * total_width) / maximum_spacing
                )

                mesh_bulk = np.linspace(
                    0.2 * total_width,
                    junction_depth[0] - 0.1 * total_width,
                    n_points,
                    endpoint=False,
                )

                end_dense = (
                    junction_depth[0] + 0.1 * total_width
                    if (junction_depth[0] + 0.1 * total_width < total_width)
                    else total_width
                )

                include_end = True if end_dense == total_width else False

                n_points = int(
                    (end_dense - (junction_depth[0] - 0.1 * total_width))
                    / minimum_spacing
                )

                mesh_dense = np.linspace(
                    junction_depth[0] - 0.1 * total_width,
                    end_dense,
                    n_points,
                    endpoint=include_end,
                )

                if not include_end:
                    n_points = int((total_width - end_dense) / maximum_spacing)

                    mesh_rest = np.linspace(
                        end_dense, total_width, n_points, endpoint=True
                    )

                    mesh = np.unique(
                        np.concatenate(
                            (front, mesh_dense_front, mesh_bulk, mesh_dense, mesh_rest)
                        )
                    )

                else:
                    mesh = np.unique(
                        np.concatenate((front, mesh_dense_front, mesh_bulk, mesh_dense))
                    )

    junction.mesh_cm = mesh
    junction.mesh = junction.mesh_cm * 1e-2  # cm to m


def update_mesh(junction: Junction, layer_width: list[float], options: State):
    """
    Update mesh after the doping profile is known.

    :param junction: a Junction object
    :param layer_width: a list of layer widths in cm
    :param options: a State object containing options for the solver
    """

    total_width = np.sum(layer_width)

    # always do 1 nm spacing near front surface?

    if "minimum_spacing" in options.keys():
        minimum_spacing = options.minimum_spacing * 1e2

    else:
        minimum_spacing = total_width / 1e4  # maximum of 10000 points

    current_mesh = junction.mesh_cm

    diff_rho = np.abs(np.gradient(junction.sesame_sys.rho, current_mesh))

    doping_change = np.where(diff_rho / np.max(diff_rho) > 1e-4)[0]

    # identify continuous stretches of high doping change, and where there are gaps:

    change_point = np.where(np.diff(doping_change) != 1)[0]

    last_points = doping_change[change_point]
    first_points = doping_change[change_point + 1]

    first_points = np.insert(first_points, 0, doping_change[0])
    last_points = np.append(last_points, doping_change[-1])

    first_x = current_mesh[first_points]
    last_x = current_mesh[last_points]

    for i1 in range(len(first_points)):
        # want to make mesh denser where there are large changes in doping
        # print('denser mesh between', first_x[i1], last_x[i1])

        mesh_above = current_mesh[current_mesh < first_x[i1]]
        mesh_below = current_mesh[current_mesh >= last_x[i1]]
        mesh_dense = np.arange(
            first_x[i1], last_x[i1] - 1e-12, minimum_spacing
        )  # last point not included

        current_mesh = np.unique(np.concatenate((mesh_above, mesh_dense, mesh_below)))

    junction.mesh_cm = current_mesh
    junction.mesh = junction.mesh_cm * 1e-2  # cm to m


def get_srv(junction: Junction):
    """
    Extract approriate surface recombination velocities for Sesame from a Junction object.
    The labelling sn and sp in the junction object refers to the doping of the region, NOT the
    carrier type. For the Sesame solver, it is also possible to specify the recombination velocities
    for the carrier types separately in the defintion of the Junction as sn_e/sn_h, sp_e/sp_h,
    being respectively the electron/hole SRVs in the n-type region and the electron/hole SRVs in the
    p_type region. Note that the majority carrier surface recombination (sn_e and sp_h) generally
    do not affect the solution.

    :param junction: a Junction object
    """

    if junction.sesame_sys.rho[junction.sesame_sys.nx - 1] < 0:
        # p-type at back
        polarity = "np"

    else:
        polarity = "pn"

    if hasattr(junction, "sn_e"):
        Sn_e = junction.sn_e * 100

    elif hasattr(junction, "sn"):
        Sn_e = junction.sn * 100

    else:
        Sn_e = 1e4  # same as Fortran PDD default (this one is in cm s-1, not m s-1

    if hasattr(junction, "sp_e"):
        Sp_e = junction.sp_e * 100

    elif hasattr(junction, "sp"):
        Sp_e = junction.sp * 100

    else:
        Sp_e = 1e4

    if hasattr(junction, "sn_h"):
        Sn_h = junction.sn_h * 100

    elif hasattr(junction, "sn"):
        Sn_h = junction.sn * 100

    else:
        Sn_h = 1e4  # is this the same as the Fortran PDD?

    if hasattr(junction, "sp_h"):
        Sp_h = junction.sp_h * 100

    elif hasattr(junction, "sp"):
        Sp_h = junction.sp * 100

    else:
        Sp_h = 1e4

    if polarity == "pn":
        return Sp_e, Sp_h, Sn_e, Sn_h

    else:
        return Sn_e, Sn_h, Sp_e, Sp_h

