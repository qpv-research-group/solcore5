from pytest import raises, approx


def test_validate():
    import xarray as xr
    import numpy as np
    from solcore.future.mesh import InvalidMesh, Mesh

    p = [1, 2, 3]

    # points not a dimension name
    z = xr.DataArray(
        p,
        dims=["y"],
        coords={
            "interval": ("y", np.ones_like(p).astype(int)),
            "group": ("y", np.ones_like(p).astype(int)),
        },
    )
    with raises(InvalidMesh):
        Mesh._validate(z)

    # No interval information
    z = xr.DataArray(
        p, dims=["points"], coords={"group": ("points", np.ones_like(p).astype(int))}
    )
    with raises(InvalidMesh):
        Mesh._validate(z)

    # No group information
    z = xr.DataArray(
        p,
        dims=["points"],
        coords={
            "interval": ("points", np.ones_like(p).astype(int)),
        },
    )
    with raises(InvalidMesh):
        Mesh._validate(z)

    # All OK
    z = xr.DataArray(
        p,
        dims=["points"],
        coords={
            "interval": ("points", np.ones_like(p).astype(int)),
            "group": ("points", np.ones_like(p).astype(int)),
        },
    )
    Mesh._validate(z)


def test_init_mesh():
    from solcore.future.mesh import Mesh
    import xarray as xr
    import numpy as np

    p = [1, 2, 3]
    z = xr.DataArray(
        p,
        dims=["points"],
        coords={
            "interval": ("points", np.ones_like(p).astype(int)),
            "group": ("points", np.ones_like(p).astype(int)),
        },
    )
    with raises(ValueError):
        Mesh(z, [0, 10])

    mesh = Mesh(z)
    assert all(mesh.nodes == [z[0].item(), z[-1].item()])


def test_uniform():
    from solcore.future.mesh import uniform

    nodes = [4, 40, 400]
    mesh = uniform(42, nodes)

    assert min(nodes) == min(mesh.z).item()
    assert max(nodes) == max(mesh.z).item()


def test_piecewise_uniform():
    from solcore.future.mesh import piecewise_uniform

    nodes = [4, 40, 400]

    # Insufficient number of points -> need to be equal to number of nodes minus 1
    with raises(ValueError):
        piecewise_uniform([42], nodes)

    npoints = [42, 100]
    mesh = piecewise_uniform(npoints, nodes)
    assert all(n in mesh.z for n in nodes)
    assert len(mesh.z) == sum(npoints)
    assert mesh.ninterval == 2
    assert mesh.ngroup == 1
    assert sum(mesh.z == nodes[1]) == 2


def test_constrained():
    from solcore.future.mesh import constrained
    import numpy as np

    nodes = np.array([0, 4, 200, 40000]) * 1e-9
    npoints = 10
    mini = 1e-9
    maxi = 100e-9
    mesh = constrained(npoints, nodes, mini, maxi)

    for m in mesh.split_in_intervals():
        assert all(np.diff(m.z.data) > mini)
        assert all(np.diff(m.z.data) < maxi)


def test_select_interval():
    from solcore.future.mesh import piecewise_uniform

    nodes = [4, 40, 400]
    npoints = [42, 100]
    mesh = piecewise_uniform(npoints, nodes)

    interval = mesh.select_interval(0)
    assert interval.ninterval == 1
    assert len(interval.z) == npoints[0]
    assert all(interval.z.interval == 0)
    assert all(interval.nodes == nodes[:2])

    interval = mesh.select_interval(1)
    assert interval.ninterval == 1
    assert len(interval.z) == npoints[1]
    assert all(interval.z.interval == 1)
    assert all(interval.nodes == nodes[1:])


def test_split_in_intervals():
    from solcore.future.mesh import piecewise_uniform

    nodes = [4, 40, 400, 4000]
    npoints = [42, 100, 200]
    mesh = piecewise_uniform(npoints, nodes)

    intervals = mesh.split_in_intervals()
    for i, interval in enumerate(intervals):
        assert interval.ninterval == 1
        assert len(interval.z) == npoints[i]
        assert all(interval.z.interval == i)
        assert all(interval.nodes == nodes[i : i + 2])


def test_concat2mesh():
    from solcore.future.mesh import uniform, piecewise_uniform, concat2mesh
    import numpy as np

    nodes1 = [4, 40, 400]
    mesh1 = piecewise_uniform([4, 40], nodes1)

    nodes2 = [4, 20]
    mesh2 = uniform(100, nodes2)

    mesh = concat2mesh(mesh1, mesh2)
    assert all(np.unique(mesh.z.interval) == (0, 1, 2))
    assert all(np.unique(mesh.z.group) == (0, 1))
    assert all(mesh.nodes == [4, 40, 400, nodes1[-1] + nodes2[-1] - nodes2[0]])
    assert len(mesh.z) == len(mesh1.z) + len(mesh2.z)
    assert max(mesh.z) == max(mesh1.z) + max(mesh2.z) - min(mesh2.z)


def test_concat_mesh():
    from solcore.future.mesh import piecewise_uniform, concat_mesh
    import numpy as np
    from functools import reduce

    rng = np.random.default_rng()
    int_per_mesh = rng.integers(2, 5)
    num_meshes = rng.integers(3, 7)

    all_nodes = []
    all_meshes = []

    for i in range(num_meshes):
        all_nodes.append(np.cumsum(rng.integers(5, 50, int_per_mesh + 1)))
        npoints = rng.integers(10, 100, int_per_mesh)
        all_meshes.append(piecewise_uniform(npoints, all_nodes[-1]))

    mesh = concat_mesh(all_meshes)
    assert all(np.unique(mesh.z.interval) == tuple(range(num_meshes * int_per_mesh)))
    assert all(np.unique(mesh.z.group) == tuple(range(num_meshes)))
    assert len(mesh.z) == reduce(lambda a, b: a + len(b.z), all_meshes, 0)
    assert all(
        mesh.nodes
        == list(
            reduce(
                lambda a, b: a + (b[1:] - b[0] + a[-1]).tolist(),
                all_nodes[1:],
                all_nodes[0].tolist(),
            )
        )
    )


def test_select_group():
    from solcore.future.mesh import uniform, piecewise_uniform, concat2mesh
    import xarray as xr

    nodes1 = [5, 40, 400]
    mesh1 = piecewise_uniform([4, 40], nodes1)

    nodes2 = [4, 20]
    mesh2 = uniform(100, nodes2)

    mesh = concat2mesh(mesh1, mesh2)
    selected = mesh.select_group(0)
    xr.testing.assert_equal(selected.z, mesh1.z)
    assert selected.nodes == approx(mesh1.nodes)

    # Only the number and relative separation of nodes and points is preserved.
    selected = mesh.select_group(1)
    assert selected.z.data == approx((mesh2.z - mesh2.z[0] + mesh1.z[-1]).data)
    assert selected.nodes == approx(mesh2.nodes - mesh2.nodes[0] + mesh1.nodes[-1])


def test_split_in_groups():
    from solcore.future.mesh import uniform, piecewise_uniform, concat2mesh
    import xarray as xr

    nodes1 = [4, 40, 400]
    mesh1 = piecewise_uniform([4, 40], nodes1)

    nodes2 = [4, 20]
    mesh2 = uniform(100, nodes2)

    mesh = concat2mesh(mesh1, mesh2)

    sel1, sel2 = mesh.split_in_groups()
    xr.testing.assert_equal(sel1.z, mesh1.z)
    assert sel1.nodes == approx(mesh1.nodes)
    assert sel2.z.data == approx((mesh2.z - mesh2.z[0] + mesh1.z[-1]).data)
    assert sel2.nodes == approx(mesh2.nodes - mesh2.nodes[0] + mesh1.nodes[-1])
