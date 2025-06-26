from solcore.spice.model import Header, Diodes, Metal, Bus, Device, Base, RearContact


def test_model_header():
    """Test Header class contains correct netlist values
    """
    # Test parameter maps to the netlist
    model = Header(temperature=100.0)
    net = model.netlist()
    degC = f"{100 - 273.15}"
    assert degC in net

    # Test default value maps to netlist
    model = Header()
    net = model.netlist()
    degC = f"{300 - 273.15}"
    assert degC in net


def test_model_diodes():
    """Test Diodes class init
    """
    model = Diodes(Eg=1.4, j01=1.0, j02=1.0, area=1.0)
    net = model.netlist()
    assert "D1" in net
    assert "D2" in net


def test_model_metal():
    """Test net labels on the metal class
    """

    # test node labels
    idx = i, j, k =  (0, 0, 0)
    model = Metal(idx, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert model.left == f"NX_{i}_{j}_{k}"
    assert model.right == f"NX_{i+1}_{j}_{k}"
    assert model.near == f"NY_{i}_{j}_{k}"
    assert model.near == f"NY_{i}_{j}_{k}"
    assert model.top == f"N_{i}_{j}_{k}"
    assert model.bottom == f"N_{i}_{j}_{k+1}"
    assert model.centre == f"N_{i}_{j}_{0}"

    # test element prefix
    assert model.element == "0_0_0"


def test_model_bus():
    """Test net labels on the Bus class
    """

    # test node labels
    idx = i, j, k =  (0, 0, 0)
    model = Bus(idx, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert model.left == f"NX_{i}_{j}_{k}"
    assert model.right == f"NX_{i+1}_{j}_{k}"
    assert model.near == f"NY_{i}_{j}_{k}"
    assert model.near == f"NY_{i}_{j}_{k}"
    assert model.top == f"in"
    assert model.bottom == f"N_{i}_{j}_{k+1}"
    assert model.centre == f"in"

    # test element prefix
    assert model.element == "0_0_0"


def test_model_device():
    """Test net labels on the Device class
    """
    idx = i, j, k =  (0, 0, 0)
    model = Device(idx, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    assert model.left == f"NX_{i}_{j}_{k}"
    assert model.right == f"NX_{i+1}_{j}_{k}"
    assert model.near == f"NY_{i}_{j}_{k}"
    assert model.near == f"NY_{i}_{j}_{k}"
    assert model.top == f"N_0_0_0"
    assert model.bottom == f"N_{i}_{j}_{k+1}"
    assert model.centre == f"N_0_0_0"

    # test element prefix
    assert model.element == "0_0_0"


def test_model_base():
    idx = (0, 0, 0)
    model = Base(idx, 1.0, 1.0)
    # netlist: connections are correct
    assert "R_BASE_Z N_0_0_0 N_0_0_1" in model.netlist() 


def test_model_rearcontact():
    idx = (0, 0, 0)
    model = RearContact(idx, 1.0, 1.0)
    # netlist: connections are correct
    assert "R_REAR_CONTACT_Z N_0_0_0 0" in model.netlist() 

def test_example_calculation():
    import numpy as np
    from solcore.structure import Junction
    from solcore.solar_cell import SolarCell
    from solcore.solar_cell_solver import solar_cell_solver
    from solcore.light_source import LightSource
    from solcore.spice.grid import HGridPattern
    from solcore.spice.netlist import generate_netlist, solve_netlist
    from solcore.spice.result import get_maximum_power_point

    temperature = 300.0

    def get_jsc(concentrationX):
        junction_model = Junction(
            kind='2D',
            T=temperature,
            reff=1,
            jref=300,
            Eg=1.4,
            A=1,
            R_sheet_top=100,
            R_sheet_bot=1e-16,
            R_shunt=1e16,
            n=3.5
        )

        solar_cell_model = SolarCell([junction_model], T=temperature)
        wl = np.linspace(350, 2000, 301) * 1e-9
        light_source = LightSource(
            source_type="standard",
            version="AM1.5g",
            x=wl,
            output_units="photon_flux_per_m",
            concentration=concentrationX
        )

        options = {
            "light_iv": True,
            "wavelength": wl,
            "light_source": light_source,
            "optics_method": "BL"
        }
        solar_cell_solver(solar_cell_model, 'iv', user_options=options)

        jsc = solar_cell_model(0).jsc
        return jsc

    def get_efficiency(concentrationX, power_in=1000.0):
        
        bus_px = 3
        fingers_px =  [2, 2, 2, 2]
        offset_px = 1
        nx, ny = 12, 12
        grid = HGridPattern(bus_px, fingers_px, offset_px=offset_px, nx=nx, ny=ny)

        # Homogeneous illumination
        illumination_map = np.ones(nx * ny).reshape((nx, ny))

        # The size of the solar is 3mm x 3mm
        size = (0.003, 0.003) # meters

        # Define a list of properies that describe each junction in the solar cell.
        # NB: currently only one junction is working.
        junctions = [
            {
                "jsc": get_jsc(concentrationX),  # solcore is calculating this for us!
                "emitter_sheet_resistance": 100.0,
                "j01": 4e-16,
                "j02": 2e-7,
                "Eg": 1.41,
                "n1": 1.0,
                "n2": 2.0
            }
        ]

        temperature = 300.0

        netlist = generate_netlist(
            grid,
            illumination_map,
            size,
            junctions,
            temperature=temperature
        )

        result = solve_netlist(netlist, temperature, 0.0, 1.5, 0.01)

        vmax, pmax, maxidx = get_maximum_power_point(result)
        
        p_per_m2 = pmax / size[0] / size[1]
        efficiency = p_per_m2 / (concentrationX * power_in)
        return efficiency

    # Get the JSC for 100x concentration
    pin = 1000 # W / m2
    jsc = get_jsc(100)
    eta = get_efficiency(jsc, power_in=pin)
    
    # The actually efficiency value is non-sense because the grid is too small 
    # to make the test run quick!
    assert eta > 0.0007
    assert eta < 0.0008


