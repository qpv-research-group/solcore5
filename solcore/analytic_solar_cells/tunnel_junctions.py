import numpy as np
from solcore.constants import kb, q


def resistive_tunnel_junction(junction, options):
    """ Calculates the IV curve of a tunnel junction when it is modelled as a simple resistor. The minimum resistance of the junction is always 1e-16.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    def iv(v):
        return v / junction.R

    def vi(j):
        return j * junction.R

    junction.voltage = options.internal_voltages
    junction.current = iv(junction.voltage)
    junction.iv = iv
    junction.vi = vi


def parametric_tunnel_junction(junction, options):
    """ Calculates the IV curve of a tunnel junction when it is modelled using a set of empirical parameters, such as peak curent and voltage, valley current and voltage, etc. The total current of the tunnel junction is summ of 3 components: the tunnel current, the excess current and the difussion current.

    :param junction: A junction object.
    :param options: Solver options.
    :return: None.
    """

    T = options.T
    junction.voltage = options.internal_voltages

    try:
        # First we calculate the tunnel current
        jp = junction.j_peak
        vp = junction.v_peak

        def tunnel_current(v):
            return jp * v / vp * np.exp(1 - v / vp)

        junction.tunnel_current = tunnel_current

        # Now we calculate the excess current due to carrier tunnelling by way of states within the forbidden gap
        jva = junction.j_valley
        vv = junction.v_valley
        pref = junction.prefactor

        def excess_current(v):
            return jva * np.exp(pref * (v - vv))

        junction.excess_current = excess_current

        # Finally we calculate the diffusion current
        j01 = junction.j01

        def diffusion_current(v):
            return j01 * np.exp(q * v / kb / T)

        junction.diffusion_current = diffusion_current
    except AttributeError:
        raise

    # And we put everything together in an IV curve function
    def iv(v):
        return junction.tunnel_current(v) + junction.excess_current(v) + junction.diffusion_current(v)

    junction.current = iv(junction.voltage)

    # Also, we calculate the inverse, needed for the calculation of the IV curve in MJ solar cells
    MJ_current = np.where(junction.voltage <= junction.v_peak, junction.current,
                          np.maximum(iv(junction.v_peak), junction.current))

    # And the corresponding function doing that
    def vi(j):
        return np.interp(j, MJ_current, junction.voltage)

    junction.iv = iv
    junction.vi = vi


def example_resistive_tunnel_junction(show=True):
    from solcore.structure import TunnelJunction
    from solcore.solar_cell_solver import default_options
    import matplotlib.pyplot as plt

    # Parametric tunnel example
    my_tunnel = TunnelJunction(R=0.05)
    resistive_tunnel_junction(my_tunnel, default_options)

    if show:
        v = my_tunnel.voltage

        plt.plot(v, my_tunnel.current, 'k', linewidth=2, label='Total')

        plt.legend(fontsize=12)
        plt.ylim(0, 1.5)
        plt.xlim(0, 1)
        plt.ylabel('Current Density(A/$m^2$)', fontsize=12)
        plt.xlabel('Voltage(V)', fontsize=12)
        plt.tick_params(labelsize=12)
        plt.tight_layout()
        plt.show()

    return my_tunnel


def example_parametric_tunnel_junction(show=True):
    from solcore.structure import TunnelJunction
    from solcore.solar_cell_solver import default_options
    import matplotlib.pyplot as plt

    # Parametric tunnel example
    my_tunnel = TunnelJunction(v_peak=0.1, j_peak=1, v_valley=0.5, j_valley=0.1, prefactor=10, j01=1e-11)
    parametric_tunnel_junction(my_tunnel, default_options)

    if show:
        v = my_tunnel.voltage

        plt.plot(v, my_tunnel.tunnel_current(v), 'r--', label='Tunnel')
        plt.plot(v, my_tunnel.excess_current(v), 'g--', label='Excess')
        plt.plot(v, my_tunnel.diffusion_current(v), 'b--', label='Diffusion')
        plt.plot(v, my_tunnel.current, 'k', linewidth=2, label='Total')

        plt.legend(fontsize=12)
        plt.ylim(0, 1.5)
        plt.xlim(0, 1)
        plt.ylabel('Current Density(A/$m^2$)', fontsize=12)
        plt.xlabel('Voltage(V)', fontsize=12)
        plt.tick_params(labelsize=12)
        plt.tight_layout()
        plt.show()

    return my_tunnel


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from solcore.structure import TunnelJunction
    from solcore.solar_cell_solver import default_options

    # tunnel = example_parametric_tunnel_junction(False)
    # example_resistive_tunnel_junction()

    tunnel = TunnelJunction(v_peak=0.5, j_peak=100, v_valley=0.5, j_valley=0, prefactor=1, j01=1e-23, kind='parametric')
    parametric_tunnel_junction(tunnel, default_options)

    v = tunnel.voltage

    I = np.linspace(0, 20, 10)

    plt.plot(v, tunnel.tunnel_current(v), 'r--', label='Tunnel')
    plt.plot(v, tunnel.excess_current(v), 'g--', label='Excess')
    plt.plot(v, tunnel.diffusion_current(v), 'b--', label='Diffusion')
    plt.plot(v, tunnel.current, 'k', linewidth=2, label='Total')
    # plt.plot(v, MJ_current, 'k-o', linewidth=2, label='MJ current')
    plt.plot(tunnel.vi(I), I, linewidth=2, label='Test')

    plt.legend(fontsize=12)
    plt.ylim(0, 150)
    plt.xlim(0, 2.5)
    plt.ylabel('Current Density(A/$m^2$)', fontsize=12)
    plt.xlabel('Voltage(V)', fontsize=12)
    plt.tick_params(labelsize=12)
    plt.tight_layout()
    plt.show()
