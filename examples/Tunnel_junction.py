import matplotlib.pyplot as plt

from solcore.structure import TunnelJunction
from solcore.solar_cell_solver import default_options
from solcore.analytic_solar_cells import parametric_tunnel_junction

# We define the tunnel junction and solve its properties
tunnel = TunnelJunction(v_peak=0.2, j_peak=100, v_valley=0.9, j_valley=10, prefactor=10, j01=1e-21, kind='parametric')
parametric_tunnel_junction(tunnel, default_options)

v = tunnel.voltage

plt.plot(v, tunnel.tunnel_current(v), 'r--', label='Tunnel')
plt.plot(v, tunnel.excess_current(v), 'g--', label='Excess')
plt.plot(v, tunnel.diffusion_current(v), 'b--', label='Diffusion')
plt.plot(v, tunnel.current, 'k', linewidth=3, color='DimGray', label='Total')
plt.plot((0.2, 0.9), (100, 10), 'ko')

plt.annotate('V$_P$, J$_P$', xy=(0.2, 110), fontsize=12)
plt.annotate('V$_V$, J$_V$', xy=(0.6, 10), fontsize=12)

plt.legend(fontsize=12, frameon=False)
plt.ylim(0, 150)
plt.xlim(0, 2)
plt.ylabel('Current Density(A/$m^2$)', fontsize=12)
plt.xlabel('Voltage(V)', fontsize=12)
plt.tick_params(labelsize=12)
plt.tight_layout()
plt.show()
