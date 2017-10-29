# Example 1 - Virtual experiments
# -------------------------------
# We use the device created in PDD_example_1.py in a couple of virtual experiments

import solcore.poisson_drift_diffusion as PDD
import matplotlib.pyplot as plt
import os

# First we load the device structure stored in MyDevice.json
this_dir = os.path.split(__file__)[0]
my_file = 'MyDevice'
full_path = os.path.join(this_dir, my_file)
MyDevice = PDD.Load(full_path)

# Then we run two virtual experiments, just calculate the structure under equilibrum and under short circuit conditions.
# We use the default settings of the solver.
EQ = PDD._Equilibrium(MyDevice)
SC = PDD._ShortCircuit(MyDevice)

# Finally, we plot the carrier densities in both cases using the information stored in the output dictionaries
plt.semilogy(EQ['Bandstructure']['x'] * 1e9, EQ['Bandstructure']['n'], 'b', label=' n equilibrium')
plt.semilogy(EQ['Bandstructure']['x'] * 1e9, EQ['Bandstructure']['p'], 'r', label=' p equilibrium')
plt.semilogy(SC['Bandstructure']['x'] * 1e9, SC['Bandstructure']['n'], 'c', label=' n short circuit')
plt.semilogy(SC['Bandstructure']['x'] * 1e9, SC['Bandstructure']['p'], 'm', label=' p short circuit')
plt.legend(loc='lower right')
plt.ylabel('Carrier densities (m-3)')
plt.xlabel('Position (nm)')

plt.show()
