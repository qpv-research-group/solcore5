# Example 3 - Virtual experiments
# -------------------------------
# Now we calculate the light IV curve under the AM1.5d spectrum assuming Rs = 0 Ohm m2 and Rs = 0.001 Ohm m2

import solcore.poisson_drift_diffusion as PDD
import matplotlib.pyplot as plt
import os

# First we load the device structure stored in MyDevice.json
this_dir = os.path.split(__file__)[0]
my_file = 'MyDevice'
full_path = os.path.join(this_dir, my_file)
MyDevice = PDD.Load(full_path)

# We use the default settings of the solver to calculate the light IV.
IV = PDD.IV(MyDevice, vfin=1.2, vstep=0.05, light=True)
IVrs = PDD.IV(MyDevice, vfin=1.2, vstep=0.05, light=True, rs=0.001)

# Finally, we plot the two curves
plt.plot(IV['IV']['V'], -IV['IV']['J'], label='Rs = 0 Ohm m2')
plt.plot(IVrs['IV']['V'], -IVrs['IV']['J'], label='Rs= 0.001 Ohm m2')
plt.ylim(0, 1.1 * max(-IV['IV']['J']))
plt.legend(loc='lower left')
plt.ylabel('Current density (A/m2)')
plt.xlabel('Voltage (V)')

plt.show()
