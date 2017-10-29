# Example 2 - Virtual experiments
# -------------------------------
# Using the same structure that before, we calculate the dark current, plotting the different components that contribute to it.

import solcore.poisson_drift_diffusion as PDD
import matplotlib.pyplot as plt
import os

# First we load the device structure stored in MyDevice.json
this_dir = os.path.split(__file__)[0]
my_file = 'MyDevice'
full_path = os.path.join(this_dir, my_file)
MyDevice = PDD.Load(full_path)

# We use the default settings of the solver to calculate the dark IV.
IV = PDD._IV(MyDevice, vfin=1.2, vstep=0.05)

# Finally, we plot the different components of the dark current using the information stored in the output dictionaries
plt.semilogy(IV['IV']['V'], IV['IV']['J'], 'o', label='Jtot')
plt.semilogy(IV['IV']['V'], IV['IV']['Jrad'], label='Jrad')
plt.semilogy(IV['IV']['V'], IV['IV']['Jsrh'], label='Jsrh')
plt.semilogy(IV['IV']['V'], IV['IV']['Jsur'], label='Jsur')
plt.legend(loc='lower right')
plt.ylabel('Current density (A/m2)')
plt.xlabel('Voltage (V)')

plt.show()
