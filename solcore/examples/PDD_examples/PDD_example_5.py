# Example 4 - Virtual experiments
# -------------------------------
# Finally, we calculate the quantum efficiency, plotting the internal and the external ones.

import solcore.poisson_drift_diffusion as PDD
import matplotlib.pyplot as plt
import os

# First we load the device structure stored in MyDevice.json
this_dir = os.path.split(__file__)[0]
my_file = 'MyDevice'
full_path = os.path.join(this_dir, my_file)
MyDevice = PDD.Load(full_path)

# We use the default settings of the solver to calculate the QE.
QE = PDD._QE(MyDevice)

# Finally, we plot the internal and external quantum efficiencies using the information stored in the output dictionaries
plt.plot(QE['QE']['wavelengths'] / 1e-9, QE['QE']['IQE'] * 100, label='IQE')
plt.plot(QE['QE']['wavelengths'] / 1e-9, QE['QE']['EQE'] * 100, label='EQE')
plt.ylim(0, 100)
plt.legend(loc='lower left')
plt.ylabel('QE (%)')
plt.xlabel('Wavelength (nm)')


plt.show()
