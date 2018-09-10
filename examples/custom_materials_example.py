from solcore import material
from solcore import si
import matplotlib.pyplot as plt
import numpy as np

SiGeSn = material('SiGeSn')()

plt.figure()
plt.plot(si(np.arange(300, 1700, 5), 'nm')*1e9, SiGeSn.n(si(np.arange(300, 1700, 5), 'nm')))
plt.plot(si(np.arange(300, 1700, 5), 'nm')*1e9, SiGeSn.k(si(np.arange(300, 1700, 5), 'nm')))
plt.show()

