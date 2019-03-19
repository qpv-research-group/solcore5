import matplotlib.pyplot as plt
import numpy as np
from solcore.absorption_calculator import calculate_rat
from solcore.absorption_calculator.dielectric_constant_models import DielectricConstantModel, Drude, Poles

x = 2 * np.logspace(3, 4, 200)

# Model parameters for three ITO layers grown at different temperatures. They have been obtained after fitting the
# ellipsometry data.
e_inf = [3.7883, 3.8915, 3.8982]
A = [16.038, 36.556, 36.806]
Br = [0.11112, 0.10413, 0.088618]
label = ['150', '250', '350']
ls = ['solid', 'dashed', 'dashdot']

for i in range(len(e_inf)):

    # We create the oscillators for each layer
    drud = Drude(An=A[i], Brn=Br[i])

    # Then we put them together inside a dielectric model
    model = DielectricConstantModel(e_inf=e_inf[i], oscillators=[drud])

    # We might want to calculate the RAT of the films
    out = calculate_rat([[300, model]], x)

    plt.figure(1)
    plt.plot(x/1000, out['R'], 'b', label='R ' + label[i], ls=ls[i])
    plt.plot(x/1000, out['A'], 'r', label='A ' + label[i], ls=ls[i])
    plt.plot(x/1000, out['T'], 'g', label='T ' + label[i], ls=ls[i])

    # And also want to know the n and k data
    n = model.n_and_k(x)

    plt.figure(2)
    plt.plot(x/1000, np.real(n), 'b', label='n ' + label[i], ls=ls[i])
    plt.plot(x/1000, np.imag(n), 'g', label='k ' + label[i], ls=ls[i])


plt.figure(1)
plt.xlabel('Wavelength (µm)')
plt.ylabel('RAT')
plt.xlim((2, 20))
plt.legend(loc=0)

plt.figure(2)
plt.xlabel('Wavelength (µm)')
plt.ylabel('n & k')
plt.xlim((2, 20))
plt.legend(loc=2)

plt.show()
