import numpy as np
from scipy.integrate import quad


class Epsi1:
    def __init__(self, epsi_2, mini=0., maxi=5.):

        self.epsi_2 = epsi_2
        self.mini = max(mini, 0)
        self.maxi = maxi

    def __call__(self, E):
        eps = 1e-2
        fun = lambda w, q: w * self.epsi_2(w) / (w ** 2 - q ** 2)

        try:
            val = np.zeros_like(E)
            for i, ee in enumerate(E):
                out1 = quad(fun, a=self.mini, b=ee - eps, args=(ee,))
                out2 = quad(fun, a=ee + eps, b=self.maxi, args=(ee,))

                val[i] = out1[0] + out2[0]
        except:
            out1 = quad(fun, a=self.mini, b=E - eps, args=(E,))
            out2 = quad(fun, a=E + eps, b=self.maxi, args=(E,))

            val = out1[0] + out2[0]

        val *= 2 / np.pi

        return val


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from solcore.absorption_calculator.dielectric_constant_models import Lorentz


    def epsi2(x):
        En = 2.5
        Br = 0.5
        A = 1

        s = 0.5 * Br / np.sqrt(np.log(2))

        out1 = np.exp(-((x - En) / s) ** 2)
        out2 = np.exp(-((x + En) / s) ** 2)

        out = A * (out1 - out2)

        return out


    epsi1 = Epsi1(epsi2)

    E = np.linspace(0, 5, 50)

    e1 = epsi1(E)
    e2 = epsi2(E)

    n = Lorentz(0.2, 2.5, 0.5).dielectric(1240 / E)
    e1_l = np.real(n)
    e2_l = np.imag(n)

    # print(e1, e2)

    plt.plot(E, e1, 'b', label='e1')
    plt.plot(E, e2, 'r', label='e2')
    plt.plot(E, e1_l, 'b--', label='e1_l')
    plt.plot(E, e2_l, 'r--', label='e2_l')
    plt.legend()
    plt.show()
