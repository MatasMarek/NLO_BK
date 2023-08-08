import numpy as np
import matplotlib.pyplot as plt
from kernel import alpha_s


def plot_alpha_s():
    rs = np.logspace(-7, 2, 225)
    alpha_s_values = alpha_s(rs)
    plt.plot(rs, alpha_s_values)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('r [GeV$^{-1}$]')
    plt.ylabel(r'$\alpha_s$')
    plt.show()
    plt.close()

plot_alpha_s()


