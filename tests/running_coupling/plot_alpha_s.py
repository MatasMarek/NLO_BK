import numpy as np
import matplotlib.pyplot as plt
from kernel import alpha_s
import const as c

c.C = 2.586
# c.C = 9.
c.alpha_freeze = 1.
c.ml = 0.14  # GeV
c.mc = 1.27  # GeV
c.mb = 4.2  # GeV

c.lambdaQCD5 = c.mZ * np.exp(-2. * np.pi / (c.alpha_mz * c.beta5))
c.lambdaQCD4 = np.power(c.mb, 1. - (c.beta5 / c.beta4)) * np.power(c.lambdaQCD5, c.beta5 / c.beta4)
c.lambdaQCD3 = np.power(c.mc, 1. - (c.beta4 / c.beta3)) * np.power(c.lambdaQCD4, c.beta4 / c.beta3)

c.rsat_sq = 4. * c.C ** 2 / (c.lambdaQCD3 ** 2 * np.exp(4. * np.pi / (c.beta3 * c.alpha_freeze)))

c.rc_sq = 4. * c.C ** 2 / (c.mc * 2)
c.rb_sq = 4. * c.C ** 2 / (c.mb * 2)

def plot_alpha_s():

    stepan_y = np.load('/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/tests/running_coupling/alpha_rc_y.npy')
    stepan_x = np.load('/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/tests/running_coupling/alpha_rc_x.npy')
    plt.plot(stepan_x, stepan_y, label='Stepan Python')

    matej = np.loadtxt('/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/tests/running_coupling/alpha_s.csv', delimiter=',')
    plt.plot(matej[:, 0], matej[:, 1], label='Matej')

    rs_sq = np.logspace(-5, 2, 225)**2
    alpha_s_values = alpha_s(rs_sq)
    plt.plot(np.sqrt(rs_sq), alpha_s_values, label='Marek')
    plt.legend()

    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('r [GeV$^{-1}$]')
    plt.ylabel(r'$\alpha_s$')
    # set y scale from 0.1 to 1
    plt.ylim(0., 1.1)
    plt.show()
    plt.close()

plot_alpha_s()


def plot_k_ci1():
    stepan = np.load('/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/tests/running_coupling/lo_stepan.npy')
    stepan_x = np.load('/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/tests/running_coupling/lo_stepan_x.npy')
    plt.plot(stepan_x, stepan*np.pi, label='Stepan Python')

    ja_digitized_again = np.loadtxt('/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/tests/running_coupling/kc1_digitized_dissertation.txt', delimiter=',')
    plt.plot(10**ja_digitized_again[:, 0], 10**ja_digitized_again[:, 1], 'o', label='Digitized again')

    r1s = np.logspace(-5, 2, 225)
    r = 1.
    theta = np.pi/2
    K_ci1 = np.zeros(len(r1s))
    min_r = np.zeros(len(r1s))
    for i in range(len(r1s)):
        r1 = r1s[i]
        r2 = np.sqrt(r**2 + r1**2 - 2*r*r1*np.cos(theta))
        min_r[i] = min(r, r1, r2)
        K_ci1[i] = 1./(2.*np.pi) * r**2/(r1**2 * r2**2)
    alpha = alpha_s(min_r ** 2)
    K_ci1 *= alpha
    # K_ci1 *= c.Nc/np.pi
    plt.plot(r1s, K_ci1, label='Marek Pythonic version')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('r1 [GeV$^{-1}$]')
    plt.ylabel(r'K$_{ci}$ LO * $\alpha_s$ term')
    plt.show()
    plt.close()

plot_k_ci1()



