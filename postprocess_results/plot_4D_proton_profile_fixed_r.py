import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps


def plot_4d_proton_profile_fixed_r(y_ind, calculation):
    grid_in_r = calculation['grid']['grid_in_r']
    grid_in_b = calculation['grid']['grid_in_b']
    grid_in_theta = calculation['grid']['grid_in_theta']
    grid_in_phi = calculation['grid']['grid_in_phi']
    N = calculation['N'][y_ind]
    y = calculation['grid']['grid_in_Y'][y_ind]

    # Find index of r closest to a set value
    r_set = 10**-2
    r_ind = np.argmin(np.abs(grid_in_r - r_set))


    path = '../output/' + calculation['run_name'] + '/plots/proton_profile'
    if not os.path.isdir(path):
        os.mkdir(path)
    grid_in_r = np.log10(grid_in_r)
    grid_in_b = np.log10(grid_in_b)
    # plot the color map
    radial_N = np.zeros((len(grid_in_b), len(grid_in_phi)))
    for b_ind in range(len(grid_in_b)):
        for phi_ind in range(len(grid_in_phi)):
            # integrate N over all thetas
            radial_N[b_ind, phi_ind] = simps(N[r_ind, b_ind, :, phi_ind], grid_in_theta)
    # plot a 2D plot of the radial_N in the polar projection

    # print(radial_N)
    plt.subplot(projection="polar")
    # Add text to the plot
    # plt.text(0, -2.5, r'log(b) [GeV$^{-1}$]', rotation=90, horizontalalignment='center', verticalalignment='center', fontsize=12)
    # plt.ylabel(r'log(b) [GeV$^{-1}$]')
    plt.xlabel(r'$\phi$ [deg], log(b) [GeV$^{-1}$]')
    plt.title(r'$\int$ N(r=' + str(r_set) + r', log(b), $\theta$, $\phi$) d$\theta$ at Y=' + str(round(y, 2)))
    plt.pcolor(grid_in_phi, grid_in_b, radial_N, shading='auto', cmap='coolwarm', vmin=0., vmax=6*10**-5)
    # plt.pcolor(grid_in_phi, grid_in_b, radial_N, shading='auto', cmap='coolwarm')
    # plt.colorbar(orientation='vertical', shrink=0.9, label=r'$\int$ N(r, b, $\phi$)d$\theta$')
    plt.colorbar(orientation='vertical', shrink=0.9)

    plt.savefig(path + '/N_y_' + str(round(y, 2)) + '.pdf', format='pdf')
    plt.show()
    plt.close()
    return
