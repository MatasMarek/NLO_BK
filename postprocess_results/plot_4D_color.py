import os
import matplotlib.pyplot as plt
import numpy as np


def plot_4d_color(y_ind, calculation):
    grid_in_r = calculation['grid']['grid_in_r']
    grid_in_b = calculation['grid']['grid_in_b']
    grid_in_theta = calculation['grid']['grid_in_theta']
    grid_in_phi = calculation['grid']['grid_in_phi']
    N = calculation['N'][y_ind]
    y = calculation['grid']['grid_in_Y'][y_ind]

    path = '../output/' + calculation['run_name'] + '/plots/3D_color'
    if not os.path.isdir(path):
        os.mkdir(path)
    grid_in_r = np.log10(grid_in_r)
    grid_in_b = np.log10(grid_in_b)
    # plot the color map
    for phi_ind in [len(grid_in_phi) // 4,len(grid_in_phi) // 4 + len(grid_in_phi) // 2]:
        phi = round(grid_in_phi[phi_ind], 2)
        path = '../output/' + calculation['run_name'] + '/plots/3D_color/phi_' + str(phi)
        if not os.path.isdir(path):
            os.mkdir(path)
        for theta_ind in [0, len(grid_in_theta) // 4]:
            N_for_theta_and_phi = N[:, :, theta_ind, phi_ind]
            plt.imshow(N_for_theta_and_phi, interpolation='nearest', origin='lower', aspect='auto', extent=[grid_in_b[0], grid_in_b[-1], grid_in_r[0], grid_in_r[-1]])

            plt.colorbar()
            title = str(calculation['dimensionality_of_N']) + 'D BK equation at ' + str(calculation['order_of_BK']) + ' order Y=' + str(round(y,2)) + r' $\theta$=' + str(round(grid_in_theta[theta_ind], 2)) + r' $\phi$=' + str(phi)
            plt.title(title)

            # plt.xscale('log')
            # plt.yscale('log')
            plt.xlabel('b [GeV$^{-1}$]')
            plt.ylabel('r [GeV$^{-1}$]')
            plt.savefig(path + '/N_y_' + str(round(y, 2)) + '_theta' + str(round(grid_in_theta[theta_ind], 2)) + '.pdf', format='pdf')
            plt.show()
            plt.close()
    return
