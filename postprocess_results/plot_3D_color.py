import os
import matplotlib.pyplot as plt
import numpy as np


def plot_3d_color(y_ind, calculation):
    grid_in_r = calculation['grid']['grid_in_r']
    grid_in_b = calculation['grid']['grid_in_b']
    grid_in_theta = calculation['grid']['grid_in_theta']
    N = calculation['N'][y_ind]
    y = calculation['grid']['grid_in_Y'][y_ind]

    path = '../output/' + calculation['run_name'] + '/plots/3D_color'
    if not os.path.isdir(path):
        os.mkdir(path)
    grid_in_r = np.log10(grid_in_r)
    grid_in_b = np.log10(grid_in_b)
    # plot the color map
    for theta_ind in [0, len(grid_in_theta) // 4]:
        N_for_theta = N[:, :, theta_ind]
        plt.imshow(N_for_theta, interpolation='nearest', origin='lower', aspect='auto', extent=[grid_in_b[0], grid_in_b[-1], grid_in_r[0], grid_in_r[-1]])

        plt.colorbar()
        title = str(calculation['dimensionality_of_N']) + 'D BK equation at ' + str(
            calculation['order_of_BK']) + ' order Y=' + str(round(y, 2)) + r' $\theta$=' + str(
            round(grid_in_theta[theta_ind], 2))
        plt.title(title)

        # plt.xscale('log')
        # plt.yscale('log')
        plt.xlabel('b [GeV$^{-1}$]')
        plt.ylabel('r [GeV$^{-1}$]')
        plt.savefig(path + '/N_y_' + str(round(y, 2)) + '_theta' + str(round(grid_in_theta[theta_ind], 2)) + '.pdf', format='pdf')
        plt.show()
        plt.close()
    return
