import os
import matplotlib.pyplot as plt
import numpy as np


def plot_2d_color(y_ind, calculation):
    grid_in_r = calculation['grid']['grid_in_r']
    grid_in_b = calculation['grid']['grid_in_b']
    N = calculation['N'][y_ind]
    y = calculation['grid']['grid_in_Y'][y_ind]
    y = round(y, 2)
    path = '../output/' + calculation['run_name'] + '/plots/2D_color'
    if not os.path.isdir(path):
        os.mkdir(path)
    grid_in_r = np.log10(grid_in_r)
    grid_in_b = np.log10(grid_in_b)
    # plot the color map
    plt.imshow(N, interpolation='nearest', origin='lower', aspect='auto', vmin = 0., vmax=0.4, extent=[grid_in_b[0], grid_in_b[-1], grid_in_r[0], grid_in_r[-1]])
    plt.colorbar()
    title = str(calculation['dimensionality_of_N']) + 'D BK equation at ' + str(calculation['order_of_BK']) + ' order Y=' + str(y)
    plt.title(title)

    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('b [GeV$^{-1}$]')
    plt.ylabel('r [GeV$^{-1}$]')
    plt.savefig(path + '/N_y_' + str(y) + '.pdf', format='pdf')
    plt.show()
    plt.close()
    return