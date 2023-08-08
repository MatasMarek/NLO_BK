import os
import matplotlib.pyplot as plt
import numpy as np


def plot_2d_fixed_b(y_ind, calculation):
    grid_in_r = calculation['grid']['grid_in_r']
    grid_in_b = calculation['grid']['grid_in_b']
    N = calculation['N'][y_ind]
    y = calculation['grid']['grid_in_Y'][y_ind]
    path = '../output/' + calculation['run_name'] + '/plots/2D_fixed_b'
    if not os.path.isdir(path):
        os.mkdir(path)
    # find index closest to b = 1
    b_ind_1 = (np.abs(grid_in_b - 1.)).argmin()
    b_ind_2 = (np.abs(grid_in_b - 1.5)).argmin()

    for b_ind in [0, b_ind_1, b_ind_2, -5]:
        b = grid_in_b[b_ind]
        # round b
        if b > 0.01:
            b = round(b, 2)
        else:
            b = round(b, 5)
        N_in_r = N[:, b_ind]
        plt.plot(grid_in_r, N_in_r, label='Y = ' + str(y) + ', b = ' + str(b) + ' GeV$^{-1}$')

    plt.legend(loc='upper left', shadow=False, frameon=False)
    # set the x-axis to log scale
    title = str(calculation['dimensionality_of_N']) + 'D BK equation at ' + str(calculation['order_of_BK']) + ' order'
    plt.title(title)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('N(r, b, Y)')
    plt.xlabel('r [GeV$^{-1}$]')
    plt.savefig(path + '/N_y_' + str(y) + '.pdf', format='pdf')
    plt.show()
    plt.close()
    return