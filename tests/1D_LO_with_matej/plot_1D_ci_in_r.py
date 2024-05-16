import matplotlib.pyplot as plt
import os
from io_bk import load_calculation
import numpy as np
from initial_conds import MV_1D, GBW_1D
from matplotlib.lines import Line2D
from initial_conds import MV_1D

def plot_1d_matej_check():
    # grid = {'grid_in_r': np.logspace(-4., 2., 40)}
    # initial_cond_test = MV_1D(grid, Qs0_sq=0.165, gamma=1.135, Lambda=0.241)
    # plt.plot(grid['grid_in_r'], initial_cond_test, label='Initial Cond MV', color='orange')

    initial_cond = np.loadtxt('data_ci/Y=00.00_honzas_method.csv', delimiter=',')
    lo_y1 = np.loadtxt('data_ci/Y=01.00_honzas_method.csv', delimiter=',')

    initial_cond_y = np.loadtxt('data_ci/Y=00.00.csv', delimiter=',')
    lo_y1_y = np.loadtxt('data_ci/Y=01.00.csv', delimiter=',')

    initial_cond_y_after_revision = np.loadtxt('data_ci/Y=00.00_after_revision.csv', delimiter=',')
    lo_y1_y_after_revision = np.loadtxt('data_ci/Y=01.00_after_revision.csv', delimiter=',')


    # for run_name in ['pilot_run_ci_1D_rk4', 'pilot_run_ci_1D']:
    for run_name in ['pilot_run_ci_1D']:

        calculation = load_calculation(run_name)
        for y_ind in range(len(calculation['grid']['grid_in_Y'])):

            if y_ind > calculation['y_ind']:
                break
            y = round(calculation['grid']['grid_in_Y'][y_ind], 2)
            print(run_name, y_ind, y)
            if 'NLO' in run_name:
                linestyle = ':'
                linewidth = 3
            else:
                linestyle = '-'
                linewidth = 1

            # if y == 0.:
            #     plt.plot(calculation['grid']['grid_in_r'], calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='black', linestyle=linestyle, linewidth=linewidth)
            #
            # if y == 1.:
            #     plt.plot(calculation['grid']['grid_in_r'], calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='green', linestyle=linestyle, linewidth=linewidth)
            #
            # if y == 10.:
            #     plt.plot(calculation['grid']['grid_in_r'], calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='red', linestyle=linestyle, linewidth=linewidth)


    plt.plot(initial_cond_y[:, 1], initial_cond_y[:, 0], color='black', label='Matej Y=0', linestyle='--')
    plt.plot(lo_y1_y[:, 1], lo_y1_y[:, 0], label='Matej Y=1', linestyle='--', color='green')

    plt.plot(initial_cond[:, 1], initial_cond[:, 0], color='black', label='Matej old method Y=0', linestyle=':', linewidth=2)
    plt.plot(lo_y1[:, 1], lo_y1[:, 0], label='Matej old method Y=1', linestyle=':', color='green', linewidth=2)

    plt.plot(initial_cond_y_after_revision[:, 1], initial_cond_y_after_revision[:, 0], color='black', label='Matej after revision Y=0', linestyle='-', linewidth=2)
    plt.plot(lo_y1_y_after_revision[:, 1], lo_y1_y_after_revision[:, 0], label='Matej after revision Y=1', linestyle='-', color='green', linewidth=2)


    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('N(r, Y)')
    plt.xlabel(r'r$\Lambda_{QCD}$')
    # set x range
    plt.xlim(4*10 **-6, 10**1)
    plt.legend(loc='upper left', shadow=False, frameon=False)
    plt.title(run_name)
    # plt.savefig('dasa_comparison_1d.pdf', format='pdf')
    plt.show()
    plt.close()
    return


plot_1d_matej_check()








