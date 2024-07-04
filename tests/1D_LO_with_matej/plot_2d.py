
import matplotlib.pyplot as plt
import os
from io_bk import load_calculation
import numpy as np
from initial_conds import MV_1D, GBW_1D
from matplotlib.lines import Line2D
from initial_conds import MV_1D_guillermo
from postprocess_results.plot_2D_color import plot_2d_color

# run_name = 'NLO_2D_linear_grid'
# run_name = 'NLO_2D_identical_init'
# run_name = 'NLO_2D_heikki_cutoff'
# run_name = 'NLO_2D_parallel'
# run_name = 'NLO_2D_same_init_fewsmallbs_nonparallel'
# run_name = 'NLO_2D_same_init_fewsmallbs'
# run_name = 'NLO_2D_same_init_allbs'
# run_name = 'NLO_2D_same_init_just_zero'
run_name = 'NLO_2D'
# run_name = 'convergence_lin_' + key + str(value)
calculation = load_calculation(run_name)
# plot only three plots, first one is the initial condition,
if calculation['y_ind'] > 8:
    ys_to_iterate_over = [0, int(calculation['y_ind']/3.), int(calculation['y_ind']/2.), int(calculation['y_ind']/1.5), calculation['y_ind']-1]
else:
    ys_to_iterate_over = range(calculation['y_ind'])
linestyles = ['-', '--', '-.', ':', '-']

#DEBUG
ys_to_iterate_over = range(calculation['y_ind'])

for y_ind in ys_to_iterate_over:
    # Lambda = 0.241
    # # for b in [0]:
    # style = 0
    # for b in range(len(calculation['grid']['grid_in_b'])):
    # # for b in [len(calculation['grid']['grid_in_b'])-1]:
    #     plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind, :, b], label='b = ' + str(calculation['grid']['grid_in_b'][b]), linestyle=linestyles[style])
    #     style += 1
    # plt.xscale('log')
    # plt.ylabel('N(r, Y)')
    # plt.xlabel(r'r$\Lambda_{QCD}$')
    # # set x range
    # plt.xlim(10 ** -4, 2.)
    # plt.title(calculation['run_name'] + ' at Y = ' + str(calculation['grid']['grid_in_Y'][y_ind]))
    # plt.legend()
    # plt.show()
    # plt.close()

    plot_2d_color(y_ind, calculation)


