import os
import numpy as np
from io_bk import load_calculation
from plot_1D_in_r import plot_1d_in_r
from plot_2D_fixed_b import plot_2d_fixed_b
from plot_2D_color import plot_2d_color
from plot_3D_color import plot_3d_color
from plot_4D_color import plot_4d_color


def run_all():
    print('--- Postprocessing has started ---')
    for dir in os.listdir('../output'):
        if not dir.startswith('.') and os.path.isfile('../output/' + dir + '/calculation.npy'):
            print(dir)
            run(dir)


def run(run_name):
    if not os.path.isdir('../output/' + run_name + '/plots'):
        os.mkdir('../output/' + run_name + '/plots')
    calculation = load_calculation(run_name)
    for y_ind in range(len(calculation['grid']['grid_in_Y'])):
        if y_ind > calculation['y_ind']:
            break
        y = round(calculation['grid']['grid_in_Y'][y_ind], 2)
        if calculation['dimensionality_of_N'] == 1 and y % 1. == 0.:
            plot_1d_in_r(calculation['N'][y_ind], calculation['grid']['grid_in_r'], y, calculation)
        if calculation['dimensionality_of_N'] == 2 and y % 1. == 0:
            plot_2d_color(y_ind, calculation)
            plot_2d_fixed_b(y_ind, calculation)
        if calculation['dimensionality_of_N'] == 3:
            plot_3d_color(y_ind, calculation)
            print('plotting 3D for y = ' + str(y))
        if calculation['dimensionality_of_N'] == 4:
            plot_4d_color(y_ind, calculation)
            print('plotting 4D for y = ' + str(y))


if __name__ == '__main__':
    # run_all()
    # run('pilot_run_LO_2D')
    # run('pilot_run_ci_1D')
    # run('pilot_run_ci_2D')
    # run('pilot_run_NLO_1D')
    # run('pilot_run_NLO_2D')
    # run('pilot_run_ci_3D')
    run('pilot_run_ci_4D')



