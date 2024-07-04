import matplotlib.pyplot as plt
import os
from io_bk import load_calculation
import numpy as np
from initial_conds import MV_1D, GBW_1D
from matplotlib.lines import Line2D
from initial_conds import MV_1D_guillermo

def plot_1d_dasa_check():
    initial_cond = np.loadtxt('data/initial_cond.txt', delimiter=',')
    initial_cond_no_lambda = np.loadtxt('data/initial_mv_not_multipled_by_lambda.txt', delimiter=',')
    lo_y5 = np.loadtxt('data/lo_y5.txt', delimiter=',')
    nlo_y5 = np.loadtxt('data/nlo_y5.txt', delimiter=',')
    lo_y10 = np.loadtxt('data/lo_y10.txt', delimiter=',')
    nlo_y10 = np.loadtxt('data/nlo_y10.txt', delimiter=',')

    nlo_y0_guillermo = np.loadtxt('data/guillermo_nlo_y0.txt')
    nlo_y5_guillermo = np.loadtxt('data/guillermo_nlo_y5.txt')
    nlo_y10_guillermo = np.loadtxt('data/guillermo_nlo_y10.txt')

    # for run_name in ['pilot_run_ci_1D', 'pilot_run_ci_1D_simps']:
    # for run_name in ['pilot_run_NLO_1D', 'pilot_run_NLO_1D_simpson']:
    # for run_name in ['pilot_run_ci_1D', 'pilot_run_ci_1D_more_angles', 'pilot_run_ci_1D_less_angles', 'pilot_run_ci_1D_more_radius', 'pilot_run_ci_1D_more_r']:
    # for run_name in ['pilot_run_NLO_1D_more_angles', 'pilot_run_NLO_1D_less_angles', 'pilot_run_NLO_1D_more_radius', 'pilot_run_NLO_1D_less_radius', 'pilot_run_NLO_1D_more_r']:
    # for run_name in ['pilot_run_NLO_1D_less_radius', 'pilot_run_NLO_1D_evenless_radius', 'pilot_run_NLO_1D_less_radius_more_angles', 'pilot_run_NLO_1D_less_radius_less_angles']:
    # for run_name in ['pilot_run_NLO_1D_evenmore_samples', 'pilot_run_NLO_1D_evenmore_samples2', 'pilot_run_NLO_1D_evenmore_samples3', 'pilot_run_NLO_1D_evenmore_samples4', 'pilot_run_NLO_1D_evenmore_samples5', 'pilot_run_NLO_1D_evenmore_samples6', 'pilot_run_NLO_1D_evenmore_samples7']:
    # for run_name in ['pilot_run_NLO_1D_evenevenmore_samples', 'pilot_run_NLO_1D_evenevenmore_samples2', 'pilot_run_NLO_1D_evenevenmore_samples_parallel']:
    # for run_name in ['pilot_run_NLO_1D_evenmore_samples_dropdoublelog_guillermo_parallel']:
    # for run_name in ['pilot_run_NLO_1D_evenmore_samples_dropdoublelog', 'pilot_run_NLO_1D_evenevenmore_samples_dropdoublelog', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_more_angles', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_less_angles', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_more_r', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_less_r', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_less_radius', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_more_radius', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_lower_r']:
    # for run_name in ['pilot_run_NLO_1D_evenmore_samples_dropdoublelog_test', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_test_dasa']:
    # for run_name in ['pilot_run_NLO_1D_evenmore_samples_dropdoublelog_test_eta', 'pilot_run_NLO_1D_evenmore_samples_dropdoublelog_test']:
    for run_name in ['pilot_run_NLO_1D_evenmore_samples_dropdoublelog_test_lin']:
    # for run_name in ['pilot_run_NLO_1D_evenmore_samples_dropdoublelog_test']:
    # for run_name in ['pilot_run_NLO_1D_evenevenmore_samples_parallel', 'pilot_run_NLO_1D_evenevenmore_samples']:

        calculation = load_calculation(run_name)
        for y_ind in range(len(calculation['grid']['grid_in_Y'])):

            if y_ind > calculation['y_ind']:
                break
            y = round(calculation['grid']['grid_in_Y'][y_ind], 2)
            print(run_name, y_ind, y)
            Lambda = 0.241

            if y == 0.:
                plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='black')
                np.savetxt('data/marek_nlo_y0.txt', np.c_[calculation['grid']['grid_in_r'], calculation['N'][y_ind]])
            # if y == 1.:
            #     plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='blue')
            #
            # if y == 3.:
            #     plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='blue')

            if y == 5.:
                plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='green')
                np.savetxt('data/marek_nlo_y5.txt', np.c_[calculation['grid']['grid_in_r'], calculation['N'][y_ind]])

            # if y == 8.:
            #     plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='green')

            if y == 10.:
                plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='red')
                np.savetxt('data/marek_nlo_y10.txt', np.c_[calculation['grid']['grid_in_r'], calculation['N'][y_ind]])

        plt.plot(10 ** initial_cond[:, 0], initial_cond[:, 1], color='black', label='Initial Cond Dasa', linestyle='--')
        # plt.plot(Lambda * 10 ** initial_cond[:, 0], MV_1D({'grid_in_r': 10 ** initial_cond[:, 0]}, (19.*0.241)**2, Lambda=0.241, gamma=0.6), color='blue', label='MV 1D Qsq high Marek')

        # plt.plot(10 ** initial_cond_no_lambda[:, 0], initial_cond_no_lambda[:, 1], color='grey', label='Initial Cond MV no Lambda Dasa', linestyle='--')
        # plt.plot(10 ** initial_cond_no_lambda[:, 0], MV_1D({'grid_in_r': 10 ** initial_cond_no_lambda[:, 0]}, Qs0_sq=0.165, Lambda=0.241, gamma=1.135), color='blue', label='MV 1D no Lambda Marek')


    if 'NLO' in run_name:
        # plot the initial condition of Guillermo
        # grid = {'grid_in_r': np.logspace(-4., 2., 40)}
        # initial_cond_test = MV_1D_guillermo(grid, Qs0_sq=1., gamma=1.)
        # plt.plot(Lambda * grid['grid_in_r'], initial_cond_test, label='Initial Cond MV', linestyle='--', color='orange')

        plt.plot(10 ** nlo_y5[:, 0], nlo_y5[:, 1], label='NLO Y=5', linestyle='--', color='green')
        plt.plot(10 ** nlo_y10[:, 0], nlo_y10[:, 1], label='NLO Y=10', linestyle='--', color='red')
        plt.plot(Lambda * nlo_y0_guillermo[:, 0], nlo_y0_guillermo[:, 1], label='NLO Y=0 Guillermo', linestyle=':', color='black')
        plt.plot(Lambda * nlo_y5_guillermo[:, 0], nlo_y5_guillermo[:, 1], label='NLO Y=5 Guillermo', linestyle=':', color='green')
        plt.plot(Lambda * nlo_y10_guillermo[:, 0], nlo_y10_guillermo[:, 1], label='NLO Y=10 Guillermo', linestyle=':', color='red')
    else:
        plt.plot(10 ** lo_y5[:, 0], lo_y5[:, 1], label='ci Y=5', linestyle='--', color='green')
        plt.plot(10 ** lo_y10[:, 0], lo_y10[:, 1], label='ci Y=10', linestyle='--', color='red')

    # create a custom legend
    ax = plt.gca()



    legend_elements = [Line2D([0], [0], color='black', linestyle=':', label='Guillermo', lw=1.5),
                       Line2D([0], [0], color='black', linestyle='--', label='Dasa', lw=1.5),
                       Line2D([0], [0], color='black', linestyle='-', label='Marek', lw=1.5),
                       ]
    leg2 = ax.legend(handles=legend_elements, loc='upper left', shadow=False, frameon=False)
    ax.add_artist(leg2)

    plt.xscale('log')
    plt.ylabel('N(r, Y)')
    plt.xlabel(r'r$\Lambda_{QCD}$')
    # set x range
    plt.xlim(10 ** -4, 2.)

    # plt.plot(grid_in_r, N_in_r, label='Y = ' + str(y))
    # plt.legend(loc='upper left', shadow=False, frameon=False)
    # set the x-axis to log scale
    plt.title(run_name)
    if 'NLO' in run_name:
        plt.savefig('dasa_comparison_1d_nlo.pdf', format='pdf')
    else:
        plt.savefig('dasa_comparison_1d.pdf', format='pdf')
    plt.show()
    plt.close()
    return




plot_1d_dasa_check()








