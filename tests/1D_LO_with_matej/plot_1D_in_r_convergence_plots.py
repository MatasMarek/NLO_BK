import matplotlib.pyplot as plt
import os
from io_bk import load_calculation
import numpy as np
from initial_conds import MV_1D, GBW_1D
from matplotlib.lines import Line2D
from initial_conds import MV_1D_guillermo

def plot_1d_dasa_check():
    params_adjusted = {
        'samples_in_r': [201, 251, 281, 341],
        'samples_in_integrand_radius': [10000, 50000, 100000],
        'samples_in_integrand_angle': [400, 600, 800, 1400, 1800],
        'no_of_samples': [10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8]
    }

# It seems that good combo is 10**7 samples, 140 angles, 50000 radius, 131 r samples

    for key in params_adjusted:
        opacity_counter = 0.
        for value in params_adjusted[key]:
            opacity_counter += 0.14
            run_name = 'convergence_' + key + str(value)
            # run_name = 'convergence_lin_' + key + str(value)
            calculation = load_calculation(run_name)

            for y_ind in range(len(calculation['grid']['grid_in_Y'])):
                if y_ind > calculation['y_ind']:
                    break
                y = round(calculation['grid']['grid_in_Y'][y_ind], 2)
                print(run_name, y_ind, y)
                Lambda = 0.241

                if y == 0.:
                    plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='black')
                if y == 1.:
                    plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='blue', alpha=opacity_counter)
                if y == 5.:
                    plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='green', alpha=opacity_counter)

                if y == 10.:
                    plt.plot(calculation['grid']['grid_in_r']*Lambda, calculation['N'][y_ind], label=run_name + ' Y = ' + str(y), color='red', alpha=opacity_counter)

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

        plt.title(key)
        plt.savefig(key + '_convergence.pdf', format='pdf')
        plt.show()
        plt.close()
    return




plot_1d_dasa_check()








