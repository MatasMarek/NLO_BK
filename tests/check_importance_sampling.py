import matplotlib.pyplot as plt
from run import run_calculation
import numpy as np
from integrand import integrand
from run import set_up_grids
from integrate import rs_bs_and_variables_for_N, get_Ns
from scipy.integrate import simps

print('Importance sampling test initiated')

steps_in_integrand_theta = 20
shift = 2.*np.pi/steps_in_integrand_theta/2.  # to avoid double counting and y-axis with z and w.
grid = {
    'grid_in_Y': np.linspace(0., 10., 201),
    'grid_in_r': np.logspace(-8., 2., 10),
    'grid_in_b': np.logspace(-1., 2., 10),
    'grid_in_integrand_radius': np.logspace(-7., 2., 200),
    'grid_in_integrand_angle': np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta),  # not to include 2pi to avoid double counting
}

dimensionality_of_N = 2  # r, b

# integration_method = 'MC'
integration_method = 'Simps'

no_of_samples = 10**4

order_of_rk = 1
# order_of_BK = 'ci'
order_of_BK = 'NLO'
number_of_cores = 1

from initial_conds import mareks_N as cond
initial_cond = cond(grid, dimensionality_of_N)
run_name = 'pilot_run_ci_2D'


calculation = {
    'dimensionality_of_N': dimensionality_of_N,  # Dimensionality of N
    'initial_cond': initial_cond,  # Form of the initial condition
    'integration_method': integration_method,  # Integration method (MC/Simpson)
    'order_of_rk': order_of_rk, # Order of RK method?
    'grid': grid,  # The grid
    'order_of_BK': order_of_BK,  # Order of the BK equation (LO/NLO)
    'number_of_cores': number_of_cores,  # number of cores for the outermost parallelization
    'run_name': run_name,  # name of the run
    'no_of_samples': no_of_samples,  # for the stochastical MC case
}

calculation = set_up_grids(calculation)
calculation['N'][0] = calculation['initial_cond']

grid_in_r = calculation['grid']['grid_in_r']
calculation['variables_for_N'] = {'xy': {}}
calculation['Ns'] = {}
grid_in_b = calculation['grid']['grid_in_b']
calculation['y_ind'] = 1
y_ind = calculation['y_ind']
z_coods = calculation['integrand_cartesian_coods']

integrand_angles = calculation['grid']['grid_in_integrand_angle']
integrand_radii = calculation['grid']['grid_in_integrand_radius']

print('Loop started')
for r_ind in range(len(grid_in_r)):
    r = grid_in_r[r_ind]
    calculation['variables_for_N']['xy']['r'] = r
    calculation['variables_for_N']['xy']['rsq'] = r ** 2
    # for b_ind in range(len(grid_in_b)):
    for b_ind in [0]:
        b = grid_in_b[b_ind]
        calculation['Ns']['xy'] = calculation['N'][y_ind - 1][r_ind][b_ind]
        calculation['variables_for_N']['xy']['b'] = b
        if calculation['dimensionality_of_N'] == 2:
            calculation['Ns']['xy'] = calculation['N'][y_ind - 1][r_ind][b_ind]
            x = np.array([-r / 2., b])
            y = np.array([r / 2., b])
            calculation = rs_bs_and_variables_for_N(calculation, x, y, 0)
            calculation = get_Ns(calculation)
            if calculation['order_of_BK'] == 'NLO' or calculation['order_of_BK'] == 'NLO_LOcut':
                integrand_simps_over_z, integrand_simps_over_wz = integrand(calculation)
                # Integral over z
                integrand_simps_over_z = integrand_simps_over_z[:len(calculation['integrand_cartesian_coods'])]  # Make it only for one integrand of z instead of the repeating tiled version
                integrand_simps_over_z = integrand_simps_over_z.reshape((len(integrand_radii), len(integrand_angles)))
                integrand_over_z = simps(integrand_simps_over_z, integrand_angles) * integrand_radii

                # Integral over both w and z
                integrand_simps_over_wz = integrand_simps_over_wz.reshape((len(integrand_radii), len(integrand_angles), len(integrand_radii), len(integrand_angles)))
                integral_over_angles = simps(integrand_simps_over_wz, integrand_angles) * integrand_radii
                integral_over_w = simps(simps(integrand_simps_over_wz, integrand_angles) * integrand_radii,integrand_radii)
                integrand_over_wz = simps(integral_over_w, integrand_angles) * integrand_radii
                plt.plot(integrand_radii, integrand_over_z, label='NLO - integrand over z')
                plt.plot(integrand_radii, integrand_over_wz, label='NLO - integrand over wz integrated over z')
                # plt.plot(integrand_radii, integral_over_angles[1, 1, :], label='NLO - integrand over wz at one point')
                plt.legend()
            else:
                integrand_simps = integrand(calculation)
                integrand_simps = integrand_simps.reshape((len(integrand_radii), len(integrand_angles)))
                integrand_simps = (integrand_simps.T * integrand_radii).T
                integrated_over_angles = simps(integrand_simps, integrand_angles)

                # distribution_importance_sampling = integrand_radii * integrand_radii
                distribution_importance_sampling = 1.

                # plt.plot(integrand_radii, distribution_importance_sampling * np.abs(integrated_over_angles))
                plt.plot(integrand_radii, distribution_importance_sampling * integrated_over_angles)
            # y-axis log scale
            # plt.yscale('log')
            plt.xscale('log')
            plt.title('r = ' + str(r) + ', b = ' + str(b))
            plt.show()
            plt.close()

            # Dont worry about angles. Worry only about r-dependence. Do the same for NLO with two r-dependences
            # TODO: Better importance sampling would be a gaussian from 10**-3 to 2.*10**1 - in log-scale
            # TODO: That means from -3 to ~1.5 in log-scale


