import numpy as np
from io_bk import save_calculation, print_calculation_stats, get_time_in_appropriate_units
from runge_kutta import runge_kutta
import time
from joblib import Parallel, delayed
import sys


def get_grids_of_N(calculation):
    list_of_grids_of_N = []
    grid_in_r = calculation['grid']['grid_in_r']
    list_of_grids_of_N.append(grid_in_r)
    if calculation['dimensionality_of_N'] >= 2:
        grid_in_b = calculation['grid']['grid_in_b']
        list_of_grids_of_N.append(grid_in_b)
    if calculation['dimensionality_of_N'] >= 3:
        grid_in_theta = calculation['grid']['grid_in_theta']
        list_of_grids_of_N.append(grid_in_theta)
    if calculation['dimensionality_of_N'] == 4:
        grid_in_phi = calculation['grid']['grid_in_phi']
        list_of_grids_of_N.append(grid_in_phi)
    return list_of_grids_of_N


def set_up_grids(calculation):
    # Get the input grids
    grid_in_y = calculation['grid']['grid_in_Y']
    list_of_grid_dimensions = [len(grid_in_y)]
    list_of_grids_of_N = get_grids_of_N(calculation)

    for grid_of_N in list_of_grids_of_N:
        list_of_grid_dimensions.append(len(grid_of_N))

    #  set up N and N_to_add
    tuple_of_grid_dimensions = tuple(list_of_grid_dimensions)
    N = np.zeros(tuple_of_grid_dimensions)
    calculation['N'] = N
    tuple_of_n_to_add_dimensions = tuple(list_of_grid_dimensions[1:])
    N_to_add = np.zeros(tuple_of_n_to_add_dimensions)
    calculation['N_to_add'] = N_to_add

    #  set up integrand coods - cartesian coordinates of each point in integrand from its polar coordinates
    integrand_angles = calculation['grid']['grid_in_integrand_angle']
    integrand_radii = calculation['grid']['grid_in_integrand_radius']
    integrand_cartesian_coods = np.zeros((len(integrand_angles) * len(integrand_radii), 2))
    for radius_ind in range(len(integrand_radii)):
        for theta_ind in range(len(integrand_angles)):
            x_cood = integrand_radii[radius_ind] * np.cos(integrand_angles[theta_ind])
            y_cood = integrand_radii[radius_ind] * np.sin(integrand_angles[theta_ind])
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][0] = x_cood
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][1] = y_cood

    calculation['integrand_cartesian_coods'] = integrand_cartesian_coods
    return calculation


def make_step_in_r(r_ind, calculation):
    y_ind = calculation['y_ind']
    r = calculation['grid']['grid_in_r'][r_ind]
    calculation['variables_for_N']['xy']['r'] = r
    calculation['variables_for_N']['xy']['rsq'] = r ** 2
    N_to_add_r_ind = calculation['N_to_add'][r_ind]
    if calculation['dimensionality_of_N'] == 1:
        x = np.array([-r / 2., 0.])
        y = np.array([r / 2., 0.])
        calculation['Ns']['xy'] = calculation['N'][y_ind - 1][r_ind]
        N_to_add_r_ind = runge_kutta(calculation, x, y)
    else:
        grid_in_b = calculation['grid']['grid_in_b']
        for b_ind in range(len(grid_in_b)):
            b = grid_in_b[b_ind]
            calculation['variables_for_N']['xy']['b'] = b
            if calculation['dimensionality_of_N'] == 2:
                calculation['Ns']['xy'] = calculation['N'][y_ind - 1][r_ind][b_ind]
                x = np.array([-r / 2., b])
                y = np.array([r / 2., b])
                N_to_add_r_ind[b_ind] = runge_kutta(calculation, x, y)
            else:
                grid_in_theta = calculation['grid']['grid_in_theta']
                for theta_ind in range(len(grid_in_theta)):
                    theta = grid_in_theta[theta_ind]
                    calculation['variables_for_N']['xy']['theta'] = theta
                    if calculation['dimensionality_of_N'] == 3:
                        calculation['Ns']['xy'] = calculation['N'][y_ind - 1][r_ind][b_ind][theta_ind]
                        x = np.array([-r / 2. * np.cos(theta), b - r / 2. * np.sin(theta)])
                        y = np.array([r / 2. * np.cos(theta), b + r / 2. * np.sin(theta)])
                        N_to_add_r_ind[b_ind][theta_ind] = runge_kutta(calculation, x, y)
                    else:
                        grid_in_phi = calculation['grid']['grid_in_phi']
                        for phi_ind in range(len(grid_in_phi)):
                            phi = grid_in_phi[phi_ind]
                            calculation['variables_for_N']['xy']['phi'] = phi
                            calculation['Ns']['xy'] = calculation['N'][y_ind - 1][r_ind][b_ind][theta_ind][phi_ind]
                            x = np.array([-r / 2. * np.cos(theta) + b * np.cos(phi),
                                        -r / 2. * np.sin(theta) + b * np.sin(phi)])
                            y = np.array([r / 2. * np.cos(theta) + b * np.cos(phi),
                                           r / 2. * np.sin(theta) + b * np.sin(phi)])
                            N_to_add_r_ind[b_ind][theta_ind][phi_ind] = runge_kutta(calculation, x, y)

    # if calculation['grid']['grid_in_r'][r_ind] >= 10**-2:
    #     if abs(N_to_add_r_ind[0] / N_to_add_r_ind[1] - 1) > 10**-4 or abs(N_to_add_r_ind[0] / N_to_add_r_ind[2] - 1) > 10**-4:
    #         print('Here be trouble r_ind', r_ind, 'N_to_add_r_ind', N_to_add_r_ind)
    return N_to_add_r_ind


def make_step(calculation):
    y_ind = calculation['y_ind']
    grid_in_r = calculation['grid']['grid_in_r']
    calculation['variables_for_N'] = {'xy': {}}
    calculation['Ns'] = {}

    number_of_cores = calculation['number_of_cores']
    if number_of_cores > 1:
        N_to_add = Parallel(n_jobs=number_of_cores)(delayed(make_step_in_r)(r_ind, calculation) for r_ind in range(len(grid_in_r)))
        calculation['N_to_add'] = np.array(N_to_add).squeeze()
    else:
        for r_ind in range(len(grid_in_r)):
            N_to_add_r_ind = make_step_in_r(r_ind, calculation)
            calculation['N_to_add'][r_ind] = N_to_add_r_ind

    # CONVERGENCE CONDITION - if N_to_add jumps too much, I zero it
    condition_one = np.abs(calculation['N_to_add']/calculation['N'][y_ind]) > 100.
    condition_two = calculation['N'][y_ind] > 0.
    condition_three = np.abs(calculation['N_to_add']) > 0.04
    # combine the conditions
    N_to_add_to_zero_out = (condition_one & condition_two) | condition_three
    calculation['N_to_add'][N_to_add_to_zero_out] = 0.
    # log it
    print('max(N_to_add)', np.max(calculation['N_to_add']), 'max(N)', np.max(calculation['N'][y_ind]), 'max(N/N_to_add)', np.max(calculation['N_to_add']))
    if np.sum(N_to_add_to_zero_out) > 0:
        print('N_to_add jumps too much at y_ind', y_ind, 'r_inds', np.where(N_to_add_to_zero_out)[0], 'for', np.sum(N_to_add_to_zero_out), 'points')
    # CONVERGENCE CONDITION

    calculation['N'][y_ind] = calculation['N'][y_ind - 1] + calculation['N_to_add']

    #  CONVERGENCE CONDITION
    # CHECK FOR NEGATIVE VALUES IN N
    if (calculation['N'][y_ind] < 0.).any():
        print(100.*np.sum(calculation['N'][y_ind] < 0.)/calculation['N'][y_ind].size, '% of negative values in N with a min of', np.min(calculation['N'][y_ind]))
        calculation['N'][y_ind][calculation['N'][y_ind] < 0.] = 0.
    if (calculation['N'][y_ind] > 1.).any():
        print(100.*np.sum(calculation['N'][y_ind] > 1.)/calculation['N'][y_ind].size, ' % of values above 1 in N with a max of', np.max(calculation['N'][y_ind]))
        calculation['N'][y_ind][calculation['N'][y_ind] > 1.] = 1.
    # CHECK FOR NEGATIVE VALUES IN N
    return calculation


def run_calculation(calculation):
    print_calculation_stats(calculation)
    stdoutOrigin = sys.stdout
    output_path = '../output/' + calculation['run_name'] + '/log.txt'
    sys.stdout = open(output_path, "w")
    print_calculation_stats(calculation)

    calculation = set_up_grids(calculation)
    start_time = time.time()
    calculation['N'][0] = calculation['initial_cond']

    sys.stdout = stdoutOrigin
    print('Y', 0.00, "--- %s seconds ---" % round(time.time() - start_time, 1))
    sys.stdout = open(output_path, "a")
    print('Y', 0.00, "--- %s seconds ---" % round(time.time() - start_time, 1))

    for y_ind in range(1, len(calculation['grid']['grid_in_Y'])):
        sys.stdout.flush()
        calculation['y_ind'] = y_ind
        if (y_ind-1) % 10 == 0. or calculation['order_of_BK'] == 'NLO' or calculation['dimensionality_of_N'] >= 2 or calculation['order_of_BK'] == 'NLO_LOcut':
            save_calculation(calculation)
        calculation = make_step(calculation)
        sys.stdout = stdoutOrigin
        run_time, time_unit = get_time_in_appropriate_units(time.time() - start_time)
        print('Y', round(calculation['grid']['grid_in_Y'][y_ind], 2), "--- " + str(run_time) + time_unit + " ---")
        sys.stdout = open(output_path, "a")
        print('Y', round(calculation['grid']['grid_in_Y'][y_ind], 2), "--- " + str(run_time) + time_unit + " ---")
    save_calculation(calculation)
    return

# TODO: Porozdelovat kod a ucesat ho

# TODO Conditions:
#  * r1 not the same size as r
#  * fraction inside kernel not infinite - done for ci, not for LO
#  * It might be that the interpolation method should extrapolate rather than cut off
#  * It might be that I ask the interpolation method for really small rs that I do not have in the grid
#  * Check with Guillermo's conditions

# TODO:
#  * Get some convergence estimates for 2D ciBK

# It might very well be, that the unwanted unstable behavior for NLO comes from points when w and z are too close. Kernels diverge then.
# TODO: I DID NOT HAVE ENOUGH POINTS IN INTEGRATION OVER r!, Change that in all inputs and add a check for that
# TODO: Check the LO and NLO integrands and come up with a better importance sampling for the two cases

