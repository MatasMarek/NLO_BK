import numpy as np

def print_out_calculation(calculation):
    print('run_name = ', calculation['run_name'])
    print('order_of_BK = ', calculation['order_of_BK'])
    print('integration_method = ', calculation['integration_method'])
    print('order_of_rk = ', calculation['order_of_rk'])
    print('number_of_cores = ', calculation['number_of_cores'])
    print('dimensionality_of_N = ', calculation['dimensionality_of_N'])
    print('N = ', calculation['N'])
    print('N_to_add = ', calculation['N_to_add'])
    print('integrand_values = ', calculation['integrand_values'])
    print('grid_in_Y = ', calculation['grid']['grid_in_Y'])
    print('grid_in_r = ', calculation['grid']['grid_in_r'])
    print('grid_in_b = ', calculation['grid']['grid_in_b'])
    print('grid_in_integrand_radius = ', calculation['grid']['grid_in_integrand_radius'])
    print('grid_in_integrand_angle = ', calculation['grid']['grid_in_integrand_angle'])
    if 'y_ind' in calculation:
        print('y_ind = ', calculation['y_ind'])

    return


def save_calculation(calculation):
    run_name = calculation['run_name']
    np.save('../output/' + run_name + '/calculation', calculation, allow_pickle=True)


def load_calculation(run_name):
    calculation = np.load('../output/' + run_name + '/calculation.npy', allow_pickle=True).item()
    return calculation


def print_calculation_stats(calculation):
    print('#################################################')
    print('############ CALCULATION STARTED ################')
    print('Dimensionality of N:', calculation['dimensionality_of_N'])
    print('Order of the BK:', calculation['order_of_BK'])
    print('RK order:', calculation['order_of_rk'])
    print('Integration method:', calculation['integration_method'])
    print('Number of cores:', calculation['number_of_cores'])
    print('Number od samples:', calculation['no_of_samples'])
    print('Run name:', calculation['run_name'])
    print('#################################################')
    print('############ RUNTIME ESTIMATES #################')

    size_of_n = calculation['initial_cond'].size
    n_of_samples = calculation['no_of_samples']
    steps_in_Y = len(calculation['grid']['grid_in_Y'])
    order_of_BK = calculation['order_of_BK']
    dimensionality_of_N = calculation['dimensionality_of_N']
    calculation_units = size_of_n * n_of_samples
    time_per_unit_LO_1D = 1.5 / 10**4 / 400
    time_per_unit_ci = 121. / 40 / 40 / 10 / 10**4
    time_per_unit_NLO = 17 / 10**4 / 40 / 40
    if order_of_BK == 'LO':
        time_per_unit = time_per_unit_LO_1D
    elif order_of_BK == 'ci':
        time_per_unit = time_per_unit_ci
    elif order_of_BK == 'NLO':
        time_per_unit = time_per_unit_NLO
    else:
        raise ValueError('Order of BK is not recognized.')
    time_per_rapidity = time_per_unit * calculation_units
    total_time = time_per_rapidity * steps_in_Y

    # convert the time such that it is in the most appropriate units
    if total_time < 60:
        total_time = total_time
        time_unit = 's'
    elif total_time < 60 * 60:
        total_time = total_time / 60
        time_unit = 'min'
    elif total_time < 60 * 60 * 24:
        total_time = total_time / 60 / 60
        time_unit = 'h'
    elif total_time < 60 * 60 * 24 * 365:
        total_time = total_time / 60 / 60 / 24
        time_unit = 'days'
    else:
        total_time = total_time / 60 / 60 / 24 / 365
        time_unit = 'years'
    print('Estimates for one core:')
    print('Time per unit of rapidity', round(time_per_rapidity, 2), 's')
    print('Total time:', round(total_time, 2), time_unit)
    print('#################################################')
    print('#################################################')



