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

