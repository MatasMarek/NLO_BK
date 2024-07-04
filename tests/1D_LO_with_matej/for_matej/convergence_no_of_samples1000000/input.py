import shutil
import sys

from joblib import delayed, Parallel

sys.path.append('../')
from run import run_calculation
import numpy as np
import os

def run_parallel_calc(calculation):
    run_calculation(calculation)
    return calculation['run_name']

threads = []
params_adjusted = {
    'samples_in_r': [131, 161, 201, 251, 281],
    'samples_in_integrand_radius': [10000, 50000, 100000],
    'samples_in_integrand_angle': [140, 200, 400, 600, 800],
    'no_of_samples': [10**5, 10**6, 10**7, 10**8]
}

calculations = []
for key in params_adjusted:
    for value in params_adjusted[key]:
        steps_in_integrand_theta = 400
        shift = 2.*np.pi/steps_in_integrand_theta/2.  # to avoid double counting and y-axis with z and w.

        grid = {
            'grid_in_Y': np.linspace(0., 10., 101),
            'grid_in_r': np.logspace(-7., 2., 201),
            'grid_in_integrand_radius': np.logspace(-7., 2., 50000),
            'grid_in_integrand_angle': np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta),  # not to include 2pi to avoid double counting
        }

        dimensionality_of_N = 1  # r, b

        integration_method = 'MC'

        no_of_samples = 10**7

        order_of_rk = 1
        order_of_BK = 'NLO'
        number_of_cores = 1
        drop_double_log = True

        if key == 'samples_in_r':
            grid['grid_in_r'] = np.logspace(-7., 2., value)
        elif key == 'samples_in_integrand_radius':
            grid['grid_in_integrand_radius'] = np.logspace(-7., 2., value)
        elif key == 'samples_in_integrand_angle':
            steps_in_integrand_theta = value
            shift = 2.*np.pi/steps_in_integrand_theta/2.
            grid['grid_in_integrand_angle'] = np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta)
        elif key == 'no_of_samples':
            no_of_samples = value

        from initial_conds import MV_1D_guillermo

        initial_cond = MV_1D_guillermo(grid, Qs0_sq=1., gamma=1.)

        run_name = 'convergence_' + key + str(value)

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
            'drop_double_log': drop_double_log,
        }

        if not os.path.isdir('../output/' + run_name):
            os.mkdir('../output/' + run_name)

        shutil.copyfile(__file__, '../output/' + run_name + '/input.py')
        shutil.copyfile('../const.py', '../output/' + run_name + '/const.py')
        calculations.append(calculation)

result = (Parallel(n_jobs=len(calculations))(delayed(run_parallel_calc)(calculation) for calculation in calculations))
print(result)
print("All threads have finished.")