import shutil
from run import run_calculation
import numpy as np
import os

# DEBUG YOU SET THE LOWER BOUND TO BE -5
np.seterr(over='raise')
# EDEBUG

steps_in_integrand_theta = 200
shift = 2.*np.pi/steps_in_integrand_theta/2.  # to avoid double counting and y-axis with z and w.

steps_in_grid_theta = 10
theta_shift = 2.*np.pi/steps_in_integrand_theta/2.  # to avoid double counting and y-axis with z and w.

steps_in_grid_phi = 16
phi_shift = 2.*np.pi/steps_in_grid_phi/2.  # to avoid double counting and y-axis with z and w.

grid = {
    'grid_in_Y': np.linspace(0., 10., 201),
    'grid_in_r': np.logspace(-6., 2., 40),
    'grid_in_b': np.logspace(-1., 2., 20),
    'grid_in_theta': np.linspace(-np.pi + theta_shift, np.pi - theta_shift, steps_in_grid_theta),  # not to include 2pi to avoid double counting
    'grid_in_phi': np.linspace(-np.pi + phi_shift, np.pi - phi_shift, steps_in_grid_phi),  # not to include 2pi to avoid double counting
    'grid_in_integrand_radius': np.logspace(-7., 2., 15000),
    'grid_in_integrand_angle': np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta),  # not to include 2pi to avoid double counting
}

dimensionality_of_N = 4  # r, b, theta, phi

integration_method = 'MC'
# integration_method = 'Simps'

no_of_samples = 10**4

order_of_rk = 1
order_of_BK = 'NLO'
number_of_cores = 3

from initial_conds import mareks_N_4D as cond
initial_cond = cond(grid)
run_name = 'NLO_4D_smaller_initial_dipole'


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

if not os.path.isdir('../output/' + run_name):
    os.mkdir('../output/' + run_name)

shutil.copyfile(__file__, '../output/' + run_name + '/input.py')
shutil.copyfile('../const.py', '../output/' + run_name + '/const.py')

run_calculation(calculation)

