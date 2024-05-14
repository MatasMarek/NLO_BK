import os.path
import shutil
import sys
sys.path.append('../')
from run import run_calculation
import numpy as np
import os

#DEBUG YOU SET THE LOWER BOUND TO BE -5
np.seterr(over='raise')
# EDEBUG

steps_in_integrand_theta = 500
shift = 2.*np.pi/steps_in_integrand_theta/2.  # to avoid double counting and y-axis with z and w.
grid = {
    'grid_in_Y': np.linspace(0., 10., 201),
    'grid_in_r': np.logspace(-4., 2., 40),
    'grid_in_integrand_radius': np.logspace(-7., 2., 50000),
    'grid_in_integrand_angle': np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta),  # not to include 2pi to avoid double counting
}


dimensionality_of_N = 1  # r, b

integration_method = 'MC'

no_of_samples = 10**7

order_of_rk = 1
order_of_BK = 'NLO'
number_of_cores = 1


from initial_conds import MV_1D as cond
initial_cond = cond(grid, Qs0_sq=(19.*0.241)**2, Lambda=0.241, gamma=0.6)
run_name = 'pilot_run_NLO_1D_less_radius_more_angles'


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

