import numpy as np
from scipy import ndimage
from integrand import integrand
import math
from scipy.integrate import simps


def integrate(calculation, x, y):
    # This function is used to calculate the integral in the BK equation.
    if calculation['integration_method'] == 'MC':
        return integrate_MC(calculation, x, y)

    elif calculation['integration_method'] == 'Simpson':
        return integrate_Simpson(calculation, x, y)

    else:
        print('Integration method not implemented.')
        exit(1)


def get_r_from_cartesian_coods(cartesian_coods_one, cartesian_coods_two):
    norm = np.linalg.norm((cartesian_coods_one - cartesian_coods_two), axis=1, keepdims=True)
    return np.squeeze(norm)


def get_b_from_cartesian_coods(cartesian_coods_one, cartesian_coods_two):
    norm = np.linalg.norm((cartesian_coods_one + cartesian_coods_two)/2., axis=1, keepdims=True)
    return np.squeeze(norm)


def get_phi_from_cartesian_coods(cartesian_coods_one, cartesian_coods_two):
    center_of_dipole_vector = (cartesian_coods_one + cartesian_coods_two)/2.
    phi = np.arctan2(center_of_dipole_vector[:, 0], center_of_dipole_vector[:, 1])
    return phi


def get_theta_from_cartesian_coods(cartesian_coods_one, cartesian_coods_two):
    # TODO: Check this routine; Does it work for all cases?
    dipole_vector = (cartesian_coods_one - cartesian_coods_two)
    theta = np.arctan2(dipole_vector[:, 0], dipole_vector[:, 1])
    return theta


def get_cood_combination(x, y, z, w=None):
    cood_combination = {}
    cood_combination['zx'] = [z, x]
    cood_combination['zy'] = [z, y]
    if w is not None:
        cood_combination['wx'] = [w, x]
        cood_combination['wz'] = [w, z]
        cood_combination['wy'] = [w, y]
    return cood_combination


def rs_bs_and_variables_for_N(calculation, x, y, no_of_samples, array_of_z=None):
    integrand_cartesian_coods = calculation['integrand_cartesian_coods']
    random_indexes_z = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    if calculation['integration_method'] == 'MC':
        array_of_z = integrand_cartesian_coods[random_indexes_z]
    calculation['quark_positions'] = {}
    calculation['quark_positions']['x'] = x
    calculation['quark_positions']['y'] = y
    calculation['quark_positions']['z'] = array_of_z
    if calculation['order_of_BK'] == 'NLO':  # get quark positions based on the order of the BK equation
        if calculation['integration_method'] != 'MC':
            print('Simpson for NLO not implemented')
            exit(1)
        random_indexes_w = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # w

        # CONVERGENCE CONDITION: Fix that you dont have w and z the same since the integral then does not work
        same_elements_in_z_and_w = random_indexes_z == random_indexes_w
        while (same_elements_in_z_and_w).any():
            random_indexes_w[same_elements_in_z_and_w] = np.random.randint(0, len(integrand_cartesian_coods), len(random_indexes_w[same_elements_in_z_and_w]))
            same_elements_in_z_and_w = random_indexes_z == random_indexes_w

        array_of_w = integrand_cartesian_coods[random_indexes_w]

        calculation['quark_positions']['w'] = array_of_w
        cood_combination = get_cood_combination(x, y, array_of_z, w=array_of_w)
    else:
        cood_combination = get_cood_combination(x, y, array_of_z)

    for cood_key in cood_combination:  # loop over all combinations of quark positions to get the arguments of all Ns
        cood_one = cood_combination[cood_key][0]
        cood_two = cood_combination[cood_key][1]
        calculation['variables_for_N'][cood_key] = {}
        calculation['variables_for_N'][cood_key]['r'] = get_r_from_cartesian_coods(cood_one, cood_two)
        calculation['variables_for_N'][cood_key]['rsq'] = calculation['variables_for_N'][cood_key]['r']**2

        if calculation['dimensionality_of_N'] >= 2:
            calculation['variables_for_N'][cood_key]['b'] = get_b_from_cartesian_coods(cood_one, cood_two)
        if calculation['dimensionality_of_N'] >= 3:
            calculation['variables_for_N'][cood_key]['theta'] = get_theta_from_cartesian_coods(cood_one, cood_two)
        if calculation['dimensionality_of_N'] == 4:
            calculation['variables_for_N'][cood_key]['phi'] = get_phi_from_cartesian_coods(cood_one, cood_two)
    return calculation


def find_fractional(array, values, log=False):
    if log:
        array = np.log10(array)
        values = np.log10(values + 10**-30)  # Add small number to avoid log(0)
    step = (array[-1] - array[0])/(len(array) - 1)
    ixs = (values - array[0])/step
    return ixs


def interpolate_N(calculation, rs, bs=None, thetas=None, phis=None):
    N = calculation['N'][calculation['y_ind'] - 1]  # N at the previous grid point in y
    index_r = find_fractional(calculation['grid']['grid_in_r'], rs, log=True)
    if bs is not None:
        index_b = find_fractional(calculation['grid']['grid_in_b'], bs, log=True)
        if thetas is not None:
            index_theta = find_fractional(calculation['grid']['grid_in_theta'], thetas, log=False)
            if phis is not None:
                index_phi = find_fractional(calculation['grid']['grid_in_phi'], phis, log=False)
                indexes = np.c_[index_r, index_b, index_theta, index_phi]
            else:
                indexes = np.c_[index_r, index_b, index_theta]
        else:
            indexes = np.c_[index_r, index_b]
    else:
        indexes = np.c_[index_r]
    N_interpolated = ndimage.map_coordinates(N, indexes.T, order=1, mode='nearest')
    return N_interpolated


def get_Ns(calculation):
    rs_bs_and_variables_for_N = calculation['variables_for_N']
    for cood_key in rs_bs_and_variables_for_N:
        if cood_key != 'xy': # This one I already have stored from the loop in make_step
            rs = rs_bs_and_variables_for_N[cood_key]['r']
            bs, thetas, phis = None, None, None
            if calculation['dimensionality_of_N'] >= 2:
                bs = rs_bs_and_variables_for_N[cood_key]['b']
            if calculation['dimensionality_of_N'] >= 3:
                thetas = rs_bs_and_variables_for_N[cood_key]['theta']
            if calculation['dimensionality_of_N'] == 4:
                phis = rs_bs_and_variables_for_N[cood_key]['phi']
            calculation['Ns'][cood_key] = interpolate_N(calculation, rs, bs, thetas, phis)
    return calculation


def jacobian(distances_from_origin):
    return distances_from_origin * np.log(10.) * distances_from_origin  # One r for polar integration and ln(10)*r for log


def integrate_MC(calculation, x, y):
    no_of_samples = calculation['no_of_samples']
    calculation = rs_bs_and_variables_for_N(calculation, x, y, no_of_samples)
    calculation = get_Ns(calculation)


    # Probability distribution normalization
    probability_normalization_polar = 2. * np.pi
    probability_normalization_log = np.log10(calculation['grid']['grid_in_integrand_radius'][-1]) - np.log10(calculation['grid']['grid_in_integrand_radius'][0])

    # Jacobian
    distances_from_origin_z_quark = np.linalg.norm(calculation['quark_positions']['z'], axis=1, keepdims=True)
    jacobian_for_z_integration = jacobian(distances_from_origin_z_quark)
    normalization = probability_normalization_polar * probability_normalization_log
    jacobian_for_z_integration = np.squeeze(jacobian_for_z_integration)

    if calculation['order_of_BK'] == 'NLO':
        distances_from_origin_w_quark = np.linalg.norm(calculation['quark_positions']['w'], axis=1, keepdims=True)
        jacobian_for_wz_integration = jacobian(distances_from_origin_z_quark) * jacobian(distances_from_origin_w_quark)
        jacobian_for_wz_integration = np.squeeze(jacobian_for_wz_integration)

        evaluated_points_for_z_integration, evaluated_points_for_wz_integration = integrand(calculation)
        integral_over_z = normalization * np.sum(jacobian_for_z_integration * evaluated_points_for_z_integration) / no_of_samples
        integral_over_wz = normalization**2 * np.sum(jacobian_for_wz_integration * evaluated_points_for_wz_integration) / no_of_samples
        return integral_over_z + integral_over_wz
    else:
        if calculation['order_of_rk'] == 1:
            evaluated_points = integrand(calculation)
            integral = normalization * np.sum(jacobian_for_z_integration * evaluated_points) / no_of_samples
            return integral
        elif calculation['order_of_rk'] == 4:
            total, kernel, split = integrand(calculation)
            integral_total = normalization * np.sum(jacobian_for_z_integration * total) / no_of_samples
            integral_kernel = normalization * np.sum(jacobian_for_z_integration * kernel) / no_of_samples
            integral_split = normalization * np.sum(jacobian_for_z_integration * split) / no_of_samples
            return integral_total, integral_kernel, integral_split


def integrate_Simpson(calculation, x, y):
    if calculation['order_of_rk'] != 1:
        print('Simpson method implemented only for Runge Kutta method of the order 1')
        exit(1)
    if calculation['order_of_BK'] != 'NLO':
        integrand_angles = calculation['grid']['grid_in_integrand_angle']
        integrand_radii = calculation['grid']['grid_in_integrand_radius']
        integrand_in_radius = []
        for radius_ind in range(len(integrand_radii)):
            z_coods = []
            for theta_ind in range(len(integrand_angles)):
                z_coods.append(calculation['integrand_cartesian_coods'][radius_ind * len(integrand_angles) + theta_ind])
            calculation = rs_bs_and_variables_for_N(calculation, x, y, 0, np.array(z_coods))
            calculation = get_Ns(calculation)

            integrand_in_radius.append(simps(integrand(calculation), integrand_angles))
        return simps(integrand_in_radius * integrand_radii, integrand_radii)  # Jacobian for the polar integral, log is there by scipy.simps
    else:
        print('Simpson integration for NLO not implemented.')
        exit(1)




