import numpy as np
from scipy import ndimage
from integrand import integrand
import math
from scipy.integrate import simps


def integrate(calculation, x, y):
    # This function is used to calculate the integral in the BK equation.
    if calculation['integration_method'] == 'MC':
        return integrate_MC(calculation, x, y)

    elif calculation['integration_method'] == 'Simps':
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
    phi = np.arctan2(center_of_dipole_vector[:, 1], center_of_dipole_vector[:, 0])
    return phi


def get_theta_from_cartesian_coods(cartesian_coods_one, cartesian_coods_two):
    # TODO: Check this routine; Does it work for all cases?
    dipole_vector = (cartesian_coods_one - cartesian_coods_two)
    theta = np.arctan2(dipole_vector[:, 1], dipole_vector[:, 0])
    return theta

# MATEJ BACKLOG
# TODO: To get to the real theta eveybody uses, you have to do theta = moje_thete - phi, I can do that in postpr

def get_cood_combination(x, y, z, w=None):
    cood_combination = {}
    cood_combination['zx'] = [z, x]
    cood_combination['zy'] = [z, y]
    if w is not None:
        cood_combination['wx'] = [w, x]
        cood_combination['wz'] = [w, z]
        cood_combination['wy'] = [w, y]
    return cood_combination


def rs_bs_and_variables_for_N(calculation, x, y, no_of_samples):
    integrand_cartesian_coods = calculation['integrand_cartesian_coods']
    random_indexes_z = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    calculation['quark_positions'] = {}

    if calculation['integration_method'] == 'MC':
        array_of_z = integrand_cartesian_coods[random_indexes_z]
    else:
        array_of_z = integrand_cartesian_coods
    calculation['quark_positions']['z'] = array_of_z

    calculation['quark_positions']['x'] = x
    calculation['quark_positions']['y'] = y

    if calculation['order_of_BK'] == 'NLO' or calculation['order_of_BK'] == 'NLO_LOcut':  # get quark positions based on the order of the BK equation
        if calculation['integration_method'] == 'MC':
            random_indexes_w = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # w
            # CONVERGENCE CONDITION: Fix that you dont have w and z the same since the integral then does not work
            same_elements_in_z_and_w = random_indexes_z == random_indexes_w
            while (same_elements_in_z_and_w).any():
                random_indexes_w[same_elements_in_z_and_w] = np.random.randint(0, len(integrand_cartesian_coods), len(random_indexes_w[same_elements_in_z_and_w]))
                same_elements_in_z_and_w = random_indexes_z == random_indexes_w

            # # CONVERGENCE CONDITION 2.0: Fix that w and z are not too close to each other
            # r_wz = get_r_from_cartesian_coods(integrand_cartesian_coods[random_indexes_w], integrand_cartesian_coods[random_indexes_z])
            # elements_w_and_z_too_close = r_wz < 10**-6
            # while (elements_w_and_z_too_close).any():
            #     random_indexes_w[elements_w_and_z_too_close] = np.random.randint(0, len(integrand_cartesian_coods), len(random_indexes_w[elements_w_and_z_too_close]))
            #     r_wz = get_r_from_cartesian_coods(integrand_cartesian_coods[random_indexes_w], integrand_cartesian_coods[random_indexes_z])
            #     elements_w_and_z_too_close = r_wz < 10 ** -6

            array_of_w = integrand_cartesian_coods[random_indexes_w]
        else:
            # Get a combination of a w for each point in z and vice versa
            array_of_z = np.vstack([integrand_cartesian_coods]*len(integrand_cartesian_coods))
            array_of_w = np.repeat(integrand_cartesian_coods, len(integrand_cartesian_coods), axis=0)
            print(len(array_of_w), len(array_of_z))

            calculation['quark_positions']['z'] = array_of_z

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

    if 'eta_rapidity' in calculation:
        if calculation['eta_rapidity']:
            N = calculation['N']
            r_xy = calculation['variables_for_N']['xy']['r']
            etas = calculation['grid']['grid_in_Y'][calculation['y_ind'] - 1] - np.maximum(np.zeros(len(rs)), 2.*np.log(r_xy/rs))
            index_eta = find_fractional(calculation['grid']['grid_in_Y'], etas, log=False)
            indexes = np.c_[index_eta, indexes]
            return ndimage.map_coordinates(N, indexes.T, order=1, mode='nearest')

    N = calculation['N'][calculation['y_ind'] - 1]  # N at the previous grid point in y
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

    if calculation['order_of_BK'] == 'NLO' or calculation['order_of_BK'] == 'NLO_LOcut':
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
    integrand_angles = calculation['grid']['grid_in_integrand_angle']
    integrand_radii = calculation['grid']['grid_in_integrand_radius']
    calculation = rs_bs_and_variables_for_N(calculation, x, y, 0)
    calculation = get_Ns(calculation)

    if calculation['order_of_BK'] == 'NLO' or calculation['order_of_BK'] == 'NLO_LOcut':
        integrand_simps_over_z, integrand_simps_over_wz = integrand(calculation)

        # Integral over z
        integrand_simps_over_z = integrand_simps_over_z[:len(calculation['integrand_cartesian_coods'])]  # Make it only for one integrand of z instead of the repeating tiled version
        integrand_simps_over_z = integrand_simps_over_z.reshape((len(integrand_radii), len(integrand_angles)))
        integral_over_z = simps(simps(integrand_simps_over_z, integrand_angles)*integrand_radii, integrand_radii)

        # Integral over both w and z
        integrand_simps_over_wz = integrand_simps_over_wz.reshape((len(integrand_radii), len(integrand_angles), len(integrand_radii), len(integrand_angles)))
        integral_over_w = simps(simps(integrand_simps_over_wz, integrand_angles)*integrand_radii, integrand_radii)
        integral_over_wz = simps(simps(integral_over_w, integrand_angles)*integrand_radii, integrand_radii)
        return integral_over_z + integral_over_wz

    else:
        integrand_simps = integrand(calculation)
        integrand_simps = integrand_simps.reshape((len(integrand_radii), len(integrand_angles)))
        # Add the Jacobian for the polar integral, log is there by scipy.simps - You added it in the Simps method
        # integrand_simps = (integrand_simps.T * integrand_radii).T
        # Integrate
        integral_value = simps(simps(integrand_simps, integrand_angles) * integrand_radii, integrand_radii)
        # TODO: Test this and the NLO version above
        return integral_value




