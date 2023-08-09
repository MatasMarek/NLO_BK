import numpy as np
import const as c
from scipy.special import j1


def Ka(calculation):
    # This function is used to calculate the Ka kernel in the BK equation.
    var = calculation['variables_for_N']
    # TODO: Check with Honza that all couplings are at the right scale!
    fraction_first = alpha_s_squared(np.array([var['xy']['rsq']])) * c.Nc**2 / (8. * np.pi**3)
    fraction_second = var['xy']['rsq'] / (var['zx']['rsq'] * var['zy']['rsq'])

    bracket = 67./9. - (np.pi**2 / 3.) - 10. * c.nf / (9.*c.Nc) - 2. * np.log(var['zx']['rsq'] / var['xy']['rsq']) * np.log(var['zy']['rsq'] / var['xy']['rsq'])

    kernel = Krc(calculation) + fraction_first * fraction_second * bracket

    #  CONVERGENCE CONDITION
    # Dasa's cutoff of the infinities
    if np.isinf(kernel).any():
        kernel[np.isinf(kernel)] = 0.
    if np.isnan(kernel).any():
        kernel[np.isnan(kernel)] = 0.
    # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    kernel[same_rxy_rzx_or_rzy] = 0.
    return kernel


def Krc(calculation):
    var = calculation['variables_for_N']
    # TODO: Is the power of this alpha correct?
    fraction = np.sqrt(alpha_s(np.array([var['xy']['rsq']]))) * c.Nc / (2. * np.pi**2)
    first = var['xy']['rsq']/(var['zx']['rsq']*var['zy']['rsq'])
    second = (1./var['zx']['rsq']) * (alpha_s(var['zx']['rsq'])/alpha_s(var['zy']['rsq']) - 1.)
    third =  (1./var['zy']['rsq']) * (alpha_s(var['zy']['rsq'])/alpha_s(var['zx']['rsq']) - 1.)
    kernel = fraction * (first + second + third)
    if abs(kernel).any() >= 10**5:
        print('diverging kernel')
        exit(1)
    return kernel


def Kb(calculation):
    var = calculation['variables_for_N']
    fraction = alpha_s_squared(np.array([var['xy']['rsq']])) * c.Nc**2/(8. * np.pi**4)
    in_front_of_bracket = -2./var['wz']['rsq']**2
    bracket_first_up = var['zx']['rsq']*var['wy']['rsq'] + var['wx']['rsq']*var['zy']['rsq'] - 4.*var['xy']['rsq']*var['wz']['rsq']

    bracket_first_down = var['wz']['rsq']**2 * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'])
    bracket_first = bracket_first_up/bracket_first_down
    # TODO Ask Honza; This Kernel diverges for z and w on the y axis. Then the zx = zy and wx = wy
    #  (Only for the 2D case) How to fix this? I have rotated the integrand to miss those points.
    # TODO This goes wrong when zx = zy and wx = wy (both z and w are on the y axis - which I do not see; is this below Numpy's precision?)
    # TODO: Also when zx = wx and zy = wy

    bracket_second = var['xy']['rsq']**2 / (var['zx']['rsq']*var['wy']['rsq'] * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq']))
    bracket_third = var['xy']['rsq'] / (var['zx']['rsq']*var['wy']['rsq']*var['wz']['rsq'])
    bracket = bracket_first + bracket_second + bracket_third
    log_term = np.log(var['zx']['rsq'] * var['wy']['rsq'] / (var['wx']['rsq'] * var['zy']['rsq']))
    kernel = fraction * (in_front_of_bracket + bracket * log_term)

    #  CONVERGENCE CONDITION
    # Dasa's cutoff of the infinities
    if np.isinf(kernel).any():
        kernel[np.isinf(kernel)] = 0.
    if np.isnan(kernel).any():
        kernel[np.isnan(kernel)] = 0.
    # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    kernel[same_rxy_rzx_or_rzy] = 0.
    return kernel


def Kf(calculation):
    var = calculation['variables_for_N']
    fraction = alpha_s_squared(np.array([var['xy']['rsq']])) * c.Nc * c.nf/(8. * np.pi**4)
    first = 2./var['wz']['rsq']**2
    second_up = var['wx']['rsq']*var['zy']['rsq'] + var['wy']['rsq']*var['zx']['rsq'] - var['xy']['rsq']*var['wz']['rsq']
    second_down = var['wz']['rsq']**2 * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'])

    second = second_up/second_down
    log_term = np.log(var['zx']['rsq'] * var['wy']['rsq'] / (var['wx']['rsq'] * var['zy']['rsq']))
    kernel = fraction*(first - second * log_term)

    #  CONVERGENCE CONDITION
    # Dasa's cutoff of the infinities
    if np.isinf(kernel).any():
        kernel[np.isinf(kernel)] = 0.
    if np.isnan(kernel).any():
        kernel[np.isnan(kernel)] = 0.

    # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    kernel[same_rxy_rzx_or_rzy] = 0.
    return kernel


def Kci(calculation):
    var = calculation['variables_for_N']

    A1 = 11./12.
    Lrxz_rxy = np.log(var['zx']['rsq']/var['xy']['rsq'])
    Lrzy_rxy = np.log(var['zy']['rsq']/var['xy']['rsq'])
    rho_sq = np.abs(Lrxz_rxy * Lrzy_rxy)
    alpha_bar = c.Nc/np.pi * alpha_s(get_min_rs(var))

    factor = alpha_bar/(2.*np.pi)
    prvni = var['xy']['rsq']/(var['zx']['rsq']*var['zy']['rsq'])
    r_zy_zx = np.c_[var['zy']['rsq'], var['zx']['rsq']]
    min_zy_zx = np.min(r_zy_zx, axis=1)

    sign = np.zeros(len(min_zy_zx)) - 1.
    sign[var['xy']['rsq'] < min_zy_zx] = 1.

    bracket = var['xy']['rsq']/min_zy_zx
    exponent = A1 * alpha_bar * sign
    druhy = np.power(bracket, exponent)

    # TODO: Ask Honza; Dasa sometimes uses I1 instead of J1. Should I do that too?

    #  CONVERGENCE CONDITION
    # Setting KDLA to one when rho_sq is zero
    kdla_array = np.zeros(len(rho_sq)) + 1.
    kdla_array[rho_sq > 0.] = kdla(alpha_bar[rho_sq > 0.], rho_sq[rho_sq > 0.])

    Kci = factor * prvni * druhy * kdla_array

    #  CONVERGENCE CONDITION
    # Setting the whole kernel to zero when r_zy or r_zx is zero
    zero_rzx_rzy = min_zy_zx < 10**-10
    Kci[zero_rzx_rzy] = 0.

    # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    Kci[same_rxy_rzx_or_rzy] = 0.
    if abs(Kci).any() >= 10**7 or np.isnan(Kci).any():
        print('diverging kernel')
        exit(1)
    return Kci


def kdla(alpha_bar, rho_sq):
    return j1(2.*np.sqrt(alpha_bar * rho_sq))/(np.sqrt(alpha_bar * rho_sq))


def get_alpha_s_squared_min_r(var):
    min_rs = get_min_rs(var)
    return alpha_s_squared(min_rs)


def get_min_rs(var):
    # all_rs = np.c_[var['xy']['rsq'] + np.zeros(len(var['zy']['rsq'])), var['zy']['rsq'], var['zx']['rsq'], var['wy']['rsq'], var['wx']['rsq'], var['wz']['rsq']]
    all_rs = np.c_[var['xy']['rsq'] + np.zeros(len(var['zy']['rsq'])), var['zy']['rsq'], var['zx']['rsq']]
    min_rs = np.min(all_rs, axis=1)
    return min_rs


def alpha_s_squared(values):
    return alpha_s(values)**2


def beta(nf):
    return 11. - nf * 2. / 3.


def alpha_s(values):
    # TODO: Ask Honza; Should the alphas be bared or not? Do I multiply with Nc/pi?
    alpha = np.zeros(len(values))
    alpha[values >= c.rsat_sq] = c.alpha_freeze
    condition_lambda_3 = np.logical_and(c.rsat_sq > values, values >= c.rc_sq)
    alpha[condition_lambda_3] = alpha_final_term(values[condition_lambda_3], 3., c.lambdaQCD3)
    condition_lambda_4 = np.logical_and(c.rc_sq > values, values >= c.rb_sq)
    alpha[condition_lambda_4] = alpha_final_term(values[condition_lambda_4], 4., c.lambdaQCD4)
    alpha[c.rb_sq > values] = alpha_final_term(values[c.rb_sq > values], 5., c.lambdaQCD5)
    return alpha


def alpha_final_term(values, nf, lambda_qcd):
    return 4. * np.pi / ((11. - nf * 2. / 3.) * np.log(4. * c.C ** 2. / (values * lambda_qcd ** 2)))





