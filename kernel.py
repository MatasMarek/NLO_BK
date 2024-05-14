import numpy as np
import const as c
from scipy.special import i1, j1


# def Ka(calculation):
#     # This function is used to calculate the Ka kernel in the BK equation.
#     var = calculation['variables_for_N']
#     # TODO: Check with Honza that all couplings are at the right scale!
#     fraction_first = alpha_s_squared(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.Nc**2 / (8. * np.pi**3)
#     fraction_second = var['xy']['rsq'] / (var['zx']['rsq'] * var['zy']['rsq'])
#     bracket = 67. / 9. - (np.pi ** 2 / 3.) - 10. * c.nf / (9. * c.Nc) - 2. * np.log(var['zx']['rsq'] / var['xy']['rsq']) * np.log(var['zy']['rsq'] / var['xy']['rsq'])
#
#     if calculation['order_of_BK'] == 'NLO_LOcut':
#         kernel = Krc(calculation)
#     else:
#         kernel = Krc(calculation) + fraction_first * fraction_second * bracket
#
#     #  CONVERGENCE CONDITION
#     # Dasa's cutoff of the infinities
#     if np.isinf(kernel).any():
#         kernel[np.isinf(kernel)] = 0.
#     if np.isnan(kernel).any():
#         kernel[np.isnan(kernel)] = 0.
#     # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
#     # same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
#     # same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
#     # same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
#     # kernel[same_rxy_rzx_or_rzy] = 0.
#     return kernel

def K1(calculation):
    var = calculation['variables_for_N']

    Lrxz_rxy = np.log(var['zx']['rsq']/var['xy']['rsq'])
    Lrzy_rxy = np.log(var['zy']['rsq']/var['xy']['rsq'])
    rho_sq = Lrxz_rxy * Lrzy_rxy
    alpha_bar = np.zeros(len(rho_sq)) + c.Nc/np.pi * alpha_s(np.array([var['xy']['rsq']]), alpha_fixed=True)

    kdla_array = kdla_with_limits(alpha_bar, rho_sq)
    kernel = Krc(calculation, alpha_fixed=True)*Kstl(calculation)*kdla_array - Ksub(calculation) + Kfin(calculation)

    # CONVERGENCE CONDITION - ZERO OUT INFINITELY SMALL rs
    zero_terms = (var['zx']['rsq'] ==  0.) | (var['zy']['rsq'] == 0.) | (var['wx']['rsq'] == 0.) | (var['wy']['rsq'] == 0.)
    kernel[zero_terms] = 0.
    return kernel

def Kfin(calculation):
    # This function is used to calculate the Kfin kernel in the BK equation.
    var = calculation['variables_for_N']
    fraction_first = alpha_s_bar_squared(np.array([var['xy']['rsq']]), alpha_fixed=True) / (8. * np.pi)
    fraction_second = var['xy']['rsq'] / (var['zx']['rsq'] * var['zy']['rsq'])
    bracket = 67. / 9. - (np.pi ** 2 / 3.) - 10. * c.nf / (9. * c.Nc)
    kernel = fraction_first * fraction_second * bracket
    return kernel


def Ksub(calculation):
    # This function is used to calculate the Kfin kernel in the BK equation.
    var = calculation['variables_for_N']
    fraction_first = alpha_s_bar(np.array([var['xy']['rsq']]), alpha_fixed=True) / (2. * np.pi)
    fraction_second = var['xy']['rsq'] / (var['zx']['rsq'] * var['zy']['rsq'])

    r_zy_zx = np.c_[var['zy']['rsq'], var['zx']['rsq']]
    min_zy_zx = np.min(r_zy_zx, axis=1)

    bracket = -alpha_s_bar(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.A1*np.abs(np.log((c.C_sub*var['xy']['rsq'])/(min_zy_zx)))
    kernel = fraction_first * fraction_second * bracket
    return kernel


def Kstl(calculation):
    # This function is used to calculate the Kfin kernel in the BK equation.
    var = calculation['variables_for_N']

    r_zy_zx = np.c_[var['zy']['rsq'], var['zx']['rsq']]
    min_zy_zx = np.min(r_zy_zx, axis=1)

    bracket = -alpha_s_bar(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.A1*np.abs(np.log((c.C_sub*var['xy']['rsq'])/(min_zy_zx)))
    kernel = np.exp(bracket)
    return kernel


def Krc(calculation, alpha_fixed=False):
    var = calculation['variables_for_N']
    fraction = alpha_s(np.array([var['xy']['rsq']]), alpha_fixed) * c.Nc / (2. * np.pi**2)
    first = var['xy']['rsq']/(var['zx']['rsq']*var['zy']['rsq'])
    second = (1./var['zx']['rsq']) * (alpha_s(var['zx']['rsq'], alpha_fixed)/alpha_s(var['zy']['rsq'], alpha_fixed) - 1.)
    third =  (1./var['zy']['rsq']) * (alpha_s(var['zy']['rsq'], alpha_fixed)/alpha_s(var['zx']['rsq'], alpha_fixed) - 1.)
    kernel = fraction * (first + second + third)
    # if abs(kernel).any() >= 10**5:
    #     print('diverging kernel')
    #     exit(1)
    return kernel


def replace_problematic_terms(problematic_terms_bool, problematic_values, not_problematic_values):
    final_thing = np.zeros(len(problematic_terms_bool))
    final_thing[problematic_terms_bool] = problematic_values
    final_thing[np.invert(problematic_terms_bool)] = not_problematic_values
    return final_thing

def Kb(calculation):
    var = calculation['variables_for_N']
    fraction = alpha_s_squared(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.Nc**2/(8. * np.pi**4)
    in_front_of_bracket = -2./var['wz']['rsq']**2
    bracket_first_up = var['zx']['rsq']*var['wy']['rsq'] + var['wx']['rsq']*var['zy']['rsq'] - 4.*var['xy']['rsq']*var['wz']['rsq']

    bracket_first_down = var['wz']['rsq']**2 * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'])

    problematic_terms_bool = var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'] < 10**-20
    not_problematic_terms_bool = np.invert(problematic_terms_bool)


    problematic_values = (bracket_first_up/var['wz']['rsq']**2)[problematic_terms_bool]
    not_problematic_values = bracket_first_up[not_problematic_terms_bool]/bracket_first_down[not_problematic_terms_bool]
    bracket_first = replace_problematic_terms(problematic_terms_bool, problematic_values, not_problematic_values)

    # TODO Ask Honza; This Kernel diverges for z and w on the y axis. Then the zx = zy and wx = wy
    #  (Only for the 2D case) How to fix this? I have rotated the integrand to miss those points.
    # TODO This goes wrong when zx = zy and wx = wy (both z and w are on the y axis - which I do not see; is this below Numpy's precision?)
    # TODO: Also when zx = wx and zy = wy

    bracket_second_down = (var['zx']['rsq']*var['wy']['rsq'] * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq']))
    problematic_values = (var['xy']['rsq']**2/(var['zx']['rsq']*var['wy']['rsq']))[problematic_terms_bool]
    not_problematic_values = (var['xy']['rsq']**2)/ bracket_second_down[not_problematic_terms_bool]
    bracket_second = replace_problematic_terms(problematic_terms_bool, problematic_values, not_problematic_values)

    bracket_third = var['xy']['rsq'] / (var['zx']['rsq']*var['wy']['rsq']*var['wz']['rsq'])
    bracket_third[problematic_terms_bool] = 0.
    bracket = bracket_first + bracket_second + bracket_third

    # DEBUG
    log_term = np.ones(len(bracket))
    # log_term = np.zeros(len(bracket))
    log_term[not_problematic_terms_bool] = np.log((var['zx']['rsq'] * var['wy']['rsq'])[not_problematic_terms_bool] / (var['wx']['rsq'] * var['zy']['rsq'])[not_problematic_terms_bool])
    kernel = fraction * (in_front_of_bracket + bracket * log_term)

    # OLD VERSION OF KERNEL
    fraction = alpha_s_squared(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.Nc**2/(8. * np.pi**4)
    in_front_of_bracket = -2./var['wz']['rsq']**2
    bracket_first_up = var['zx']['rsq']*var['wy']['rsq'] + var['wx']['rsq']*var['zy']['rsq'] - 4.*var['xy']['rsq']*var['wz']['rsq']
    bracket_first_down = var['wz']['rsq']**2 * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'])
    bracket_first = bracket_first_up/bracket_first_down
    bracket_second = var['xy']['rsq']**2 / (var['zx']['rsq']*var['wy']['rsq'] * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq']))
    bracket_third = var['xy']['rsq'] / (var['zx']['rsq']*var['wy']['rsq']*var['wz']['rsq'])
    bracket = bracket_first + bracket_second + bracket_third
    log_term = np.log(var['zx']['rsq'] * var['wy']['rsq'] / (var['wx']['rsq'] * var['zy']['rsq']))
    kernel = fraction * (in_front_of_bracket + bracket * log_term)
    # print('Kb Old and new kernel check', (kernel == kernel_new)[not_problematic_terms_bool].all())
    # OLD VERSION OF KERNEL

    #  CONVERGENCE CONDITION
    # Dasa's cutoff of the infinities
    if np.isinf(kernel).any():
        kernel[np.isinf(kernel)] = 0.
        print('Kb Inf')
    if np.isnan(kernel).any():
        kernel[np.isnan(kernel)] = 0.
        print('Kb Nan')
    # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    # same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    # same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    # same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    # kernel[same_rxy_rzx_or_rzy] = 0.

    # CONVERGENCE CONDITION - ZERO OUT INFINITELY SMALL rs
    zero_terms = (var['zx']['rsq'] ==  0.) | (var['zy']['rsq'] == 0.) | (var['wx']['rsq'] == 0.) | (var['wy']['rsq'] == 0.)
    kernel[zero_terms] = 0.

    if calculation['order_of_BK'] == 'NLO_LOcut':
        return np.zeros(len(kernel))
    return kernel


def Kf(calculation):
    var = calculation['variables_for_N']
    fraction = alpha_s_squared(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.Nc * c.nf/(8. * np.pi**4)
    first = 2./var['wz']['rsq']**2
    second_up = var['wx']['rsq']*var['zy']['rsq'] + var['wy']['rsq']*var['zx']['rsq'] - var['xy']['rsq']*var['wz']['rsq']
    second_down = var['wz']['rsq']**2 * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'])

    problematic_terms_bool = var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'] < 10**-20
    not_problematic_terms_bool = np.invert(problematic_terms_bool)

    # PERFORMING THE LIMIT ON THE LOG-TERM
    problematic_values = (second_up/var['wz']['rsq']**2)[problematic_terms_bool]
    not_problematic_values = second_up[not_problematic_terms_bool]/second_down[not_problematic_terms_bool]
    second = replace_problematic_terms(problematic_terms_bool, problematic_values, not_problematic_values)

    # DEBUG
    log_term = np.ones(len(second))
    # log_term = np.zeros(len(second))
    log_term[not_problematic_terms_bool] = np.log((var['zx']['rsq'] * var['wy']['rsq'])[not_problematic_terms_bool] / (var['wx']['rsq'] * var['zy']['rsq'])[not_problematic_terms_bool])

    kernel = fraction*(first - second * log_term)

    # OLD KERNEL
    fraction = alpha_s_squared(np.array([var['xy']['rsq']]), alpha_fixed=True) * c.Nc * c.nf/(8. * np.pi**4)
    first = 2./var['wz']['rsq']**2
    second_up = var['wx']['rsq']*var['zy']['rsq'] + var['wy']['rsq']*var['zx']['rsq'] - var['xy']['rsq']*var['wz']['rsq']
    second_down = var['wz']['rsq']**2 * (var['zx']['rsq']*var['wy']['rsq'] - var['wx']['rsq']*var['zy']['rsq'])
    second = second_up/second_down
    log_term = np.log(var['zx']['rsq'] * var['wy']['rsq'] / (var['wx']['rsq'] * var['zy']['rsq']))
    kernel = fraction*(first - second * log_term)
    # print('Kf Old and new kernel check', (kernel == kernel_new)[not_problematic_terms_bool].all())
    # OLD KERNEL

    #  CONVERGENCE CONDITION
    # Dasa's cutoff of the infinities
    if np.isinf(kernel).any():
        kernel[np.isinf(kernel)] = 0.
        print('Kf Inf')
    if np.isnan(kernel).any():
        kernel[np.isnan(kernel)] = 0.
        print('Kf Nan')

    # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    # same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    # same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    # same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    # kernel[same_rxy_rzx_or_rzy] = 0.
    if calculation['order_of_BK'] == 'NLO_LOcut':
        return np.zeros(len(kernel))

    zero_terms = (var['zx']['rsq'] ==  0.) | (var['zy']['rsq'] == 0.) | (var['wx']['rsq'] == 0.) | (var['wy']['rsq'] == 0.)
    kernel[zero_terms] = 0.

    return kernel


def Kci(calculation):
    var = calculation['variables_for_N']

    Lrxz_rxy = np.log(var['zx']['rsq']/var['xy']['rsq'])
    Lrzy_rxy = np.log(var['zy']['rsq']/var['xy']['rsq'])
    rho_sq = Lrxz_rxy * Lrzy_rxy
    alpha_bar = c.Nc/np.pi * alpha_s_run(get_min_rs(var))

    factor = alpha_bar/(2.*np.pi)
    prvni = var['xy']['rsq']/(var['zx']['rsq']*var['zy']['rsq'])
    r_zy_zx = np.c_[var['zy']['rsq'], var['zx']['rsq']]
    min_zy_zx = np.min(r_zy_zx, axis=1)

    sign = np.zeros(len(min_zy_zx)) - 1.
    sign[var['xy']['rsq'] < min_zy_zx] = 1.

    bracket = var['xy']['rsq']/min_zy_zx
    exponent = c.A1 * alpha_bar * sign
    druhy = np.power(bracket, exponent)

    if calculation['order_of_BK'] == 'ci_nokdla':
        Kci_values = factor * prvni * druhy
    else:
        kdla_array = kdla_with_limits(alpha_bar, rho_sq)
        Kci_values = factor * prvni * druhy * kdla_array


    #  CONVERGENCE CONDITION
    # Setting the whole kernel to zero when r_zy or r_zx is zero
    # zero_rzx_rzy = min_zy_zx < 10**-10
    # Kci_values[zero_rzx_rzy] = 0.
    #
    # # Setting the whole kernel to zero when r_xy is the same as r_zy or r_zx
    # same_rxy_rzx = np.abs(var['xy']['rsq'] - var['zx']['rsq']) < 10**-10
    # same_rxy_rzy = np.abs(var['xy']['rsq'] - var['zy']['rsq']) < 10**-10
    # same_rxy_rzx_or_rzy = np.logical_or(same_rxy_rzx, same_rxy_rzy)
    # Kci_values[same_rxy_rzx_or_rzy] = 0.
    # if abs(Kci_values).any() >= 10**7 or np.isnan(Kci_values).any():
    #     print('diverging kernel')
    #     exit(1)
    return Kci_values


def kdla_with_limits(alpha_bar, rho_sq):
    # Setting KDLA to one when rho_sq is zero
    kdla_array = np.ones(len(rho_sq))
    cutoff = 10 ** -8
    kdla_array[rho_sq > cutoff] = kdla(alpha_bar[rho_sq > cutoff], rho_sq[rho_sq > cutoff], j_bessel=True)
    kdla_array[rho_sq < -cutoff] = kdla(alpha_bar[rho_sq < -cutoff], rho_sq[rho_sq < -cutoff], j_bessel=False)
    return kdla_array

def kdla(alpha_bar, rho_sq, j_bessel):
    if j_bessel:
        return j1(2.*np.sqrt(alpha_bar * rho_sq))/(np.sqrt(alpha_bar * rho_sq))
    else:
        return i1(2.*np.sqrt(-alpha_bar * rho_sq))/(np.sqrt(-alpha_bar * rho_sq))


def get_alpha_s_squared_min_r(var, alpha_fixed):
    min_rs = get_min_rs(var)
    return alpha_s_squared(min_rs, alpha_fixed)


def get_min_rs(var):
    # all_rs = np.c_[var['xy']['rsq'] + np.zeros(len(var['zy']['rsq'])), var['zy']['rsq'], var['zx']['rsq'], var['wy']['rsq'], var['wx']['rsq'], var['wz']['rsq']]
    all_rs = np.c_[var['xy']['rsq'] + np.zeros(len(var['zy']['rsq'])), var['zy']['rsq'], var['zx']['rsq']]
    min_rs = np.min(all_rs, axis=1)
    return min_rs


def alpha_s_squared(values, alpha_fixed):
    return alpha_s(values, alpha_fixed)**2


def alpha_s_bar_squared(values, alpha_fixed):
    return alpha_s(values, alpha_fixed)**2 * c.Nc**2/np.pi**2

def alpha_s_bar(values, alpha_fixed):
    return alpha_s(values, alpha_fixed) * c.Nc/np.pi

def beta(nf):
    return 11. - nf * 2. / 3.


def alpha_s_run(values):
    alpha = np.zeros(len(values))
    alpha[values >= c.rsat_sq] = c.alpha_freeze
    condition_lambda_3 = np.logical_and(c.rsat_sq > values, values >= c.rc_sq)
    alpha[condition_lambda_3] = alpha_final_term(values[condition_lambda_3], 3., c.lambdaQCD3)
    condition_lambda_4 = np.logical_and(c.rc_sq > values, values >= c.rb_sq)
    alpha[condition_lambda_4] = alpha_final_term(values[condition_lambda_4], 4., c.lambdaQCD4)
    alpha[c.rb_sq > values] = alpha_final_term(values[c.rb_sq > values], 5., c.lambdaQCD5)
    return alpha

def alpha_s(values, alpha_fixed):
    if alpha_fixed:
        inside = (c.mu_over_lam**2)**(1./c.c_coupl) + (c.four_e_gamma/(values * c.Lam_fixed**2))**(1./c.c_coupl)
        alpha_bottom =  (11. - 3. * 2. / 3.) * np.log(inside**c.c_coupl)
        return 4. * np.pi / alpha_bottom
    else:
        return alpha_s_run(values)

def alpha_final_term(values, nf, lambda_qcd):
    return 4. * np.pi / ((11. - nf * 2. / 3.) * np.log(4. * c.C ** 2. / (values * lambda_qcd ** 2)))





