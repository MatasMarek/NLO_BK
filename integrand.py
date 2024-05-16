from kernel import Kb, Kf, Krc, Kci, K1
# from kernel import Ka
import numpy as np


def integrand(calculation):
    # This function is used to calculate the integrand in the BK equation.
    if calculation['order_of_BK'] == 'LO':
        return integrand_LO(calculation)
    elif calculation['order_of_BK'] == 'NLO' or calculation['order_of_BK'] == 'NLO_LOcut':
        return integrand_NLO(calculation)
    elif calculation['order_of_BK'] == 'ci' or calculation['order_of_BK'] == 'ci_nokdla':
        return integrand_ci(calculation)
    return


def integrand_LO(calculation):
    Ns = calculation['Ns']
    if calculation['order_of_rk'] == 4:
        total = Krc(calculation) * (Ns['zx'] + Ns['zy'] - Ns['xy'] - Ns['zx'] * Ns['zy'])
        kernel = Krc(calculation)
        split = Krc(calculation) * (Ns['zx'] + Ns['zy'])
        return total, kernel, split
    else:
        return Krc(calculation) * (Ns['zx'] + Ns['zy'] - Ns['xy'] - Ns['zx'] * Ns['zy'])


def integrand_ci(calculation):
    Ns = calculation['Ns']
    if calculation['order_of_rk'] == 4:
        total = Kci(calculation) * (Ns['zx'] + Ns['zy'] - Ns['xy'] - Ns['zx'] * Ns['zy'])
        kernel = Kci(calculation)
        split = Kci(calculation) * (Ns['zx'] + Ns['zy'])
        return total, kernel, split
    else:
        return Kci(calculation) * (Ns['zx'] + Ns['zy'] - Ns['xy'] - Ns['zx'] * Ns['zy'])


def integrand_NLO(calculation):
    # allowed keys are 'zx', 'zy', 'wx', 'wz', 'wy'
    Ns = calculation['Ns']

    first_term = K1(calculation) * (Ns['zx'] + Ns['zy'] - Ns['xy'] - Ns['zx'] * Ns['zy'])
    second_term = Kb(calculation) * (Ns['wy'] + Ns['wz'] - Ns['zy'] - Ns['zx'] * Ns['wz'] - Ns['zx'] * Ns['wy'] - Ns['wz'] * Ns['wy'] + Ns['zx'] * Ns['zy'] + Ns['zx'] * Ns['wz'] * Ns['wy'])
    third_term = Kf(calculation) * (Ns['wx'] - Ns['zx'] - Ns['zy'] * Ns['wx'] + Ns['zx'] * Ns['zy'])

    return first_term, second_term + third_term  # First term is integrated just over z, the other two over w




