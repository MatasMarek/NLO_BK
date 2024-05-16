from kernel import Kb, K1, Kf
from kernel import alpha_s
import numpy as np


# alpha is the same, k1 is the same

def test_k():
    calculation = {}
    calculation['order_of_BK'] = 'NLO'
    calculation['variables_for_N'] = {}
    for cood_key in ['xy', 'zx', 'zy', 'wx', 'wy', 'wz']:
        calculation['variables_for_N'][cood_key] = {}

    calculation['variables_for_N']['xy']['r'] = 1.00000000e-06
    calculation['variables_for_N']['zx']['r'] = np.array([4.00000000e-07])
    calculation['variables_for_N']['zy']['r'] = np.array([6.00000000e-07])
    calculation['variables_for_N']['wx']['r'] = np.array([5.*10**-8])
    calculation['variables_for_N']['wy']['r'] = np.array([1.5*10**-7])
    calculation['variables_for_N']['wz']['r'] = np.array([0.])

    for cood_key in ['xy', 'zx', 'zy', 'wx', 'wy', 'wz']:
        calculation['variables_for_N'][cood_key]['rsq'] = calculation['variables_for_N'][cood_key]['r'] ** 2

    k1_results = K1(calculation)



    # calculation['variables_for_N']['xy']['r'] = 1.
    # calculation['variables_for_N']['zx']['r'] = np.array([2.])
    # calculation['variables_for_N']['zy']['r'] = np.array([4.])
    # calculation['variables_for_N']['wx']['r'] = np.array([3.])
    # calculation['variables_for_N']['wy']['r'] = np.array([7.])
    # calculation['variables_for_N']['wz']['r'] = np.array([5.])

    matej_data = np.loadtxt('data/out.csv', delimiter=',')
    for line in matej_data:
        calculation['variables_for_N']['xy']['r'] = line[0]
        calculation['variables_for_N']['zx']['r'] = np.array([line[1]])
        calculation['variables_for_N']['wx']['r'] = np.array([line[2]])
        calculation['variables_for_N']['zy']['r'] = np.array([line[3]])
        calculation['variables_for_N']['wz']['r'] = np.array([line[4]])
        calculation['variables_for_N']['wy']['r'] = np.array([line[5]])





        for cood_key in ['xy', 'zx', 'zy', 'wx', 'wy', 'wz']:
            calculation['variables_for_N'][cood_key]['rsq'] = calculation['variables_for_N'][cood_key]['r'] ** 2

        kb_results = Kb(calculation)
        kf_results = Kf(calculation)
        k1_results = K1(calculation)

        threshold = 10**-16
        if abs(k1_results - line[6]) > threshold and  abs(1.- k1_results/line[6]) > 0.00001:
            print('k1', k1_results/line[6], abs(k1_results - line[6]), k1_results, line)
        if abs(kb_results - line[7]) > threshold and abs(1.- kb_results/line[7]) > 0.00001:
            print('kb', kb_results/line[7], abs(kb_results - line[7]), kb_results, line)
        if abs(kf_results - line[8]) > threshold and abs(1.- kf_results/line[8]) > 0.00001:
            print('kf', kf_results/line[8], abs(kf_results - line[8]), kf_results, line)

        # print(kb_results)
        # print(kf_results)
        # print(k1_results)


# for r in [0.01, 0.1, 1., 10., 100.]:
#     print(alpha_s(r**2, alpha_fixed=True))
test_k()