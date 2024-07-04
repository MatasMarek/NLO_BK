from kernel import K1, Kci, Kb, Kf
import numpy as np
import matplotlib.pyplot as plt

def plot_k1_and_kci():
    count = 100
    rxzs = np.linspace(0.000000001, 0.0000001, count)
    ryzs = np.linspace(0.1, 0.4, count)

    # stack the arrays such that I can calculate the kernel for all combinations of rxzs and ryzs
    rxzs_tiled = np.tile(rxzs, count)
    ryzs_tiled = np.repeat(ryzs, count)

    rxy = 0.3
    calculation = {
        'variables_for_N': {
            'zx': {'r': rxzs_tiled, 'rsq': rxzs_tiled**2},
            'zy': {'r': ryzs_tiled, 'rsq': ryzs_tiled**2},
            'xy': {'r': rxy, 'rsq': rxy**2},
        },
        'order_of_BK': 'ci',

    }
    k1 = np.log10(np.abs(K1(calculation)))
    kci = np.log10(np.abs(Kci(calculation)))
    k1_kci = np.log10(np.abs(K1(calculation) / Kci(calculation)))
    # k1 = K1(calculation)
    # kci = Kci(calculation)
    # Plot the two kernels in 2D
    for k in ['k1', 'kci', 'k1/kci', 'kb', 'kf']:
    # for k in ['k1/kci']:
        plt.figure()
        if k == 'k1':
            plt.scatter(rxzs_tiled, ryzs_tiled, c=k1, cmap='viridis')
        elif k == 'kci':
            plt.scatter(rxzs_tiled, ryzs_tiled, c=kci, cmap='viridis')
        if k == 'k1/kci':
            # plt.scatter(rxzs_tiled, ryzs_tiled, c=k1_kci, cmap='viridis', vmax=5., vmin=-1.)
            plt.scatter(rxzs_tiled, ryzs_tiled, c=k1_kci, cmap='viridis')
        # add a vertical and horizontal line at rxy
        # plt.axvline(x=rxy, color='r')
        plt.axhline(y=rxy, color='r')
        plt.colorbar()
        plt.xlabel('rxz')
        plt.ylabel('ryz')
        plt.title(k)
        plt.show()
        plt.close()


plot_k1_and_kci()






