import numpy as np
import const as c


def mareks_N(grid, dimensionality_of_N):
    grid_in_r = grid['grid_in_r']
    if dimensionality_of_N == 2:
        grid_in_b = grid['grid_in_b']
        N = np.zeros((len(grid_in_r), len(grid_in_b)))
        for r_ind in range(len(grid_in_r)):
            for b_ind in range(len(grid_in_b)):
                bq_sq = grid_in_r[r_ind]**2/4. + grid_in_b[b_ind]**2
                T = np.exp(-bq_sq**2/(2.*c.B_G)) + np.exp(-bq_sq**2/(2.*c.B_G))
                N[r_ind][b_ind] = 1. - np.exp(-0.5 * c.Qs0_sq/4. * grid_in_r[r_ind]**2 * T)
    else:
        print('Initial condition not implemented')
        exit(1)
    return N


def MV_1D(grid, Qs0_sq, Lambda, gamma):
    grid_in_r = grid['grid_in_r']
    N = np.zeros(len(grid_in_r))
    for r_ind in range(len(grid_in_r)):
        N[r_ind] = 1. - np.exp(- (grid_in_r[r_ind]**2 * Qs0_sq)**gamma/4. * np.log(1./(grid_in_r[r_ind] * Lambda) + np.e))
    return N


def mareks_N_3D(grid):
    grid_in_r = grid['grid_in_r']
    grid_in_b = grid['grid_in_b']
    grid_in_theta = grid['grid_in_theta']

    N = np.zeros((len(grid_in_r), len(grid_in_b), len(grid_in_theta)))

    for r_ind in range(len(grid_in_r)):
        for b_ind in range(len(grid_in_b)):
            for theta_ind in range(len(grid_in_theta)):
                r = grid_in_r[r_ind]
                b = grid_in_b[b_ind]
                theta = grid_in_theta[theta_ind]
                x_quark = np.array([-r / 2. * np.cos(theta), b - r / 2. * np.sin(theta)])
                y_quark = np.array([r / 2. * np.cos(theta), b + r / 2. * np.sin(theta)])
                bqx_sq = np.linalg.norm(x_quark)**2
                bqy_sq = np.linalg.norm(y_quark)**2

                T = np.exp(-bqx_sq**2/(2.*c.B_G)) + np.exp(-bqy_sq**2/(2.*c.B_G))
                N[r_ind][b_ind][theta_ind] = 1. - np.exp(-0.5 * c.Qs0_sq/4. * grid_in_r[r_ind]**2 * T)
    return N


def mareks_N_4D(grid):
    grid_in_r = grid['grid_in_r']
    grid_in_b = grid['grid_in_b']
    grid_in_theta = grid['grid_in_theta']
    grid_in_phi = grid['grid_in_phi']

    N = np.zeros((len(grid_in_r), len(grid_in_b), len(grid_in_theta), len(grid_in_phi)))

    for r_ind in range(len(grid_in_r)):
        for b_ind in range(len(grid_in_b)):
            for theta_ind in range(len(grid_in_theta)):
                for phi_ind in range(len(grid_in_phi)):
                    r = grid_in_r[r_ind]
                    b = grid_in_b[b_ind]
                    theta = grid_in_theta[theta_ind]
                    phi = grid_in_phi[phi_ind]
                    x_quark = np.array([-r / 2. * np.cos(theta) + b * np.cos(phi), b - r / 2. * np.sin(theta) + b * np.sin(phi)])
                    y_quark = np.array([r / 2. * np.cos(theta) + b * np.cos(phi), b + r / 2. * np.sin(theta) + b * np.sin(phi)])

                    bqx_sq = np.linalg.norm(x_quark) ** 2
                    bqy_sq = np.linalg.norm(y_quark) ** 2

                    # Make a dipole-like deformation in the initial condition
                    if x_quark[1] < 0.:
                        bqx_sq *= 0.3
                    if y_quark[1] < 0.:
                        bqy_sq *= 0.3

                    T = np.exp(-bqx_sq ** 2 / (2. * c.B_G)) + np.exp(-bqy_sq ** 2 / (2. * c.B_G))
                    N[r_ind][b_ind][theta_ind][phi_ind] = 1. - np.exp(-0.5 * c.Qs0_sq / 4. * grid_in_r[r_ind] ** 2 * T)
    return N


