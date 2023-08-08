
import numpy as np

# Initial conditions

# bdep-BK 2D
Qs0_sq = 0.496  # GeV^2
B_G = 3.2258  # GeV^-2
Nc = 3.
nf = 5.
mZ = 91.1876  # GeV
ml = 0.1  # GeV
mc = 1.3  # GeV
mb = 4.5  # GeV
C = np.sqrt(0.3)
alpha_freeze = 0.761911
alpha_mz = 0.1189

# Internal parameters for running coupling
beta5 = 11. - nf * 2. / 3.
beta4 = 11. - (nf - 1.) * 2. / 3.
beta3 = 11. - (nf - 2.) * 2. / 3.

lambdaQCD5 = mZ * np.exp(-2. * np.pi / (alpha_mz * beta5))
lambdaQCD4 = np.power(mb, 1. - (beta5 / beta4)) * np.power(lambdaQCD5, beta5 / beta4)
lambdaQCD3 = np.power(mc, 1. - (beta4 / beta3)) * np.power(lambdaQCD4, beta4 / beta3)

rsat_sq = 4. * C ** 2 / (lambdaQCD3 ** 2 * np.exp(4. * np.pi / (beta3 * alpha_freeze)))
rc_sq = 4. * C ** 2 / (mc * 2)
rb_sq = 4. * C ** 2 / (mb * 2)






