import numpy as np


def function_to_integrate(x, y):
    return np.exp(-x**2 - y**2)


def polar_MC(polar):
    size = 100000
    integral = 0.
    integration_radius = 4.
    if polar:
        for _ in range(size):
            r = np.random.random()*integration_radius
            phi = np.random.random()*2.*np.pi
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            integral += function_to_integrate(x, y) * r
        integral = integral * 2.*np.pi * integration_radius / size
    else:
        for _ in range(size):
            length = 2. * integration_radius
            x = np.random.random()*length - length/2.
            y = np.random.random()*length - length/2.
            integral += function_to_integrate(x, y)
        integral = integral * length**2 / size
    print('POLAR: True integral should be pi ', '; MC:', integral, polar)


def log_MC(log):
    size = 100000
    integral = 0.
    if log:
        for _ in range(size):
            x = np.random.uniform(-2, 7.)
            jacobian_MC_log = (10**x * np.log(10))*9.
            integral += 10**x * jacobian_MC_log
            # x = np.random.uniform(np.log(10**-2), np.log(10**7.))
            # integral += np.e**x
            # (np.log(10**7) - np.log(10**-2))
        integral = integral / size

    else:
        for _ in range(size):
            x = np.random.uniform(10**-2, 10**7)
            integral += x
        integral = integral*10**7 / size

    print('LOG: True integral should be 0.5*10**7*10**7 = 5*10**13; MC:', integral/10**13, '* 10**13', log)


def both_combined_the_Mareks_way():
    no_of_samples = 500000
    steps_in_integrand_theta = 20
    integral = 0.

    shift = 2. * np.pi / steps_in_integrand_theta / 2.  # to avoid double counting and y-axis with z and w.
    integrand_angles = np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta)
    integrand_radii = np.logspace(-7., 2., 20000)
    integrand_cartesian_coods = np.zeros((len(integrand_angles) * len(integrand_radii), 2))
    for radius_ind in range(len(integrand_radii)):
        for theta_ind in range(len(integrand_angles)):
            x_cood = integrand_radii[radius_ind] * np.cos(integrand_angles[theta_ind])
            y_cood = integrand_radii[radius_ind] * np.sin(integrand_angles[theta_ind])
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][0] = x_cood
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][1] = y_cood

    random_indexes_z = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    array_of_z = integrand_cartesian_coods[random_indexes_z]
    for integration_point in array_of_z:
        x = integration_point[0]
        y = integration_point[1]
        jacobian = np.linalg.norm(integration_point) * np.log(10.) * np.linalg.norm(integration_point)  # One r for polar integration and ln(10)*r for log
        # integral += function_to_integrate(x, y)*jacobian
        integral += 1.*jacobian
    probability_normalization_polar = 2. * np.pi
    probability_normalization_log = np.log10(integrand_radii[-1]) - np.log10(integrand_radii[0])
    print(probability_normalization_log, probability_normalization_polar)
    integral = integral * probability_normalization_polar * probability_normalization_log / no_of_samples
    print('True integral should be pi; MC:', integral)


def both_combined_the_Mareks_way_2D():
    no_of_samples = 1000000
    steps_in_integrand_theta = 20
    integral = 0.

    shift = 2. * np.pi / steps_in_integrand_theta / 2.  # to avoid double counting and y-axis with z and w.
    integrand_angles = np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta)
    integrand_radii = np.logspace(-7., 2., 20000)
    integrand_cartesian_coods = np.zeros((len(integrand_angles) * len(integrand_radii), 2))
    for radius_ind in range(len(integrand_radii)):
        for theta_ind in range(len(integrand_angles)):
            x_cood = integrand_radii[radius_ind] * np.cos(integrand_angles[theta_ind])
            y_cood = integrand_radii[radius_ind] * np.sin(integrand_angles[theta_ind])
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][0] = x_cood
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][1] = y_cood

    random_indexes_z = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    random_indexes_w = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    array_of_z = integrand_cartesian_coods[random_indexes_z]
    array_of_w = integrand_cartesian_coods[random_indexes_w]
    # THE FUNCTION IS A UNIT FUNCTION
    function_to_integrate = 1.
    for integration_idx in range(no_of_samples):
        integration_point_z = array_of_z[integration_idx]
        integration_point_w = array_of_w[integration_idx]
        jacobian_z = np.linalg.norm(integration_point_z) * np.log(10.) * np.linalg.norm(integration_point_z)  # One r for polar integration and ln(10)*r for log
        jacobian_w = np.linalg.norm(integration_point_w) * np.log(10.) * np.linalg.norm(integration_point_w)  # One r for polar integration and ln(10)*r for log
        integral += function_to_integrate*jacobian_z*jacobian_w
    probability_normalization_polar = 2. * np.pi
    probability_normalization_log = np.log10(integrand_radii[-1]) - np.log10(integrand_radii[0])
    probability_normalization = probability_normalization_polar * probability_normalization_log

    integral = integral * probability_normalization**2 / no_of_samples
    print('True integral should be pi*10000**2 = 98696; MC:', integral/986960000.)


def function_to_integrate_importance(importance, integration_point_z, integration_point_w):
    if importance:
        # / np.exp(-np.power(radii_importance - (-1.0), 2.) / (2 * np.power(2.0, 2.)))
        return 1. / np.exp(-np.power(np.linalg.norm(integration_point_z) - 1.0, 2.) / (2 * np.power(2.0, 2.))) / np.exp(-np.power(np.linalg.norm(integration_point_w) - 1.0, 2.) / (2 * np.power(2.0, 2.)))
    else:
        return 1. * np.zeros(len(integration_point_z))

def both_combined_the_Mareks_way_2D_normal_distribution_sampling():
    importance = True
    no_of_samples = 1000000
    steps_in_integrand_theta = 20
    integral = 0.

    shift = 2. * np.pi / steps_in_integrand_theta / 2.  # to avoid double counting and y-axis with z and w.
    integrand_angles = np.linspace(-np.pi + shift, np.pi - shift, steps_in_integrand_theta)

    if importance:
        radii_importance = np.random.normal(loc=-1.0, scale=2.0, size=20000)
        integrand_radii = np.power(10, radii_importance)
    else:
        integrand_radii = np.logspace(-7., 2., 20000)


    integrand_cartesian_coods = np.zeros((len(integrand_angles) * len(integrand_radii), 2))
    for radius_ind in range(len(integrand_radii)):
        for theta_ind in range(len(integrand_angles)):
            x_cood = integrand_radii[radius_ind] * np.cos(integrand_angles[theta_ind])
            y_cood = integrand_radii[radius_ind] * np.sin(integrand_angles[theta_ind])
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][0] = x_cood
            integrand_cartesian_coods[radius_ind * len(integrand_angles) + theta_ind][1] = y_cood

    random_indexes_z = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    random_indexes_w = np.random.randint(0, len(integrand_cartesian_coods), no_of_samples)  # z
    array_of_z = integrand_cartesian_coods[random_indexes_z]
    array_of_w = integrand_cartesian_coods[random_indexes_w]

    for integration_idx in range(no_of_samples):
        integration_point_z = array_of_z[integration_idx]
        integration_point_w = array_of_w[integration_idx]
        jacobian_z = np.linalg.norm(integration_point_z) * np.log(10.) * np.linalg.norm(integration_point_z)  # One r for polar integration and ln(10)*r for log
        jacobian_w = np.linalg.norm(integration_point_w) * np.log(10.) * np.linalg.norm(integration_point_w)  # One r for polar integration and ln(10)*r for log
        functional_value = function_to_integrate_importance(importance, integration_point_z, integration_point_w)
        integral += functional_value*jacobian_z*jacobian_w
    probability_normalization_polar = 2. * np.pi
    probability_normalization_log = np.log10(integrand_radii[-1]) - np.log10(integrand_radii[0])
    probability_normalization = probability_normalization_polar * probability_normalization_log

    integral = integral * probability_normalization**2 / no_of_samples
    print('True integral should be pi*10000**2 = 98696; Analytic/MC=', integral/986960000.)



# polar_MC(polar=True)
# polar_MC(polar=False)
#
# log_MC(log=True)
# log_MC(log=False)

both_combined_the_Mareks_way()
# both_combined_the_Mareks_way_2D()
# both_combined_the_Mareks_way_2D_normal_distribution_sampling()




