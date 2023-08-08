from integrate import integrate


def runge_kutta(calculation, x, y):
    # This function is used to calculate the integral in the BK equation.
    if calculation['order_of_rk'] == 1:
        return runge_kutta_1(calculation, x, y)
    elif calculation['order_of_rk'] == 4:
        return runge_kutta_4(calculation, x, y)
    else:
        print('Runge-Kutta order not implemented.')
        exit(1)


def runge_kutta_1(calculation, x, y):
    y_ind = calculation['y_ind']
    integral_value = integrate(calculation, x, y)
    return integral_value * (calculation['grid']['grid_in_Y'][y_ind] - calculation['grid']['grid_in_Y'][y_ind - 1])


def runge_kutta_4(calculation, x, y):
    y_ind = calculation['y_ind']
    total, kernel, split = integrate(calculation, x, y)
    h = calculation['grid']['grid_in_Y'][y_ind] - calculation['grid']['grid_in_Y'][y_ind - 1]

    k1 = total
    k2 = k1 + h/2.*k1*kernel - h/2.*k1*split - h**2/4. * k1**2 * kernel
    k3 = k2 + h/2.*k2*kernel - h/2.*k2*split - h**2/4. * k2**2 * kernel
    k4 = k2 + h/2.*k3*kernel - h/2.*k3*split - h**2/4. * k3**2 * kernel

    return (h/6.) * (k1 + 2.*k2 + 2.*k3 + k4)

