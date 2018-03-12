import numpy as np


def x_c_fn(N, alpha, beta):
    """
    Calculates the critical x_c
    """
    return N ** (-(1 - alpha) / (beta - alpha))


def phi_c_fn(N, alpha, beta):
    """
    Calculates the critical phi
    """
    return N ** (-(alpha * beta - alpha) / (beta - alpha))


def dChi_fn(alpha, beta, x_c):
    """
    Returns a function that calculates the derivative of
    the chi function.
    """
    x_c = float(x_c)

    def dChi_(x):
        if x / x_c < 1:
            return (-alpha / x_c) * (x / x_c) ** (-alpha - 1)

        elif x / x_c > 1:
            return (-beta / x_c) * (x / x_c) ** (-beta - 1)
        else:
            return 0

    return dChi_


def SFS(N, mu, alpha, beta):
    """
    Calculates the SFS
    :param N: size of the tumor
    :param mu: the mutation rate
    :param alpha: scaling factor
    :param beta: scaling factor
    """
    N, alpha, beta = float(N), float(alpha), float(beta)
    x_c = x_c_fn(N, alpha, beta)
    phi_c = phi_c_fn(N, alpha, beta)
    dChi = dChi_fn(alpha, beta, x_c)
    print(np.log10(phi_c))
    print(np.log10(x_c))

    def fn_(x):
        conditional_factor = dChi(x)
        return -phi_c * N * mu * conditional_factor

    return fn_

