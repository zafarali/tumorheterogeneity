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


def dChi_fn(alpha, beta, x_c, x_c_threshold=None):
    """
    Returns a function that calculates the derivative of
    the chi function.
    """
    x_c = float(x_c)
    if x_c_threshold is None: x_c_threshold = x_c

    def dChi_(x):
        if x / x_c_threshold < 1:
            return (-alpha / x_c) * (x / x_c) ** (-alpha - 1)

        elif x / x_c_threshold > 1:
            return (-beta / x_c) * (x / x_c) ** (-beta - 1)
        else:
            return 0

    return dChi_


def x_c_prime_fn(alpha, beta, x_c):
    return x_c * (alpha/beta)**(1/(alpha-beta))

def SFS(N, mu, alpha, beta, with_correction=True):
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
    x_c_prime = x_c_prime_fn(alpha, beta, x_c) if with_correction else None
    dChi = dChi_fn(alpha, beta, x_c, x_c_threshold=x_c_prime)

    def fn_(x):
        conditional_factor = dChi(x)
        return -phi_c * N * mu * conditional_factor

    return fn_, x_c_prime

