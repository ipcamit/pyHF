import numpy as np
from .Molecule import Molecule
from .Atom import Atom
from .Basis import GTO
from scipy.special import erf, gamma, gammainc

"""
This file simply contains all the functions needed for manipulation of gaussian
basis set functions. current list include:
1. gaussian multiplication
2. gaussian integral
3. orthogonalization? (not sure of that yet)
"""


def gauss_product(gauss1: GTO, gauss2: GTO) -> GTO:
    """
    simple gaussian close form integral. takes in two gto objects
    (simple structure with N, alpha, d), and returns a new GTO object.

    :param gauss1:
    :param gauss2:
    :return: gauss_final
    """
    p = (gauss2.alpha + gauss1.alpha)
    q = gauss2.alpha * gauss1.alpha
    K = np.exp(-(q / p) * (gauss1.rp - gauss2.rp) ** 2)
    xp = (gauss2.alpha * gauss2.xp + gauss1.alpha * gauss1.xp) / p
    yp = (gauss2.alpha * gauss2.yp + gauss1.alpha * gauss1.yp) / p
    zp = (gauss2.alpha * gauss2.zp + gauss1.alpha * gauss1.zp) / p
    Rp = np.sqrt(xp * xp + yp * yp + zp * zp)
    d = gauss2.d * gauss1.d
    gauss_final = GTO(alpha=p, d=d * K)
    gauss_final.set_coord(xp, yp, zp, r=Rp)
    gauss_final.nc = gauss1.nc * gauss2.nc
    return gauss_final


def gauss_int(gauss1: GTO) -> float:
    """
    Takes in GTO type object and return its integral -infty to +infty
    :param gauss1
    :return:
    """
    return gauss1.nc * gauss1.d * np.sqrt(np.pi / gauss1.alpha) ** 3


def ke_gauss_int(gauss1: GTO, gauss2: GTO) -> float:
    """
    takes in two gaussian type orbital objects and return kinetic energy
    matrix T. szabo pp 412
    :param gauss1:
    :param gauss2:
    :return:
    """
    pre_factor = (gauss1.alpha * gauss2.alpha) / (gauss1.alpha + gauss2.alpha)
    del_R = (gauss2.xp - gauss1.xp) ** 2 + (gauss2.yp - gauss1.yp) ** 2 + (gauss2.zp - gauss1.zp) ** 2
    integral = pre_factor * (3 - 2 * pre_factor * del_R) * \
               (np.pi / (gauss1.alpha + gauss2.alpha)) ** (3 / 2) * np.exp(-pre_factor * del_R)
    return integral * gauss2.d * gauss2.nc * gauss1.d * gauss1.nc


def coulombic_integral(gauss1: GTO, gauss2: GTO, center_atom: Atom) -> float:
    """
    Takes in two basis function and center atom object to return following integral:
    <gauss1|-Z_{center_atom}/r_{center_atom}|gauss2>
    szabo pp415, A33
    :param gauss1:
    :param gauss2:
    :param center_atom:
    :return:
    """
    sum_alpha = gauss1.alpha + gauss2.alpha
    prod_alpha = gauss1.alpha * gauss2.alpha
    del_R = (gauss2.xp - gauss1.xp) ** 2 + (gauss2.yp - gauss1.yp) ** 2 + (gauss2.zp - gauss1.zp) ** 2
    gauss3 = gauss_product(gauss1, gauss2)
    del_Rc = (gauss3.xp - center_atom.xc) ** 2\
           + (gauss3.yp - center_atom.yc) ** 2\
           + (gauss3.zp - center_atom.zc) ** 2
    integral = -2 * np.pi / sum_alpha * center_atom.Z * np.exp(-prod_alpha / sum_alpha * del_R) \
               * F0(sum_alpha * del_Rc)
    return integral * gauss2.d * gauss2.nc * gauss1.d * gauss1.nc


def F0(t: float) -> float:
    """
    Fourier integral thingy for t^-1/2 \int_0^{t/2}dy exp^-{y^2}
    special case of  Kummer confluent hypergeometric function. Will be modified suitably in future
    right now just erf based implementation as in szabo.
    :param t:
    :return:
    """
    if t == 0:
        return 1.0  # internet says it. bout couldn't find proper reference for this implementation
    else:
        return 0.5 * np.sqrt(np.pi / t) * erf(np.sqrt(t))


def Fo(n: int, t: float) -> float:
    if t == 0:
        output = 1 / (2 * n + 1)
    else:
        output = gamma(n + 0.5) * gammainc(t, n + 0.5) / (2 * t ** (n + 0.5))
    return output
