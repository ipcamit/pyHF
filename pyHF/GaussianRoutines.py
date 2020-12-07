import numpy as np
from .Atom import Atom
from .Basis import GTO
"""
This file simply contains all the functions needed for manipulation of gaussian
basis set functions. current list include:
1. gaussian multiplication
2. gaussian integral
3. orthogonalization? (not sure of that yet)
"""


def gauss_product(gauss1 :GTO, gauss2 :GTO) -> GTO:
    """
    simple gaussian close form integral. takes in two gto objects
    (simple structure with N, alpha, d), and returns a new GTO object.

    :param gauss1:
    :param gauss2:
    :return: gauss_final
    """
    p = (gauss2.alpha + gauss1.alpha)
    q = gauss2.alpha*gauss1.alpha
    K = np.exp(-(q/p) * (gauss1.rp - gauss2.rp)**2)
    xp = (gauss2.alpha * gauss2.xp + gauss1.alpha * gauss1.xp)/p
    yp = (gauss2.alpha * gauss2.yp + gauss1.alpha * gauss1.yp)/p
    zp = (gauss2.alpha * gauss2.zp + gauss1.alpha * gauss1.zp)/p
    Rp = (gauss2.alpha * gauss2.rp + gauss1.alpha * gauss1.rp)/p
    d = gauss2.d*gauss1.d
    gauss_final = GTO(alpha=p,d=d*K)
    gauss_final.set_coord(xp, yp, zp, r=Rp)
    return gauss_final


def gauss_int(gauss1: GTO) -> float:
    """Takes in GTO type object and return its integral -infty to +infty"""
    return gauss1.d*np.sqrt(gauss1.alpha/np.pi)**3


def overlap(atom1: Atom, atom2: Atom) -> np.ndarray:
    n1 = 0
    n2 = 0
    g1 = 0
    g2 = 0
    N = atom1.num_gaussian + atom2.num_gaussian
    S = np.zeros((N, N))
    atom_list = [atom1, atom2]
    for n1 in range(N):
        for n2 in range(N):
            for g1 in range(atom_list[n1].contractions[n1]):
                for g2 in range(atom_list[2].contractions[n2]):
                    S[n1, n2] += gauss_int(gauss_product(atom[].gaussians[n1][g1],atom[].gaussians[n2][g2]))
    return S
