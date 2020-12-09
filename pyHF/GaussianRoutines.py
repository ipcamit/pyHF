import numpy as np
from .Molecule import Molecule
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
    d = gauss2.d * gauss1.d
    gauss_final = GTO(alpha=p, d=d*K)
    gauss_final.set_coord(xp, yp, zp, r=Rp)
    gauss_final.nc = gauss1.nc * gauss2.nc
    print(gauss_final.nc,gauss_final.alpha, gauss_final.d, gauss_final.rp)
    return gauss_final


def gauss_int(gauss1: GTO) -> float:
    """Takes in GTO type object and return its integral -infty to +infty"""
    return gauss1.d*np.sqrt(gauss1.alpha/np.pi)**3


def overlap(mol: Molecule) -> np.ndarray:
    """
    returns the overlap matrix
    :param mol:
    :return S:
    """
    n = 0
    for atom in mol:
        n += atom.num_gaussian
    S = np.zeros((n, n))
    for n1, atom1 in enumerate(mol):
        for n2, atom2 in enumerate(mol):
            for i, gauss1 in enumerate(atom1.gaussians):
                for j, gauss2 in enumerate(atom2.gaussians):
                    for prim1 in gauss1:
                        for prim2 in gauss2:
                            S[atom1.num_gaussian * n1 + i, atom2.num_gaussian * n2 + j]\
                                += gauss_int(gauss_product(prim1, prim2))
                            # print(gauss_int(gauss_product(prim1, prim2)))
    return S
