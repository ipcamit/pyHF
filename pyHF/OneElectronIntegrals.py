from .Molecule import Molecule
import numpy as np
from .GaussianRoutines import gauss_int, gauss_product, ke_gauss_int, coulombic_integral


class OneElectronIntegrals:
    """As one electron integrals remain unchanged almost always,
    this class will contain all the required routines. Majority of
    gaussian manipulation subroutines are written in gaussian routines file"""

    def __init__(self):
        pass

    @staticmethod
    def kinetic_energy(mol: Molecule) -> np.ndarray:
        """
        Calculates core hamiltonian, or Kinetic energy Operator T.

        :param mol:
        :return T:
        """
        n = 0
        for atom in mol:
            n += atom.num_gaussian
        T = np.zeros((n, n))
        for n1, atom1 in enumerate(mol):
            for n2, atom2 in enumerate(mol):
                for i, gauss1 in enumerate(atom1.gaussians):
                    for j, gauss2 in enumerate(atom2.gaussians):
                        for prim1 in gauss1:
                            for prim2 in gauss2:
                                T[atom1.num_gaussian * n1 + i, atom2.num_gaussian * n2 + j] \
                                    += ke_gauss_int(prim1, prim2)
        return T

    @staticmethod
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
                                S[atom1.num_gaussian * n1 + i, atom2.num_gaussian * n2 + j] \
                                    += gauss_int(gauss_product(prim1, prim2))
        return S

    @staticmethod
    def coulombic_potential_energy(mol: Molecule) -> np.ndarray:
        """
        Calculates core hamiltonian, or Kinetic energy Operator T.

        :param mol:
        :return T:
        """
        n = 0
        for atom in mol:
            n += atom.num_gaussian
        Vc = np.zeros((n, n))
        for center_atom in mol:
            for n1, atom1 in enumerate(mol):
                for n2, atom2 in enumerate(mol):
                    for i, gauss1 in enumerate(atom1.gaussians):
                        for j, gauss2 in enumerate(atom2.gaussians):
                            for prim1 in gauss1:
                                for prim2 in gauss2:
                                    Vc[atom1.num_gaussian * n1 + i, atom2.num_gaussian * n2 + j] \
                                                += coulombic_integral(prim1, prim2, center_atom)
                                    # print(Vc)
            # exit()
        return Vc
