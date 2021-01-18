from .Molecule import Molecule
from .Atom import Atom
import numpy as np
from .GaussianRoutines import exchange_2e_integral


class TwoElectronIntegrals:
    """
    This class contains all the functions related to calculating the
    two electron integrals and their modifications etc.
    """

    def __init__(self):
        pass

    @staticmethod
    def generate_2e_tensor(mol: Molecule) -> np.ndarray:
        """
        this function generates the two electron tensor of
        dimensions n x n x n x n, where n is number of basis functions
        each element of g_{ABCD}
        :param mol:
        :return:
        """
        n = 0
        for atom in mol:
            n += atom.num_gaussian

        Gabcd = np.zeros((n, n, n, n))

        for n1, atom1 in enumerate(mol):
            for i, gauss1 in enumerate(atom1.gaussians):
                for prim1 in gauss1:
                    for n2, atom2 in enumerate(mol):
                        for j, gauss2 in enumerate(atom2.gaussians):
                            for prim2 in gauss2:
                                for n3, atom3 in enumerate(mol):
                                    for k, gauss3 in enumerate(atom3.gaussians):
                                        for prim3 in gauss3:
                                            for n4, atom4 in enumerate(mol):
                                                for l, gauss4 in enumerate(atom4.gaussians):
                                                    for prim4 in gauss4:
                                                        Gabcd[atom1.num_gaussian * n1 + i, atom2.num_gaussian * n2 + j,
                                                              atom3.num_gaussian * n3 + k, atom4.num_gaussian * n4 + l]\
                                                                += exchange_2e_integral(prim1, prim2, prim3, prim4)
        return Gabcd

