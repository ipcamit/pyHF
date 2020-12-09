from .Molecule import Molecule
import numpy as np
class OneElectronIntegrals:
    """As one electron integrals remain unchanged almost always,
    this class will contain all the required routines. Majority of
    gaussian manipulation subroutines are written in gaussian routines file"""

    def __init__(self):
        pass

    def hamiltonian_core(self,mol: Molecule) -> np.ndarray:
        """
        Calculates core hamiltonian, or Kinetic energy Operator T.

        :param mol:
        :return T:
        """
        n = 0
        for atom in mol:
            n += atom.num_gaussian
        T = np.zeros((n, n))