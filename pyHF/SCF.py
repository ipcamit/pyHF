from .Molecule import Molecule
import numpy as np
from pyHF import OneElectronIntegrals as oe
from pyHF import TwoElectronIntegrals as te
from .NuclearNuclearIntegrals import NuclearNuclearRepulsion as nnr

import numpy.linalg as la

class SCF:
    """This will be one base class for all main loop calculations.
    other methods will inherit it and modify it. Starting simply with HF"""
    def __init__(self, molecule: Molecule):
        """
        This routine initializes an SCF cycle on a given molecule. The kind of calculations would
        be carried out by "kind" parameter. As of now only HF is supported, MP2 next? To be inherited by
        respective method?
        :param molecule:
        :param kind:
        """
        self.H0: np.ndarray = np.zeros((1, 1))
        self.T: np.ndarray = np.zeros((1, 1))
        self.Vc: np.ndarray = np.zeros((1, 1))
        self.Gabcd: np.ndarray = np.zeros((1, 1))
        self.S = np.zeros((1, 1))
        self.P: np.ndarray = np.zeros((1,1))
        self.C: np.ndarray = np.zeros((1,1))
        self.G: np.ndarray = np.zeros((1,1))
        self.molecule: Molecule = molecule
        self.initialize()

    def initialize(self):
        """
        start building all the integrals.
        :return:
        """
        self.S = oe.OneElectronIntegrals.overlap(self.molecule)
        self.T = oe.OneElectronIntegrals.kinetic_energy(self.molecule)
        self.Vc = oe.OneElectronIntegrals.coulombic_potential_energy(self.molecule)
        self.Gabcd = te.TwoElectronIntegrals.generate_2e_tensor(self.molecule)
        self.H0 = self.T + self.Vc

    def get_P(self):
        self.P = np.zeros(self.C.shape)
        for mu in range(self.C.shape[0]):
            for nu in range(self.C.shape[1]):
                for a in range(int(self.molecule.Ne/2)):
                    self.P[mu, nu] += 2*self.C[mu, a]*self.C[nu,a]

    def get_G(self):
        self.G = np.zeros(self.P.shape)
        for i in range(self.P.shape[0]):
            for j in range(self.P.shape[0]):
                for x in range(self.P.shape[0]):
                    for y in range(self.P.shape[0]):
                        self.G[i, j] += self.P[x, y]*(self.Gabcd[i, j, y, x] \
                                                        - 0.5*self.Gabcd[i, x, y, j])

    def eval_energy(self):
        E = 0.0
        for i in range(self.P.shape[0]):
            for j in range(self.P.shape[1]):
                E += self.P[i,j] * (self.H0[i, j] + self.F[i, j]) # 2H + G
        return E


class HF(SCF):
    """
    HF SCF cycle
    """
    def __init__(self,molecule: Molecule, maxcycle=20, tol=1e-4):
        self.molecule = molecule
        self.H0: np.ndarray = None
        self.Vc: np.ndarray = None
        self.Gabcd: np.ndarray = None
        self.T: np.ndarray = None
        self.S: np.ndarray = None
        self.F: np.ndarray = None
        self.C:np.ndarray = None
        self.maxcycle = maxcycle
        self.tol = tol
        super(HF, self).__init__(self.molecule)

    def energy(self):
        """get energy"""
        # STEP 3 Diagonalize overlap matrix
        s: np.ndarray = np.zeros((1, self.S.shape[0]))
        U: np.ndarray = np.zeros(self.S.shape)
        s, U = la.eig(self.S)
        # sort eigs
        # Get transformation matrix (Szabo 3.167)
        X = U.dot(np.sqrt(np.diag(1/s)).dot(U.T))
        # STEP 4 initial guess from H0 and get P
        self.F = self.H0
        print(self.H0)
        print(self.S)
        Fp = X.T.dot(self.F.dot(X))
        Cp, Cp_vec = la.eig(Fp)
        self.C = X.dot(Cp_vec)
        self.get_P()
        E0 = self.eval_energy()
        iterations = 0
        E = 0.
        self.get_G()
        print(self.G)
        while True:
            # STEP 5 calculate G from C
            self.get_G()
            # print(self.G)
            # Calculate E
            self.F = self.H0 + self.G
            Fp = X.T.dot(self.F.dot(X))
            Cp, Cp_vec = la.eig(Fp)
            self.C = X.dot(Cp_vec)
            self.get_P()
            print(E0)
            E0 = E
            E = self.eval_energy()
            if np.abs(E-E0) < self.tol:
                print(E,E0)
                break
            iterations +=1
            if iterations > self.maxcycle:
                break

        print("Final Energy after {} iterations: {}".format(iterations, E))
        En = nnr.nuclear_nuclear_energy(self.molecule)
        print("Final Nuclear:  {} Total: {}".format(En, E + En))
