import pyHF as ph
from pyHF import OneElectronIntegrals as oe
import numpy as np

mol = ph.Molecule('xyz/h2.xyz',basis_file_list=['basis/H.basis'])

T = oe.OneElectronIntegrals.coulombic_potential_energy(mol)
print(T)
