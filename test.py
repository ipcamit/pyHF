import pyHF as ph
from pyHF import OneElectronIntegrals as oe
from pyHF import TwoElectronIntegrals as te
import numpy as np

mol = ph.Molecule('xyz/h2.xyz',basis_file_list=['basis/H.basis'])

T = te.TwoElectronIntegrals.generate_2e_tensor(mol)
print(T[:,:,0,0])
