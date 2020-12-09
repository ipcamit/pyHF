import pyHF as ph
import numpy as np

mol = ph.Molecule('xyz/h2.xyz',basis_file_list=['basis/H.basis'])

S = ph.overlap(mol)
print(S)
