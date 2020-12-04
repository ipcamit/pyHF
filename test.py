from pyHF import Atom
from pyHF import Basis

atom = Atom("C",filename="./basis/C.basis")
print(atom.gaussians)

basis = Basis("./basis/C.basis")
print(basis.gaussians)