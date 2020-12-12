from .Basis import Basis
from numpy import sqrt


class Atom(Basis):
    """Contains all information regarding an atom. this call will inherit Basis class.
    Will hold atomic coefficients or other properties in future such as partial charges,
    density matrix coefficient of concerning atoms etc. Not certain where it fits
     right now so keeping it as thin as possible"""
    def __init__(self, symbol, basis_file_name=None):
        self.xc = 0
        self.yc = 0
        self.zc = 0
        if symbol is None:
            pass
        else:
            self.symbol = symbol
            if basis_file_name:
                self.basis_file = basis_file_name
            else:
                self.basis_file = symbol+".basis"
            super().__init__(self.basis_file)

    def set_coord(self, x: float, y: float, z: float):
        for i in range(self.num_gaussian):
            for j in range(self.contractions[i]):
                self.gaussians[i][j].xp = x
                self.gaussians[i][j].yp = y
                self.gaussians[i][j].zp = z
                self.xc = x
                self.yc = y
                self.zc = z
                self.gaussians[i][j].rp = sqrt(x*x + y*y + z*z)