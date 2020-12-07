import json as json
from numpy import pi, sqrt
from typing import List


class GTO:
    """simple gaussian function object"""
    def __init__(self, alpha: float, d: float):
        self.alpha = alpha
        self.d = d
        self.nc = (2*self.alpha/pi)**0.75
        self.xp = 0.0
        self.yp = 0.0
        self.zp = 0.0
        self.rp = 0.0

    def gauss_norm(self, n, alpha):
        pass

    def set_coord(self, x: float, y: float, z: float, r=None):
        self.xp = x
        self.yp = y
        self.zp = z
        if r is None:
            self.rp = sqrt(x*x + y*y + z*z)
        else:
            self.rp = r


GaussList = List[GTO]
GaussListList = List[GaussList]


class Basis:
    """This class acts like a simple container for saving basis
    for a single atom. It will have coefficients and other information
    such as number of functions etc. This object will only update what
    it reads from the file. the centers are update in Molecule object file"""
    def __init__(self, basis_file):
        json_obj = ""
        self.basis_file = basis_file
        self.num_gaussian = 0
        self.contractions = []
        self.z = 0
        self.gaussians: GaussListList = []
        self.read_basis()

    def read_basis(self):
        with open(self.basis_file) as fp:
            json_obj = json.load(fp)

        self.z = int(list(json_obj["elements"].keys())[0])

        for shell in json_obj["elements"][str(self.z)]["electron_shells"]:

            self.num_gaussian += 1
            self.contractions.append(len(shell['exponents']))
            contracted_gaussian_list: GaussList = []

            for alpha, d in zip(shell["exponents"], shell["coefficients"][0]):
                contracted_gaussian_list.append(GTO(float(alpha),float(d)))

            self.gaussians.append(contracted_gaussian_list)


if __name__ == '__main__':
    c_basis = Basis("/home/user/programms/pyHF/basis/C.basis")
    for i in range(c_basis.num_gaussian):
        for j in range(c_basis.contractions[i]):
            print(c_basis.gaussians[i][j].alpha, c_basis.gaussians[i][j].d)
