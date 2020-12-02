import json as json


class GTO:
    """simple gaussian function object"""
    def __init__(self, alpha: float, d: float):
        self.alpha = alpha
        self.d = d


class Basis:
    """This class acts like a simple container for saving basis
    for a single atom. It will have coefficients and other information
    such as number of functions etc."""
    def __init__(self, basis_file):
        json_obj = ""
        
        with open(basis_file) as fp:
            json_obj = json.load(fp)
        
        self.z = int(list(json_obj["elements"].keys())[0])
        self.num_gaussian = 0
        self.contractions = []

        self.gaussians = []

        for shell in json_obj["elements"][str(self.z)]["electron_shells"]:

            self.num_gaussian += 1
            self.contractions.append(len(shell['exponents']))
            contracted_gaussian_list = []

            for alpha,d in zip(shell["exponents"], shell["coefficients"][0]):
                contracted_gaussian_list.append(GTO(float(alpha),float(d)))

            self.gaussians.append(contracted_gaussian_list)


if __name__ == '__main__':
    c_basis = Basis("/home/user/programms/pyHF/basis/C.basis")
    for i in range(c_basis.num_gaussian):
        for j in range(c_basis.contractions[i]):
            print(c_basis.gaussians[i][j].alpha, c_basis.gaussians[i][j].d)
