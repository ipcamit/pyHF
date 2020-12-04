from .Basis import Basis


class Atom(Basis):
    """Contains all information regarding an atom. this call will inherit Basis class.
    Will hold atomic coefficients or other properties in future such as partial charges,
    density matrix coefficient of concerning atoms etc. Not certain where it fits
     right now so keeping it as thin as possible"""
    def __init__(self, symbol, basis_file_name=None):
        self.symbol = symbol
        if basis_file_name:
            self.basis_file = basis_file_name
        else:
            self.basis_file = symbol+".basis"
        super().__init__(self.basis_file)

        self.x = 0
        self.y = 0
        self.z = 0