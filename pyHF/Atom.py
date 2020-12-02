from .Basis import Basis


class Atom(Basis):
    """Contains all information regarding an atom. this call will inherit Basis class.
    Will hold atomic coefficients or other properties in future such as partial charges,
    density matrix coefficient of concerning atoms etc. Not certain where it fits
     right now so keeping it as thin as possible"""
    def __init__(self, symbol, filename=None):
        self.symbol = symbol
        if filename:
            self.basis_file = filename
        else:
            self.basis_file = symbol+".basis"
        super().__init__(self.basis_file)
