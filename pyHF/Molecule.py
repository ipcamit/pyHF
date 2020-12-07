from typing import List
from .Atom import Atom
# from .Basis import GaussList, GaussListList

AtomList = List[Atom]
StrList = List[str]


class Molecule:
    """
    Molecule class is main class that will store all necessary structural information
     regarding a molecule.
    """
    def __init__(self, xyz_file, basis_file_list: StrList = None):
        self.xyz_file = xyz_file
        self.num_atoms = None
        self.xyz_info = None
        self.basis_file_list = basis_file_list
        self.atoms: AtomList = []
        self.xyz_parser()

    def xyz_parser(self):
        """
        parses xyz file and lists all atom objects with basis file
        :return:
        """
        with open(self.xyz_file) as xyz_file:
            self.num_atoms = int(xyz_file.readline().strip('\n'))
            self.xyz_info = xyz_file.readline()
            for i in range(self.num_atoms):
                curr_line = xyz_file.readline()
                atom_sym, x, y, z = curr_line.split()[0:4]
                atom_obj = Atom(None)
                if self.basis_file_list:
                    for basis_file in self.basis_file_list:
                        if atom_sym + ".basis" in basis_file:
                            atom_obj = Atom(atom_sym,basis_file_name=basis_file)
                else:
                    atom_obj = Atom(atom_sym)
                atom_obj.set_coord(float(x), float(y), float(z))
                self.atoms.append(atom_obj)

    # def set_coords(self, x: float,y: float,z: float):
    #     for atom in self.atoms:
    #         for i in range(atom.num_gaussian):
    #             for j in range(atom.contractions[i]):
    #                 atom.gaussians[i][j].xp = x
    #                 atom.gaussians[i][j].yp = y
    #                 atom.gaussians[i][j].zp = z