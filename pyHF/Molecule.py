from .Atom import Atom


class Molecule:
    def __init__(self, xyz_file, basis_file_list=None):
        self.xyz_file = xyz_file
        self.num_atoms = None
        self.xyz_info = None
        self.basis_file_list = basis_file_list
        self.atoms = []
        self.xyz_parser()

    def xyz_parser(self):
        with open(self.xyz_file) as xyz_file:
            self.num_atoms = int(xyz_file.readline().strip('\n'))
            self.xyz_info = xyz_file.readline()
            for i in range(self.num_atoms):
                curr_line = xyz_file.readline()
                atom_sym, x, y, z = curr_line.split()[0:3]
                if self.basis_file_list:
                    for basis_file in self.basis_file_list:
                        if atom_sym + ".basis" in basis_file:
                            atom_obj = Atom(atom_sym,basis_file_name=basis_file)
                else:
                    atom_obj = Atom(atom_sym)
                atom_obj.x = float(x)
                atom_obj.y = float(y)
                atom_obj.z = float(z)
                self.atoms.append(atom_obj)