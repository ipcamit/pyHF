from .Molecule import Molecule


class NuclearNuclearRepulsion:
    @staticmethod
    def nuclear_nuclear_energy(mol:Molecule):
        En = 0
        for i in range(mol.num_atoms):
            for j in range(i+1, mol.num_atoms):
                En += mol.atoms[i].Z * mol.atoms[j].Z / (
                    (mol.atoms[i].xc - mol.atoms[j].xc)**2 +
                    (mol.atoms[i].yc - mol.atoms[j].yc)**2 +
                    (mol.atoms[i].zc - mol.atoms[j].zc)**2
                )
        return En
