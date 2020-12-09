from pyscf import gto, scf

mol = gto.M(atom='H 0 0 0; H 0 0 1.0', basis='6-31g')
mf = scf.RHF(mol)

# Overlap matrix for H2
print(mol.intor('int1e_ovlp_cart'))