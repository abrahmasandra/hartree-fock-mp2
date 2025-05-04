import numpy as np
from pyscf import gto
from src.utils import build_pyscf_molecule

def compute_1e_integrals(molecule, basis_name="sto-3g"):
    """
    Returns overlap (S), kinetic (T), and nuclear attraction (V) matrices using PySCF.
    """
    mol = build_pyscf_molecule(molecule, basis_name)

    S = mol.intor("int1e_ovlp")     # Overlap
    T = mol.intor("int1e_kin")      # Kinetic
    V = mol.intor("int1e_nuc")      # Nuclear attraction

    return np.array(S), np.array(T), np.array(V)

def compute_2e_integrals(molecule, basis_name="sto-3g"):
    """
    Compute 4-index two-electron integrals (ERI) in AO basis using PySCF.

    Returns:
        eri[p, q, r, s] = (pq|rs)
    """
    mol = build_pyscf_molecule(molecule, basis_name)
    eri = mol.intor("int2e")  # shape: (nbf, nbf, nbf, nbf)
    return np.array(eri)


if __name__ == "__main__":
    from molecule import Molecule

    symbols = ["H", "H"]
    coords = [[0.0, 0.0, 0.0],
              [0.0, 0.0, 0.74]]
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol, basis_name="sto-3g")

    print("Overlap matrix S:\n", S)
    print("Kinetic matrix T:\n", T)
    print("Nuclear attraction matrix V:\n", V)
