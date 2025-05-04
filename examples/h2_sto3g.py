from src.molecule import Molecule
from src.basis import load_basis_set_json
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.scf import run_scf
from src.mp2 import transform_eri_ao_to_mo, compute_mp2_energy

if __name__ == "__main__":
    symbols = ["H", "H"]
    coords = [[0, 0, 0], [0, 0, 0.74]]  # H₂ bond length ~0.74 Å
    mol = Molecule(symbols, coords)

    basis = load_basis_set_json("sto-3g", symbols)
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)

    # Run SCF
    E_scf, eps, C, D = run_scf(S, T, V, eri, mol.n_electrons)
    print(f"HF Energy: {E_scf:.10f} Hartree")
    print("Orbital Energies (Hartree):", eps)
    
    # MP2
    n_occ = mol.n_electrons // 2
    mo_eri = transform_eri_ao_to_mo(eri, C)
    E_mp2_corr = compute_mp2_energy(mo_eri, eps, n_occ)
    E_mp2_total = E_scf + E_mp2_corr

    print(f"MP2 Correlation Energy: {E_mp2_corr:.10f} Hartree")
    print(f"Total MP2 Energy: {E_mp2_total:.10f} Hartree")
