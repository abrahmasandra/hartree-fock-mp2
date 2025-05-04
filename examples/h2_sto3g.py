from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.scf import run_scf
from src.mp2 import transform_eri_ao_to_mo, compute_mp2_energy
from src.plot_utils import plot_scf_convergence, plot_orbital_energies, plot_density_matrix
from src.utils import build_pyscf_molecule
from src.orbital_plot import evaluate_mo_on_grid, plot_molecular_orbital

if __name__ == "__main__":
    symbols = ["H", "H"]
    coords = [[0, 0, 0], [0, 0, 0.74]]  # H₂ bond length ~0.74 Å
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol)
    eri_ao = compute_2e_integrals(mol)

    # Run SCF
    E_scf, eps, C, D, energy_history = run_scf(S, T, V, eri_ao, mol.n_electrons)
    print(f"HF Energy: {E_scf:.10f} Hartree")
    print("Orbital Energies (Hartree):", eps)
    
    # MP2
    n_occ = mol.n_electrons // 2
    eri_mo = transform_eri_ao_to_mo(eri_ao, C)
    E_mp2_corr = compute_mp2_energy(eri_mo, eps, n_occ)
    E_mp2_total = E_scf + E_mp2_corr

    print(f"MP2 Correlation Energy: {E_mp2_corr:.10f} Hartree")
    print(f"Total MP2 Energy: {E_mp2_total:.10f} Hartree")

    # Plotting MOs
    pyscf_mol = build_pyscf_molecule(mol)
    mo_grid, grid_pts = evaluate_mo_on_grid(pyscf_mol, C, mo_index=0)
    plot_molecular_orbital(
        mo_grid, grid_pts,
        atom_coords=mol.coordinates,
        atom_symbols=mol.symbols,
        iso_value=0.02,
        title="H₂ HOMO (Orbital 1)"
    )
