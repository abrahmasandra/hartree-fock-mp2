from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.scf import run_scf
from src.orbital_plot import evaluate_mo_on_grid, plot_molecular_orbital
from src.utils import build_pyscf_molecule

if __name__ == "__main__":
    # Water molecule (H2O) in Angstroms
    symbols = ["O", "H", "H"]
    coords = [
        [ 0.000000,  0.000000,  0.000000],
        [ 0.758602,  0.000000,  0.504284],
        [-0.758602,  0.000000,  0.504284]
    ]
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    E_scf, eps, C, D, _ = run_scf(S, T, V, eri, mol.n_electrons)

    # Use PySCF for orbital evaluation
    pyscf_mol = build_pyscf_molecule(mol)
    mo_grid, grid_pts = evaluate_mo_on_grid(pyscf_mol, C, mo_index=3)  # HOMO for H2O

    plot_molecular_orbital(
        mo_grid, grid_pts,
        atom_coords=mol.coordinates,
        atom_symbols=mol.symbols,
        title="H2O LUMO"
    )
