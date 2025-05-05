from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.utils import build_pyscf_molecule
from src.scf import run_scf
from src.orbital_plot import evaluate_mo_on_grid, plot_molecular_orbital

if __name__ == "__main__":
    # Carbon monoxide (CO), bond length ~1.13 Å
    symbols = ["C", "O"]
    coords = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.13]]
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    E_scf, eps, C, D, _ = run_scf(S, T, V, eri, mol.n_electrons)

    pyscf_mol = build_pyscf_molecule(mol)
    mo_grid, grid_pts = evaluate_mo_on_grid(pyscf_mol, C, mo_index=6)  # HOMO

    plot_molecular_orbital(
        mo_grid, grid_pts,
        atom_coords=mol.coordinates,
        atom_symbols=mol.symbols,
        title="CO HOMO"
    )
