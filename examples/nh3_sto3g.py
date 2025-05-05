from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.utils import build_pyscf_molecule
from src.scf import run_scf
from src.orbital_plot import evaluate_mo_on_grid, plot_molecular_orbital

if __name__ == "__main__":
    # Ammonia (NH3), idealized geometry
    symbols = ["N", "H", "H", "H"]
    coords = [
        [ 0.0000,  0.0000,  0.1173],
        [ 0.9377,  0.0000, -0.4692],
        [-0.4688,  0.8126, -0.4692],
        [-0.4688, -0.8126, -0.4692]
    ]
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    E_scf, eps, C, D, _ = run_scf(S, T, V, eri, mol.n_electrons)

    pyscf_mol = build_pyscf_molecule(mol)
    mo_grid, grid_pts = evaluate_mo_on_grid(pyscf_mol, C, mo_index=4)  # HOMO

    plot_molecular_orbital(
        mo_grid, grid_pts,
        atom_coords=mol.coordinates,
        atom_symbols=mol.symbols,
        title="NH3 HOMO"
    )
