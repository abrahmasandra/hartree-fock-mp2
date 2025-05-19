import argparse
from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.utils import build_pyscf_molecule
from src.scf import run_scf
from src.mp2 import transform_eri_ao_to_mo, compute_mp2_energy
from src.energies import compute_total_energy
from src.orbital_plot import evaluate_mo_on_grid, plot_molecular_orbital
from presets import MOLECULE_PRESETS

from pyscf import gto, scf, mp

def get_pyscf_reference(mol: Molecule):
    mol = build_pyscf_molecule(mol)

    hf = scf.RHF(mol).run()
    mp2_calc = mp.MP2(hf).run()
    return hf.e_tot, mp2_calc.e_corr, mp2_calc.e_tot

def main():
    parser = argparse.ArgumentParser(description="Run HF + MP2 and visualize MO.")
    parser.add_argument('--visualize', action='store_true', default=False, help="Visualize orbitals?")
    parser.add_argument("--mo-index", type=int, default=0, help="Molecular orbital index to visualize")
    parser.add_argument("--output", type=str, default="results.txt", help="Output file for energies")
    parser.add_argument(
        "--molecule", type=str, default="h2o", choices=MOLECULE_PRESETS.keys(),
        help="Name of the molecule to run (default: h2o)"
    )
    args = parser.parse_args()

    # Define molecule (default: H2O)
    mol = MOLECULE_PRESETS[args.molecule]
    pyscf_mol = build_pyscf_molecule(mol)
    print(f"Loaded molecule: {args.molecule.upper()}")
    print(f"Number of electrons: {mol.n_electrons}")

    # Integrals
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    print(f"Number of basis functions (orbitals): {S.shape[0]}")
    
    # SCF
    print()
    print(" > Running SCF procedure...")
    E_scf, eps, C, D, history = run_scf(S, T, V, eri, mol.n_electrons)
    E_scf_total = compute_total_energy(E_scf, mol)
    print(" > Done.")

    # MP2
    print()
    print(" > Computing MP2 correction...")
    eri_mo = transform_eri_ao_to_mo(eri, C)
    E_mp2_corr = compute_mp2_energy(eri_mo, eps, mol.n_electrons // 2)
    E_mp2_total = E_scf_total + E_mp2_corr
    print(" > Done.")

    # PySCF
    print()
    print(" > Computing PySCF reference values...")
    E_scf_total_ref, E_mp2_corr_ref, E_mp2_total_ref = get_pyscf_reference(mol)
    print(" > Done.")

    # Report energies
    print()
    print("--- Our energies ---")
    print(f"SCF Elec. Energy: {E_scf:.8f} Eh")
    print(f"SCF Total Energy: {E_scf_total:.8f} Eh")
    print(f"MP2 Corr. Energy: {E_mp2_corr:.8f} Eh")
    print(f"MP2 Total Energy: {E_mp2_total:.8f} Eh")

    # Reference energies from PySCF
    print()
    print("--- PySCF energies ---")
    print(f"SCF Total Energy: {E_scf_total_ref:.8f} Eh")
    print(f"MP2 Corr. Energy: {E_mp2_corr_ref:.8f} Eh")
    print(f"MP2 Total Energy: {E_mp2_total_ref:.8f} Eh")

    # Errors
    print()
    print("--- Our error ---")
    print(f"SCF Total Energy: {abs(E_scf_total - E_scf_total_ref):.8f} Eh")
    print(f"MP2 Corr. Energy: {abs(E_mp2_corr - E_mp2_corr_ref):.8f} Eh")
    print(f"MP2 Total Energy: {abs(E_mp2_total - E_mp2_total_ref):.8f} Eh")

    # Write to file
    with open(args.output, "w") as f:
        f.write("--- Our energies ---")
        f.write(f"SCF Elec. Energy: {E_scf:.8f} Eh")
        f.write(f"SCF Total Energy: {E_scf_total:.8f} Eh")
        f.write(f"MP2 Corr. Energy: {E_mp2_corr:.8f} Eh")
        f.write(f"MP2 Total Energy: {E_mp2_total:.8f} Eh")
    
        # Reference energies from PySCF
        f.write("\n")
        f.write("--- PySCF energies ---")
        f.write(f"SCF Total Energy: {E_scf_total_ref:.8f} Eh")
        f.write(f"MP2 Corr. Energy: {E_mp2_corr_ref:.8f} Eh")
        f.write(f"MP2 Total Energy: {E_mp2_total_ref:.8f} Eh")
    
        # Errors
        f.write("\n")
        f.write("--- Our error ---")
        f.write(f"SCF Total Energy: {abs(E_scf_total - E_scf_total_ref):.8f} Eh")
        f.write(f"MP2 Corr. Energy: {abs(E_mp2_corr - E_mp2_corr_ref):.8f} Eh")
        f.write(f"MP2 Total Energy: {abs(E_mp2_total - E_mp2_total_ref):.8f} Eh")

    # Visualize MO
    if args.visualize:
        pyscf_mol = build_pyscf_molecule(mol)
        mo_grid, grid_pts = evaluate_mo_on_grid(pyscf_mol, C, mo_index=args.mo_index)
    
        plot_molecular_orbital(
            mo_grid, grid_pts,
            atom_coords=mol.coordinates,
            atom_symbols=mol.symbols,
            title=f"Molecular Orbital {args.mo_index}"
        )

if __name__ == "__main__":
    main()