# Hartree-Fock + MP2 Quantum Chemistry Code

This repository contains a Python implementation of the Hartree-Fock (HF) method and Møller–Plesset perturbation theory to second order (MP2) for performing ab initio quantum chemistry calculations.

The goal of this project is to build a simple, modular, and educational codebase for computing electronic energies of small molecules using minimal basis sets such as STO-3G. This project is suitable for students learning computational chemistry and quantum mechanics.

---

## Features

✅ Hartree-Fock Self-Consistent Field (SCF) procedure  
✅ MP2 post-Hartree-Fock correlation energy correction  
✅ Minimal basis set support (STO-3G, others optional)  
✅ Modular code structure with unit tests

---

## File Structure

hartree-fock-mp2/
├── src/ # Core implementation code
├── tests/ # Unit tests
├── examples/ # Example HF + MP2 runs
├── data/ # Basis set data
├── README.md
├── requirements.txt

## Installation

Ensure Python 3.8+ is installed. Clone the repo and install dependencies:

```bash
git clone https://github.com/abrahmasandra/hartree-fock-mp2.git
cd hartree-fock-mp2
```

Install dependencies:

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
pip install -r requirements.txt
```

## Usage

This package implements Hartree-Fock (HF) and MP2 quantum chemistry calculations from first principles, with support for molecular orbital visualization.

To use it, follow this general workflow:

1. **Define a molecule**  
    Use the `Molecule` class to specify atomic symbols and Cartesian coordinates in Angstroms.

    ```python
    symbols = ["H", "H"]
    coords = [[0, 0, 0], [0, 0, 0.74]]  # H₂ bond length ~0.74 Å
    mol = Molecule(symbols, coordinates)
    ```

2. **Compute one and two electron integrals**
    Use the `compute_1e_integrals` and `compute_2e_integrals` functions to get the overlap (S), kinetic (T), potential (V), and electron repulsion integrals (ERI) in the AO basis.

    ```python
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    ```

3. **Run the SCF Procedure**  
   Use the `run_scf(S, T, V, eri, n_electrons)` function to perform the SCF procedure and obtain the SCF total energy, orbital energies, coefficients, and density matrix.

   ```python
   E_scf, eps, C, D, energy_history = run_scf(S, T, V, eri, mol.n_electrons)
   ```

4. **Compute MP2 Energy**
    Use the `transform_eri_ao_to_mo` to preprocess eri to the MO basis, and then use `compute_mp2_energy(eri_mo, eps, n_occ)` function to calculate the MP2 correlation energy.

    ```python
    n_occ = mol.n_electrons // 2
    eri_mo = transform_eri_ao_to_mo(eri_ao, C)
    E_mp2_corr = compute_mp2_energy(eri_mo, eps, n_occ)
    E_total = E_scf + E_mp2_corr
    ```

5. **Visualize Molecular Orbitals (Optional)**
    Use `build_pyscf_mol` to create a PySCF-compatible molecule object for visualization. Then use `evaluate_mo_on_grid` to compute the molecular orbital values on a grid. Finally, use `plot_molecular_orbital` to visualize the molecular orbitals.

    ```python
    pyscf_mol = build_pyscf_molecule(mol)
    mo_grid, grid_pts = evaluate_mo_on_grid(pyscf_mol, C, mo_index=0)
    plot_molecular_orbital(
        mo_grid, grid_pts,
        atom_coords=mol.coordinates,
        atom_symbols=mol.symbols,
        iso_value=0.02,
        title="H₂ HOMO (Orbital 1)"
    )
    ```

    To visualize other orbitals (e.g. LUMO), simply change `mo_index`.
    💡 Tip: For a closed-shell molecule with `n` electrons, the HOMO is at index `n // 2 - 1` and the LUMO is at index `n // 2`.

## Testing

Run unit tests using pytest:

```bash
python -m pytest tests/
```

Each module includes at least one unit test verifying correctness of outputs, such as integral symmetries, energy convergence, and MP2 correction values.
