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

```plaintext
hartree-fock-mp2/
├── src/ # Core implementation code
├── tests/ # Unit tests
├── examples/ # Example HF + MP2 runs
├── README.md
├── requirements.txt
```

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

## Running the Main Program

You can run the main driver script using:

```bash
python main.py --molecule h2o --mo-index 3 --output results.txt
```

This will run the Hartree-Fock calculation for the water molecule (H₂O) and output the results to `results.txt`. You can specify different molecules and molecular orbital indices as needed.

### 🔧 Available arguments

- `--molecule`  
  Select one of the predefined molecules to simulate. Options include:
  - `h2`
  - `h2o`
  - `nh3`
  - `ch4`  
  Default: `h2o`

- `--mo-index`  
  The index of the molecular orbital to visualize. Use:
  - HOMO ≈ `n_electrons // 2 - 1`
  - LUMO ≈ `n_electrons // 2`

- `--output`  
  Name of the file to save computed SCF and MP2 energy results.  
  Default: `results.txt`

---

### 📥 Example: Run MP2 calculation on ammonia and visualize the HOMO

```bash
python main.py --molecule nh3 --mo-index 4 --output nh3_energy.txt
```

This will:

- Run Hartree-Fock and MP2 energy calculations  
- Print and save the results to `nh3_energy.txt`  
- Display the 3D isosurface plot of molecular orbital #4

## Usage (Python Code)

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
   E_scf_total = compute_total_energy(E_scf, mol) # adds nuclear-nuclear repulsion term
   ```

4. **Compute MP2 Energy**
    Use the `transform_eri_ao_to_mo` to preprocess eri to the MO basis, and then use `compute_mp2_energy(eri_mo, eps, n_occ)` function to calculate the MP2 correlation energy.

    ```python
    n_occ = mol.n_electrons // 2
    eri_mo = transform_eri_ao_to_mo(eri_ao, C)
    E_mp2_corr = compute_mp2_energy(eri_mo, eps, n_occ)
    E_total = E_scf_total + E_mp2_corr
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

## Examples

The `examples/` folder contains ready-to-run scripts demonstrating how to use this package on real molecules. These examples include full Hartree-Fock workflows, optional MP2 correlation energy corrections, and 3D molecular orbital visualization using PyVista.

### 📄 Available examples

- `h2_sto3g.py` – Minimal H₂ molecule with HF and MP2
- `h2o_sto3g.py` – Water molecule, visualize HOMO/LUMO
- `nh3_sto3g.py` – Ammonia molecule with orbital plotting
- `co_sto3g.py` – Carbon monoxide, small polar molecule

To run an example:

```bash
python -m examples.h2o_sto3g
```

💡 You can modify any script in examples/ to experiment with different geometries, basis sets, or orbitals. Visualization defaults to PyVista but can be adapted.

## Testing

Run unit tests using pytest:

```bash
python -m pytest tests/
```

Each module includes at least one unit test verifying correctness of outputs, such as integral symmetries, energy convergence, and MP2 correction values.
