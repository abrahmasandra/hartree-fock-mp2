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
pip install -r requirements.txt
```

## Usage

Need to write this section

## Testing

Run unit tests using pytest:

```bash
pytest tests/
```

Each module includes at least one unit test verifying correctness of outputs, such as integral symmetries, energy convergence, and MP2 correction values.
