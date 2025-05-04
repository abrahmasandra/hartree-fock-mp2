from pyscf import gto
import numpy as np
from src.molecule import Molecule

def build_pyscf_molecule(molecule: Molecule, basis_name="sto-3g") -> gto.Mole:
    """
    Construct a PySCF Mole object from our Molecule class.

    Args:
        molecule (Molecule): Molecule object containing symbols, coordinates, charge, and multiplicity.
        basis_name (str): Name of the basis set to use.
    Returns:
        gto.Mole: PySCF Mole object.
    """
    atom_str = ""
    for sym, coord in zip(molecule.symbols, molecule.coordinates):
        atom_str += f"{sym} {coord[0]} {coord[1]} {coord[2]}; "

    mol = gto.Mole()
    mol.atom = atom_str
    mol.unit = "Angstrom"
    mol.basis = basis_name
    mol.charge = molecule.charge
    mol.spin = molecule.multiplicity - 1  # multiplicity = 2S + 1
    mol.build()
    return mol