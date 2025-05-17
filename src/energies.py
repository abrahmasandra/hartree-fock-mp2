import numpy as np
from periodictable import elements
from src.utils import build_pyscf_molecule
from src.molecule import Molecule

def compute_nuclear_repulsion_energy(mol):
    """
    Compute nuclear repulsion energy.

    Args:
        mol (Molecule): Molecule object

    Returns:
        float: nuclear repulsion energy
    """
    pyscf_mol = build_pyscf_molecule(mol)
    return pyscf_mol.energy_nuc()

def compute_total_energy(e_scf, mol: Molecule):
    """
    Compute total energy (SCF + nuclear repulsion).

    Args:
        e_scf (float): SCF electronic energy
        mol (Molecule): Molecule object

    Returns:
        float: total energy
    """
    e_nuc = compute_nuclear_repulsion_energy(mol)
    return e_scf + e_nuc
