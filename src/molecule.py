import numpy as np
import periodictable as pt

class Molecule:
    def __init__(self, symbols, coordinates, charge=0, multiplicity=1):
        """
        symbols: list of atomic symbols, e.g. ['H', 'H']
        coordinates: list of 3D coordinates in Angstroms
        charge: net molecular charge
        multiplicity: spin multiplicity (1 for singlet, 2 for doublet, etc.)
        """
        self.symbols = symbols
        self.coordinates = np.array(coordinates)
        self.charge = charge
        self.multiplicity = multiplicity
        self.n_atoms = len(symbols)
        self.n_electrons = self._get_num_electrons()

    def _get_num_electrons(self):
        total = sum(pt.elements.symbol(sym).number for sym in self.symbols)
        return total - self.charge

    def __repr__(self):
        lines = [f"{sym} {xyz[0]: .4f} {xyz[1]: .4f} {xyz[2]: .4f}"
                 for sym, xyz in zip(self.symbols, self.coordinates)]
        return "\n".join(lines)
    
if __name__ == "__main__":
    symbols = ['H', 'H']
    coords = [[0.0, 0.0, 0.0],
              [0.0, 0.0, 0.74]]  # Angstroms

    mol = Molecule(symbols, coords, charge=0)
    print(f"# electrons: {mol.n_electrons}")
    print(mol)        