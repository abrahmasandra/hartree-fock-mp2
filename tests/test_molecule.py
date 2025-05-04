from src.molecule import Molecule

def test_molecule_electron_count():
    mol = Molecule(["H", "O", "H"], [[0,0,0],[0,0,1],[0,1,0]])
    assert mol.n_electrons == 10
