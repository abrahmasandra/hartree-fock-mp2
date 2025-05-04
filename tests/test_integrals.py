from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals

def test_integral_shapes():
    mol = Molecule(["H", "H"], [[0,0,0],[0,0,0.74]])
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    
    assert S.shape == T.shape == V.shape
    assert eri.ndim == 4
    assert eri.shape[0] == eri.shape[1] == eri.shape[2] == eri.shape[3]
