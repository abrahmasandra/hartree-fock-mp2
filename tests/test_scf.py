from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.scf import run_scf

def test_scf_energy_convergence():
    mol = Molecule(["H", "H"], [[0,0,0],[0,0,0.74]])
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    E_scf, _, _, _, history = run_scf(S, T, V, eri, mol.n_electrons)

    assert abs(history[-1] - E_scf) < 1e-6
    assert len(history) <= 50
    assert E_scf < 0
