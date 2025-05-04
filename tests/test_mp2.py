from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.scf import run_scf
from src.mp2 import transform_eri_ao_to_mo, compute_mp2_energy

def test_mp2_energy_sign():
    mol = Molecule(["H", "H"], [[0,0,0],[0,0,0.74]])
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    E_scf, eps, C, D, _ = run_scf(S, T, V, eri, mol.n_electrons)

    eri_mo = transform_eri_ao_to_mo(eri, C)
    E_mp2 = compute_mp2_energy(eri_mo, eps, mol.n_electrons // 2)

    assert E_mp2 <= 0
    assert abs(E_mp2) < 0.1  # Should be small for H2 with STO-3G
