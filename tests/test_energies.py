from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals
from src.scf import run_scf
from src.mp2 import transform_eri_ao_to_mo, compute_mp2_energy

from pyscf import gto, scf, mp

def compare_energies(symbols, coords, hf_tol=1e-3, mp2_tol=1e-3):
    mol = Molecule(symbols, coords)
    S, T, V = compute_1e_integrals(mol)
    eri = compute_2e_integrals(mol)
    E_scf, eps, C, D, _ = run_scf(S, T, V, eri, mol.n_electrons)
    
    eri_mo = transform_eri_ao_to_mo(eri, C)
    E_mp2_corr = compute_mp2_energy(eri_mo, eps, mol.n_electrons // 2)
    E_mp2 = E_scf + E_mp2_corr

    # Run PySCF reference calculation
    pyscf_mol = gto.Mole()
    pyscf_mol.atom = [[sym, coord] for sym, coord in zip(symbols, coords)]
    pyscf_mol.basis = "sto-3g"
    pyscf_mol.unit = "Angstrom"
    pyscf_mol.charge = mol.charge
    pyscf_mol.spin = mol.multiplicity - 1
    pyscf_mol.build()

    hf_calc = scf.RHF(pyscf_mol).run()
    mp2_calc = mp.MP2(hf_calc).run()

    assert abs(E_scf - hf_calc.e_tot) < hf_tol, f"HF energy mismatch: {E_scf} vs {hf_calc.e_tot}"
    assert abs(E_mp2 - mp2_calc.e_tot) < mp2_tol, f"MP2 energy mismatch: {E_mp2} vs {mp2_calc.e_tot}"

def test_h2_comparison():
    symbols = ["H", "H"]
    coords = [[0, 0, 0], [0, 0, 0.74]]
    compare_energies(symbols, coords)

def test_h2o_comparison():
    symbols = ["O", "H", "H"]
    coords = [
        [0.000000,  0.000000,  0.000000],
        [0.758602,  0.000000,  0.504284],
        [-0.758602, 0.000000,  0.504284]
    ]
    compare_energies(symbols, coords)
