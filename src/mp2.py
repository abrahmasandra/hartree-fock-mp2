import numpy as np

## REFERENCE CODE: https://github.com/psi4/psi4numpy/blob/master/Tutorials/05_Moller-Plesset/5a_conventional-mp2.ipynb
def transform_eri_ao_to_mo(eri_ao, C):
    """
    Transform AO-basis ERIs (μν|λσ) to MO-basis ERIs (ij|ab).
    Assumes AO ERIs are in chemist's notation from PySCF ("int2e").

    Args:
        eri_ao (np.ndarray): AO-basis 2-electron integrals (μν|λσ), shape (n,n,n,n)
        C (np.ndarray): MO coefficient matrix, shape (n,n)

    Returns:
        np.ndarray: MO-basis 2-electron integrals (ij|ab), shape (n,n,n,n)
    """
    eri_mo = np.einsum("mnls,mi,nj,la,sb->ijab", eri_ao, C, C, C, C, optimize=True)
    return eri_mo

def compute_mp2_energy(eri_mo, mo_energies, n_occ):
    """
    Compute MP2 correlation energy using MO-basis ERIs and orbital energies.

    Args:
        eri_mo: MO-transformed ERIs, shape (n_orb, n_orb, n_orb, n_orb)
        mo_energies: orbital energies ε_i
        n_occ: number of occupied orbitals

    Returns:
        E_MP2: MP2 correlation energy
    """
    n_orb = eri_mo.shape[0]

    E_MP2 = 0.0
    for i in range(n_occ):
        for j in range(n_occ):
            for a in range(n_occ, n_orb):
                for b in range(n_occ, n_orb):
                    # Compute the numerator and denominator for the MP2 energy
                    iajb = eri_mo[i, a, j, b]
                    ibja = eri_mo[i, b, j, a]
                    numer = iajb * (2 * iajb - ibja)
                    denom = mo_energies[i] + mo_energies[j] - mo_energies[a] - mo_energies[b]
                    E_MP2 += numer / denom
    
    return E_MP2
