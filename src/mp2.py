import numpy as np

def transform_eri_ao_to_mo(eri_ao, C):
    """
    Transform AO-basis ERIs (μν|λσ) to MO-basis ERIs (pq|rs).
    Assumes AO ERIs are in chemist's notation from PySCF ("int2e").

    Args:
        eri_ao (np.ndarray): AO-basis 2-electron integrals (μν|λσ), shape (n,n,n,n)
        C (np.ndarray): MO coefficient matrix, shape (n,n)

    Returns:
        np.ndarray: MO-basis 2-electron integrals (pq|rs), shape (n,n,n,n)
    """
    return np.einsum("munv,mi,nj,uk,vl->ijkl", eri_ao, C, C, C, C, optimize=True)

def compute_mp2_energy(mo_eri, mo_energies, n_occ):
    """
    Compute MP2 correlation energy using MO-basis ERIs and orbital energies.

    Args:
        mo_eri: MO-transformed ERIs, shape (n_orb, n_orb, n_orb, n_orb)
        mo_energies: orbital energies ε_i
        n_occ: number of occupied orbitals

    Returns:
        E_MP2: MP2 correlation energy
    """
    n_orb = mo_eri.shape[0]
    E_MP2 = 0.0

    for i in range(n_occ):
        for j in range(n_occ):
            for a in range(n_occ, n_orb):
                for b in range(n_occ, n_orb):
                    ijab = mo_eri[i, j, a, b]
                    ijba = mo_eri[i, j, b, a]
                    denom = mo_energies[i] + mo_energies[j] - mo_energies[a] - mo_energies[b]
                    E_MP2 += (ijab - ijba) ** 2 / denom

    return E_MP2