import numpy as np

def transform_eri_ao_to_mo(eri_ao, C):
    """
    Transform AO-basis ERIs (μν|λσ) to MO-basis ERIs (pq|rs).
    Full 4-index transformation: (μν|λσ) -> (pq|rs)

    Args:
        eri_ao (np.ndarray): Two-electron integrals in AO basis (μν|λσ).
                             Shape: (n_orb, n_orb, n_orb, n_orb)
        C (np.ndarray): MO coefficients matrix.
                        Shape: (n_orb, n_mo)
    Returns:
        np.ndarray: Two-electron integrals in MO basis (pq|rs).
                    Shape: (n_mo, n_mo, n_mo, n_mo)
    """
    return np.einsum("pqrs,pi,qj,rk,sl->ijkl", eri_ao, C, C, C, C, optimize=True)

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