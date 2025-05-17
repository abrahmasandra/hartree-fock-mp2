import numpy as np
from scipy.linalg import eigh

def build_core_hamiltonian(T, V):
    """
    Build the core Hamiltonian matrix.
    Args:
        T: Kinetic energy matrix
        V: Nuclear attraction matrix
    Returns:
        H_core: Core Hamiltonian matrix
    """
    return T + V

def compute_g_matrix(D, eri):
    """
    Build the two-electron contribution to the Fock matrix.
    G[p,q] = sum_{r,s} D[r,s] * [ (pq|rs) - 0.5 * (pr|qs) ]

    Args:
        D: Density matrix
        eri: Two-electron integrals (pq|rs)
    Returns:
        G: Two-electron contribution to the Fock matrix
    """
    J = np.einsum('rs,pqrs->pq', D, eri)        # Coulomb term
    K = np.einsum('rs,prqs->pq', D, eri)        # Exchange term
    G = 2 * J - K
    return G

def orthogonalizer(S):
    """
    Compute the symmetric orthogonalization matrix A = S^{-1/2}.
    Args:
        S: Overlap matrix
    Returns:
        S_inv_sqrt: Symmetric orthogonalizer
    """
    eigvals, eigvecs = eigh(S)
    S_inv_sqrt = eigvecs @ np.diag(1 / np.sqrt(eigvals)) @ eigvecs.T
    return S_inv_sqrt

def build_fock_matrix(H_core, D, eri):
    """
    Build the Fock matrix.
    Args:
        H_core: Core Hamiltonian matrix
        D: Density matrix
        eri: Two-electron integrals
    Returns:
        F: Fock matrix
    """
    G = compute_g_matrix(D, eri)
    F = H_core + G
    return F

def solve_roothaan_equations(F, S_inv_sqrt):
    """
    Solve the Roothaan equations in the orthonormal basis.
    Args:
        F: Fock matrix
        S_inv_sqrt: Symmetric orthogonalizer
    Returns:
        eps: Orbital energies
        C: MO coefficients in AO basis
    """
    F_prime = S_inv_sqrt @ F @ S_inv_sqrt
    eps, C_prime = eigh(F_prime)
    C = S_inv_sqrt @ C_prime
    return eps, C

def compute_density_matrix(C, n_occ):
    """
    Compute the density matrix from occupied orbitals.
    Args:
        C: MO coefficients
        n_occ: Number of occupied orbitals
    Returns:
        D: Density matrix
    """
    C_occ = C[:, :n_occ]
    D = C_occ @ C_occ.T
    return D

def run_scf(S, T, V, eri, n_electrons, max_iter=50, convergence=1e-8):
    """
    Perform Hartree-Fock SCF calculation.

    Args:
        S: Overlap matrix
        T: Kinetic energy matrix
        V: Nuclear attraction matrix
        eri: Two-electron integrals (pq|rs)
        n_electrons: Total number of electrons
        max_iter: Maximum SCF iterations
        convergence: Energy convergence threshold

    Returns:
        E_total: Final SCF energy
        eps: Orbital energies
        C: MO coefficients
        D: Final density matrix
        energy_history: List of SCF energies per iteration
    """
    n_orb = S.shape[0]
    n_occ = n_electrons // 2  # Restricted HF, closed-shell

    H_core = build_core_hamiltonian(T, V)
    S_inv_sqrt = orthogonalizer(S)

    D = np.zeros((n_orb, n_orb))  # Initial density matrix
    E_total = 0.0
    energy_history = []
    for iteration in range(1, max_iter + 1):
        F = build_fock_matrix(H_core, D, eri)
        eps, C = solve_roothaan_equations(F, S_inv_sqrt)
        D_new = compute_density_matrix(C, n_occ)

        # Compute electronic energy
        E_elec = np.sum((D_new * (H_core + F)))
        delta_E = E_elec - E_total
        E_total = E_elec
        energy_history.append(E_total)

        print(f"SCF Iter {iteration:2d}  E = {E_total:.10f}  dE = {delta_E:.3e}")

        if abs(delta_E) < convergence:
            break

        D = D_new
    else:
        print("SCF did not converge!")

    return E_total, eps, C, D, energy_history

if __name__ == "__main__":
    from molecule import Molecule
    from integrals import compute_1e_integrals, compute_2e_integrals

    symbols = ["H", "H"]
    coords = [[0.0, 0.0, 0.0],
            [0.0, 0.0, 0.74]]
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol, basis_name="sto-3g")
    eri = compute_2e_integrals(mol, basis_name="sto-3g")

    E_scf, eps, C, D, energy_history = run_scf(S, T, V, eri, mol.n_electrons)
    print(f"\nFinal SCF Energy: {E_scf:.8f} Hartree")
