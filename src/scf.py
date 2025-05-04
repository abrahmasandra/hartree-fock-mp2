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
    """
    n_orb = S.shape[0]
    n_occ = n_electrons // 2  # Restricted HF, closed-shell

    H_core = build_core_hamiltonian(T, V)

    # Orthogonalizer A = S^{-1/2}
    eigvals, eigvecs = eigh(S)
    S_inv_sqrt = eigvecs @ np.diag(1 / np.sqrt(eigvals)) @ eigvecs.T

    D = np.zeros((n_orb, n_orb))  # Initial density matrix
    E_total = 0.0
    energy_history = []
    for iteration in range(1, max_iter + 1):
        G = compute_g_matrix(D, eri)
        F = H_core + G

        # Step 1: Transform Fock matrix to orthonormal basis
        F_prime = S_inv_sqrt @ F @ S_inv_sqrt

        # Step 2: Solve standard eigenvalue problem (this is the secular equation in orthonormal basis)
        eps, C_prime = eigh(F_prime)

        # Step 3: Transform back to original (non-orthogonal) AO basis
        C = S_inv_sqrt @ C_prime
        C_occ = C[:, :n_occ]

        D_new = C_occ @ C_occ.T

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
