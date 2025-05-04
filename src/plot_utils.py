import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_orbital_energies(orbital_energies, n_occ, title="Orbital Energies", show=True):
    """
    Plot occupied and virtual orbital energies.
    """
    x = np.arange(len(orbital_energies))
    occ = orbital_energies[:n_occ]
    virt = orbital_energies[n_occ:]

    plt.figure(figsize=(8, 5))
    plt.scatter(x[:n_occ], occ, color='blue', label='Occupied')
    plt.scatter(x[n_occ:], virt, color='red', label='Virtual')
    plt.axhline(0, color='gray', lw=0.5)
    plt.xlabel("Orbital Index")
    plt.ylabel("Energy (Hartree)")
    plt.title(title)
    plt.legend()
    if show:
        plt.show()

def plot_scf_convergence(energies, title="SCF Convergence", show=True):
    """
    Plot SCF total energy vs iteration.
    """
    iterations = np.arange(1, len(energies) + 1)
    plt.figure(figsize=(8, 5))
    plt.plot(iterations, energies, marker='o')
    plt.xlabel("SCF Iteration")
    plt.ylabel("Total Energy (Hartree)")
    plt.title(title)
    plt.grid(True)
    if show:
        plt.show()

def plot_density_matrix(D, title="Density Matrix Heatmap", show=True):
    """
    Heatmap of the density matrix.
    """
    plt.figure(figsize=(6, 5))
    sns.heatmap(D, annot=True, cmap='coolwarm', fmt=".3f")
    plt.title(title)
    if show:
        plt.show()
