from src.molecule import Molecule
from src.integrals import compute_1e_integrals, compute_2e_integrals

if __name__ == "__main__":
    """
    Example usage of the compute_1e_integrals and compute_2e_integrals functions.
    """

    # Define a simple H2 molecule
    symbols = ["H", "H"]
    coords = [[0.0, 0.0, 0.0],
            [0.0, 0.0, 0.74]]
    mol = Molecule(symbols, coords)

    S, T, V = compute_1e_integrals(mol, basis_name="sto-3g")
    eri = compute_2e_integrals(mol, basis_name="sto-3g")

    print("One-electron integrals:")
    print("S:\n", S)
    print("T:\n", T)
    print("V:\n", V)

    print("\nTwo-electron integrals shape:", eri.shape)
    print("Example eri[0,0,0,0] = (00|00):", eri[0, 0, 0, 0])
