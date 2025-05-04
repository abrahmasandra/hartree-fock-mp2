import json
from pathlib import Path
from pprint import pprint

def load_basis_set_json(basis_name: str, elements: list[str], path: str="data/basis_sets") -> dict:
    """
    Load basis set for specified elements from a JSON file.

    Args:
        basis_name (str): e.g. 'sto-3g'
        elements (list of str): e.g. ['H', 'O']
        path (str): path to basis set files

    Returns:
        dict: element symbol -> list of shells with angular_momentum, exponents, coefficients
    """
    basis_set_path = Path(path) / f"{basis_name.lower()}.json"
    
    if not basis_set_path.exists():
        raise FileNotFoundError(f"Basis set file {basis_set_path} does not exist.")
    
    with open(basis_set_path, 'r') as f:
        basis_data = json.load(f)
    
    # Filter the basis set data for the specified elements
    filtered_basis_data = {el: basis_data[el] for el in elements if el in basis_data}
    
    return filtered_basis_data

if __name__ == "__main__":
    # Example usage
    elements = ["H", "O"]
    basis = load_basis_set_json("sto-3g", elements)

    pprint(basis["H"])