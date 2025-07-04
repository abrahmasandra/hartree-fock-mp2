from src.molecule import Molecule

MOLECULE_PRESETS = {
    "h2": Molecule(
        symbols=["H", "H"],
        coordinates=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]
    ),
    "h2o": Molecule(
        symbols=["O", "H", "H"],
        coordinates=[
            [0.000000,  0.000000,  0.000000],
            [0.758602,  0.000000,  0.504284],
            [-0.758602, 0.000000,  0.504284]
        ]
    ),
    "nh3": Molecule(
        symbols=["N", "H", "H", "H"],
        coordinates=[
            [0.0000,  0.0000,  0.1173],
            [0.9377,  0.0000, -0.4692],
            [-0.4688,  0.8126, -0.4692],
            [-0.4688, -0.8126, -0.4692]
        ]
    ),
    "ch4": Molecule(
        symbols=["C", "H", "H", "H", "H"],
        coordinates=[
            [0.0000, 0.0000, 0.0000],
            [0.6291, 0.6291, 0.6291],
            [-0.6291, -0.6291, 0.6291],
            [-0.6291, 0.6291, -0.6291],
            [0.6291, -0.6291, -0.6291]
        ]
    ),
    "hf": Molecule(
        symbols=["H", "F"],
        coordinates=[
            [0.0000, 0.0000, 0.0000],
            [0.0000, 0.0000, 0.9176]
        ]
    ),
}
