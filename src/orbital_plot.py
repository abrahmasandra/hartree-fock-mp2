import numpy as np
from pyscf import gto
import pyvista as pv

def plot_molecular_orbital(
    orbital_grid,
    grid_points,
    atom_coords,
    atom_symbols,
    iso_value=0.02,
    title="Molecular Orbital"
):
    """
    Visualize molecular orbital isosurface with PyVista.

    Args:
        orbital_grid: 3D numpy array of MO values on a grid
        grid_points: (X, Y, Z) meshgrid arrays
        atom_coords: (N, 3) array of nuclear positions
        atom_symbols: list of atomic symbols
        iso_value: isosurface threshold
        title: plot window title
    """
    x, y, z = grid_points
    spacing = (x[1,0,0] - x[0,0,0], y[0,1,0] - y[0,0,0], z[0,0,1] - z[0,0,0])
    origin = (x.min(), y.min(), z.min())
    dims = orbital_grid.shape

    # Create PyVista grid
    grid = pv.ImageData()
    grid.dimensions = dims
    grid.origin = origin
    grid.spacing = spacing
    grid.point_data["MO"] = orbital_grid.flatten(order="F")

    # Start plotter
    plotter = pv.Plotter()
    plotter.add_title(title)

    # Add positive and negative isosurfaces
    pos_surf = grid.contour([iso_value], scalars="MO")
    neg_surf = grid.contour([-iso_value], scalars="MO")

    if pos_surf.n_cells > 0:
        plotter.add_mesh(pos_surf, color="blue", opacity=0.5, label="+")

    if neg_surf.n_cells > 0:
        plotter.add_mesh(neg_surf, color="red", opacity=0.5, label="-")

    # Add atoms as spheres
    for sym, pos in zip(atom_symbols, atom_coords):
        sphere = pv.Sphere(radius=0.3, center=pos)
        plotter.add_mesh(sphere, color="black")
        plotter.add_point_labels([pos], [sym], point_size=0, font_size=14)

    plotter.show()

def evaluate_mo_on_grid(mol: gto.Mole, C, mo_index=0, grid_spacing=0.1, box_size=5.0):
    """
    Evaluate a molecular orbital on a cubic grid.

    Returns:
        grid_values: 3D array
        (X, Y, Z): meshgrid arrays for plotting
    """
    # Make 3D grid
    x = y = z = np.arange(-box_size, box_size + grid_spacing, grid_spacing)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    coords = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

    # Evaluate AO basis functions
    ao_values = mol.eval_gto("GTOval_cart", coords)  # shape: (n_grid, n_basis)

    # Evaluate MO
    mo_coeff = C[:, mo_index]
    mo_values = ao_values @ mo_coeff  # (n_grid,)
    mo_grid = mo_values.reshape(X.shape)

    return mo_grid, (X, Y, Z)