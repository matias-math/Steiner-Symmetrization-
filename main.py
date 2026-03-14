

import numpy as np
import pyvista as pv
from scipy.interpolate import RegularGridInterpolator

def run_protocol():
    print("--- Steiner Symmetrization - 3D Visual Representation ---")
    f_str = input("Define Set A (Implicit Equation): ").strip()
    v_str = input("Vector a (Direction): ")
    
    # 1. Vector Normalization (Any vector 'a' is accepted and normalized)
    vec_raw = np.array([float(i) for i in v_str.split(',')], dtype=float)
    a = vec_raw / np.linalg.norm(vec_raw)
    
    # 2. Domain Expansion (The '3000' scale / Expanded Box)
    limit = 10.0 
    res = 400 # Visual resolution for a smooth manifold
    
    # 3. Setup Orthonormal Basis for Plane Pa
    # Pa := { x in Rn | x · a = 0 }
    if np.abs(a[0]) < 0.8:
        ref = np.array([1, 0, 0])
    else:
        ref = np.array([0, 1, 0])
    u = np.cross(a, ref); u /= np.linalg.norm(u)
    v = np.cross(a, u)

    # 4. Discrete Analysis of H1 Measure
    # We sweep the plane Pa at high density to find the measure of intersections
    map_res = 300 
    plane_coords = np.linspace(-limit, limit, map_res)
    t_samples = np.linspace(-limit*2, limit*2, 1000) 
    dt = t_samples[1] - t_samples[0]
    
    h1_map = np.zeros((map_res, map_res))
    
    print(f"Analyzing fibers Lb^a in a fixed domain of size [-{limit}, {limit}]...")
    for i in range(map_res):
        for j in range(map_res):
            # Point b in Pa
            b = plane_coords[i] * u + plane_coords[j] * v
            # Line Lb^a := { b + ta | t in R }
            line_pts = b + np.outer(t_samples, a)
            
            # Evaluate the indicator function of the set A
            inside = eval(f_str, {"__builtins__": None}, 
                         {"x": line_pts[:,0], "y": line_pts[:,1], "z": line_pts[:,2], "np": np}) <= 0
            
            if np.any(inside):
                # H1 is the 1D Hausdorff measure of the intersection
                h1_map[i, j] = np.sum(inside) * dt

    # 5. Construction of the Symmetrized Set Sa(A)
    lin = np.linspace(-limit, limit, res)
    x, y, z = np.meshgrid(lin, lin, lin, indexing='ij')
    grid_pts = np.stack([x, y, z], axis=-1).reshape(-1, 3)
    
    t_grid = np.dot(grid_pts, a)
    b_grid = grid_pts - np.outer(t_grid, a)
    
    # Interpolation of H1 measures over the 3D domain
    interp = RegularGridInterpolator((plane_coords, plane_coords), h1_map, 
                                     bounds_error=False, fill_value=0)
    
    voxel_h1 = interp(np.stack([np.dot(b_grid, u), np.dot(b_grid, v)], axis=-1))
    
    # Sa(A) Condition: |t| <= 1/2 * H1(A ∩ Lb^a)
    mask_sa = (np.abs(t_grid) <= (voxel_h1 / 2.0)).reshape(res, res, res)

    # 6. Formal Visualization
    plotter = pv.Plotter(shape=(1, 2), window_size=[1600, 800])
    spacing = (2*limit/(res-1), 2*limit/(res-1), 2*limit/(res-1))
    origin = (-limit, -limit, -limit)

    # Symmetry Plane Pa Reference
    plane_mesh = pv.Plane(center=(0, 0, 0), direction=a, i_size=limit*2, j_size=limit*2)

    for i, (m, title, color) in enumerate([
        (eval(f_str, {"__builtins__": None}, {"x":x, "y":y, "z":z, "np":np}) <= 0, "Original Set A", "silver"), 
        (mask_sa, "Steiner Symmetrization Sa(A)", "magenta")]):
        
        plotter.subplot(0, i)
        grid = pv.ImageData(dimensions=(res, res, res), spacing=spacing, origin=origin)
        if np.any(m):
            mesh = grid.contour(isosurfaces=[0.5], scalars=m.flatten(order="F").astype(float))
            plotter.add_mesh(mesh, color=color, opacity=0.8, smooth_shading=True)
        
        plotter.add_mesh(plane_mesh, color="cyan", opacity=0.1, show_edges=True, label="Plane Pa")
        plotter.add_arrows(np.array([[0,0,0]]), a.reshape(1,3), mag=limit/2, color='yellow')
        plotter.add_text(title, font_size=12)
        plotter.show_grid() # Shows the numerical scale of the larger domain

    plotter.link_views()
    plotter.show()

if __name__ == "__main__":
    run_protocol()