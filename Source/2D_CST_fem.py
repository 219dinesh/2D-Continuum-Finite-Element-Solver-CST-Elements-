import numpy as np

# ==========================================
# 1. DEFINE GEOMETRY, MESH, & MATERIALS
# ==========================================
# A 2m x 1m rectangular plate
# Nodes: [X, Y] in meters
nodes = np.array([
    [0.0, 0.0],  # Node 0 (Bottom Left)
    [2.0, 0.0],  # Node 1 (Bottom Right)
    [2.0, 1.0],  # Node 2 (Top Right)
    [0.0, 1.0]   # Node 3 (Top Left)
])

# Connectivity: 3 nodes per triangle (Counter-Clockwise order)
elements = np.array([
    [0, 1, 2],  # Element 0: Bottom right triangle
    [0, 2, 3]   # Element 1: Top left triangle
])

# Material Properties (Steel)
E = 200e9      # Young's Modulus (N/m^2)
nu = 0.3       # Poisson's ratio
t = 0.01       # Thickness of the plate (m)

num_nodes = len(nodes)
num_dofs = 2 * num_nodes  # u_x and u_y for each node

print(f"--- 2D CST FEM Solver ({len(elements)} Elements, {num_dofs} DOFs) ---")

# ==========================================
# 2. CONSTITUTIVE MATRIX (Plane Stress)
# ==========================================
# Relates Stress to Strain: sigma = D * epsilon
factor = E / (1 - nu**2)
D = factor * np.array([
    [1.0, nu,  0.0],
    [nu,  1.0, 0.0],
    [0.0, 0.0, (1 - nu) / 2.0]
])

# ==========================================
# 3. GLOBAL MATRICES & ASSEMBLY
# ==========================================
K_global = np.zeros((num_dofs, num_dofs))
F_global = np.zeros((num_dofs, 1))

# Apply 50,000 N pull to the right edge (Nodes 1 and 2, X-direction)
F_global[2*1, 0] = 50000.0  # Node 1 X-DOF
F_global[2*2, 0] = 50000.0  # Node 2 X-DOF

# Store B matrices and Areas to calculate stress later
B_matrices = []
Areas = []

for e, elem_nodes in enumerate(elements):
    # Get coordinates of the 3 nodes
    x1, y1 = nodes[elem_nodes[0]]
    x2, y2 = nodes[elem_nodes[1]]
    x3, y3 = nodes[elem_nodes[2]]
    
    # Calculate geometric constants (beta and gamma)
    b1, g1 = y2 - y3, x3 - x2
    b2, g2 = y3 - y1, x1 - x3
    b3, g3 = y1 - y2, x2 - x1
    
    # Calculate Area of the triangle (Determinant method)
    # Area = 0.5 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
    Area = 0.5 * (x1*b1 + x2*b2 + x3*b3)
    Areas.append(Area)
    
    # Construct the Strain-Displacement Matrix (B)
    # Size: 3 strains (ex, ey, gxy) x 6 DOFs (u1, v1, u2, v2, u3, v3)
    B = (1.0 / (2.0 * Area)) * np.array([
        [b1,  0, b2,  0, b3,  0],
        [ 0, g1,  0, g2,  0, g3],
        [g1, b1, g2, b2, g3, b3]
    ])
    B_matrices.append(B)
    
    # Calculate Element Stiffness Matrix: k = t * Area * (B^T * D * B)
    # Because B and D are constant, the integral over volume just becomes (t * Area)
    K_e = t * Area * np.dot(B.T, np.dot(D, B))
    
    # Assemble into Global Matrix
    # Get the 6 global DOF indices for these 3 nodes
    dofs = []
    for n in elem_nodes:
        dofs.extend([2*n, 2*n+1])
    
    K_global[np.ix_(dofs, dofs)] += K_e

# ==========================================
# 4. BOUNDARY CONDITIONS (Partitioning)
# ==========================================
# Fix the left wall (Node 0 and Node 3 are clamped)
fixed_dofs = [0, 1, 6, 7]  # DOFs for Node 0 (X,Y) and Node 3 (X,Y)
free_dofs = [2, 3, 4, 5]   # DOFs for Node 1 and Node 2

# Partition the Matrix
K_ff = K_global[np.ix_(free_dofs, free_dofs)]
F_f  = F_global[free_dofs]

# ==========================================
# 5. SOLVE
# ==========================================
u_global = np.zeros((num_dofs, 1))

# Solve K_ff * u_f = F_f
u_free = np.linalg.solve(K_ff, F_f)
u_global[free_dofs] = u_free

print("\nNodal Displacements (meters):")
for i in range(num_nodes):
    print(f"Node {i}: u_x = {u_global[2*i, 0]:.4e}, u_y = {u_global[2*i+1, 0]:.4e}")

# ==========================================
# 6. POST-PROCESSING: Calculate Element Stresses
# ==========================================
print("\nElement Stresses (MPa):")
# Stress = D * B * u_element
for e, elem_nodes in enumerate(elements):
    # Extract the 6 displacements for this specific element
    dofs = []
    for n in elem_nodes:
        dofs.extend([2*n, 2*n+1])
    u_e = u_global[dofs]
    
    # Calculate Strain (epsilon = B * u)
    strain = np.dot(B_matrices[e], u_e)
    
    # Calculate Stress (sigma = D * strain)
    stress = np.dot(D, strain)
    
    # Convert Pa to MPa for readable output
    stress_MPa = stress / 1e6 
    
    print(f"Element {e}:")
    print(f"  Sigma_xx = {stress_MPa[0,0]:.2f} MPa")
    print(f"  Sigma_yy = {stress_MPa[1,0]:.2f} MPa")
    print(f"  Tau_xy   = {stress_MPa[2,0]:.2f} MPa")
