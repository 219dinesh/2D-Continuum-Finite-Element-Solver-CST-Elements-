# 2D Continuum Finite Element Solver (CST Elements)

A fully-commented Python script for performing 2D Finite Element Analysis (FEA) on solid planar structures using Constant Strain Triangle (CST) elements.

Unlike 1D bars or 2D trusses (which only transfer axial loads), continuum elements model solid continuous material. This code demonstrates how to solve a 2D solid mechanics problem assuming **Plane Stress** conditions, calculating not only nodal displacements but also the internal stress tensor ($\sigma_{xx}, \sigma_{yy}, \tau_{xy}$) for each element.

## Features

* **Plane Stress Formulation:** Constructs the constitutive material matrix ($\mathbf{D}$) for thin plates.
* **CST Shape Functions:** Calculates the geometric Area, $\beta$, and $\gamma$ constants to build the Strain-Displacement matrix ($\mathbf{B}$).
* **Stiffness Integration:** Derives the element stiffness matrix for a continuum element without requiring complex numerical integration, taking advantage of the "constant strain" property of the CST.
* **Stress Recovery:** Post-processes the solved global displacements back into local element strains and stresses (MPa).
* **Zero Dependencies:** Runs entirely on standard linear algebra using `numpy`.

## The Mathematics of the CST Element

The Constant Strain Triangle is the simplest 2D solid element. It has 3 nodes, each with 2 Degrees of Freedom ($X$ and $Y$), meaning each element has a $6 \times 6$ stiffness matrix. Because the edges remain straight during deformation, the strain inside the triangle is uniform (constant) throughout its area.

### 1. The Constitutive Matrix ($\mathbf{D}$)
Hooke's law in 2D relates the stress vector $\boldsymbol{\sigma} = [\sigma_{xx}, \sigma_{yy}, \tau_{xy}]^T$ to the strain vector $\boldsymbol{\epsilon} = [\epsilon_{xx}, \epsilon_{yy}, \gamma_{xy}]^T$ via the material property matrix $\mathbf{D}$.

For a thin plate where out-of-plane stresses are zero (**Plane Stress**), the matrix is defined using Young's Modulus ($E$) and Poisson's ratio ($\nu$):

$$ \mathbf{D} = \frac{E}{1 - \nu^2} \begin{bmatrix} 1 & \nu & 0 \\ 
\nu & 1 & 0 \\ 0 & 0 & \frac{1 - \nu}{2} \end{bmatrix} $$

### 2. The Strain-Displacement Matrix ($\mathbf{B}$)
The $\mathbf{B}$ matrix relates the nodal displacements $\mathbf{u}_e = [u_1, v_1, u_2, v_2, u_3, v_3]^T$ directly to the internal strains: $\boldsymbol{\epsilon} = \mathbf{B} \mathbf{u}_e$.

For a CST element, $\mathbf{B}$ is a $3 \times 6$ matrix derived from the geometric coordinates of the three nodes $(x_1, y_1), (x_2, y_2), (x_3, y_3)$. We define geometric constants:
* $\beta_1 = y_2 - y_3$, $\quad \beta_2 = y_3 - y_1$, $\quad \beta_3 = y_1 - y_2$
* $\gamma_1 = x_3 - x_2$, $\quad \gamma_2 = x_1 - x_3$, $\quad \gamma_3 = x_2 - x_1$

Using the area of the triangle ($A$), the $\mathbf{B}$ matrix is assembled as:

$$ \mathbf{B} = \frac{1}{2A} \begin{bmatrix} \beta_1 & 0 & \beta_2 & 0 & \beta_3 & 0 \\ 
0 & \gamma_1 & 0 & \gamma_2 & 0 & \gamma_3 \\ 
\gamma_1 & \beta_1 & \gamma_2 & \beta_2 & \gamma_3 & \beta_3 \end{bmatrix} $$

### 3. The Element Stiffness Matrix ($\mathbf{K}_e$)
The general formula for element stiffness is the volume integral $\int \mathbf{B}^T \mathbf{D} \mathbf{B} dV$. 

Because $\mathbf{B}$ is comprised entirely of geometric constants (thanks to the linear shape functions of the triangle), the matrix product $\mathbf{B}^T \mathbf{D} \mathbf{B}$ is also a constant. Therefore, the integral simply evaluates to multiplying by the element's volume (Thickness $t \times$ Area $A$):

$$ \mathbf{K}_e = t A \mathbf{B}^T \mathbf{D} \mathbf{B} $$

### 4. Global Assembly and Solving
Similar to 1D problems, the $6 \times 6$ local matrices are added ("stamped") into a larger Global Matrix ($\mathbf{K}_{global}$). After applying boundary conditions (partitioning out the fixed nodes where displacement is 0), we solve the linear system for the unknown displacements:

$$ \mathbf{U}_{free} = \mathbf{K}_{ff}^{-1} \mathbf{F}_{free} $$

### 5. Stress Recovery
Once all nodal displacements are solved, we extract the 6 specific displacements for a given element ($\mathbf{u}_e$) and calculate the internal stresses. 

$$ \boldsymbol{\sigma} = \mathbf{D} \mathbf{B} \mathbf{u}_e $$

This yields the exact stress components $[\sigma_{xx}, \sigma_{yy}, \tau_{xy}]$ for that triangle.

## Prerequisites
To run this script, you only need Python and the NumPy library installed.

```Bash
pip install numpy
```
## Usage 

Run the script directly from your terminal:
```bash
python 2d_cst_fem.py
```
## Example Output
Based on the default configuration (a 2m x 1m steel plate modeled with 2 triangles, clamped on the left, pulled with 100kN total force on the right), you will see:

```Plaintext
--- 2D CST FEM Solver (2 Elements, 8 DOFs) ---

Nodal Displacements (meters):
Node 0: u_x = 0.0000e+00, u_y = 0.0000e+00
Node 1: u_x = 1.0000e-04, u_y = 7.5000e-06
Node 2: u_x = 1.0000e-04, u_y = -7.5000e-06
Node 3: u_x = 0.0000e+00, u_y = 0.0000e+00

Element Stresses (MPa):
Element 0:
  Sigma_xx = 10.00 MPa
  Sigma_yy = 0.00 MPa
  Tau_xy   = 0.00 MPa
Element 1:
  Sigma_xx = 10.00 MPa
  Sigma_yy = -0.00 MPa
  Tau_xy   = 0.00 MPa
```
Notice how the Poisson effect is visible in the nodal displacements: as the plate stretches in the X-direction (1.00e-04), it symmetrically contracts in the Y-direction (+7.50e-06 at the bottom right, -7.50e-06 at the top right)!
