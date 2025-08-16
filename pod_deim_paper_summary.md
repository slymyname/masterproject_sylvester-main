Summary of "POD-DEIM model order reduction technique for model predictive control in continuous chemical processing"

Paper Title: POD-DEIM model order reduction technique for model predictive control in continuous chemical processing
Authors: Van Bo Nguyen, Si Bui Quang Tran, Saif A. Khan, Jiawei Rong, Jing Lou
Source: Computers and Chemical Engineering 133 (2020) 106638

1. Abstract Summary

The paper proposes a model order reduction (MOR) technique using Proper Orthogonal Decomposition (POD) and the Discrete Empirical Interpolation Method (DEIM) to address the computational challenges of large-scale, nonlinear models in chemical process control, particularly for Model Predictive Control (MPC).
POD is used to project the original state variables onto a lower-dimensional subspace. DEIM is used to approximate the projected nonlinear terms by evaluating the original nonlinear function at a few selected interpolation points. This avoids the full computation of the nonlinear part of the model.
The combined POD-DEIM approach significantly reduces the computational cost of the simulation, making it suitable for real-time MPC applications. The effectiveness is demonstrated on a multi-scale chemical reactor model.

2. Core Problem Addressed

Large-scale, first-principle models in chemical engineering are often described by systems of nonlinear ordinary differential equations (ODEs) or partial differential equations (PDEs) after spatial discretization. These models can be computationally expensive, making their use in demanding applications like real-time optimization and control (e.g., MPC) challenging. MPC requires repeated simulations of the process model over a prediction horizon.

3. Proposed Solution: POD-DEIM

The paper employs a two-step MOR technique:

1.  POD:To find a low-dimensional representation of the system's state variables.
2.  DEIM: To efficiently approximate the nonlinear function in the reduced system, which arises from the original nonlinearities.

This combination allows for substantial computational savings while retaining essential model dynamics.

4. Key Methodologies Explained

4.1. Full Order Model (FOM)

The FOM is typically a system of nonlinear differential equations:
\[ \\dot{\\mathbf{x}}(t) = \\mathbf{f}(\\mathbf{x}(t), \\mathbf{u}(t), t) + \\mathbf{B}\\mathbf{u}(t) \]
\[ \\mathbf{y}(t) = \\mathbf{C}\\mathbf{x}(t) \]
(Adapted from general forms, paper uses specific forms like Eq. 2.1-2.6 for their MPC problem, e.g., \( \\dot{\\mathbf{x}} = A\\mathbf{x} + \\mathbf{f}(\\mathbf{x}) + Bu \\)).
For our project, \( \\mathbf{x} \) would primarily be the temperature vector \( \\mathbf{T} \). The nonlinear function \( \\mathbf{f}(\\mathbf{x}) \) would encompass terms dependent on \( \\mathbf{T} \), like fluid properties (density \( \\rho \), specific heat \( c_p \), viscosity \( \mu \)) and temperature-dependent heat losses.

4.2. Proper Orthogonal Decomposition (POD)

Purpose: To find an optimal set of basis vectors that can represent the system's state with a much lower dimension.
Process:
    1.  Snapshot Collection: Simulate the FOM for various inputs/conditions and collect snapshots of the state vector \( \\mathbf{x} \) at different time points: \( S = [\\mathbf{x}(t_1), \\mathbf{x}(t_2), ..., \\mathbf{x}(t_{N_s})] \).
         In our project: These are the `T_state_snapshots.npy` we are generating.
    2.  Singular Value Decomposition (SVD): Perform SVD on the snapshot matrix \( S \): \( S = U \\Sigma V^T \).
    3.  Basis Generation: The POD basis \( \\mathbf{V}_r \) (or \( \\Phi \) in the paper, denoted as \( V_x \) in their Fig 1) consists of the first \( r \) left singular vectors (columns of \( U \)) corresponding to the largest singular values. The dimension \( r \\ll N \), where \( N \) is the dimension of the full state vector.
State Approximation: The full state \( \\mathbf{x} \) is approximated by a linear combination of these basis vectors:
    \[ \\mathbf{x}(t) \\approx \\mathbf{V}_r \\mathbf{x}_r(t) \]
    (Similar to paper's Eq. 3.4: \( x = V \\tilde{x} \), where \( V \) is \( \\mathbf{V}_r \) and \( \\tilde{x} \) is \( \mathbf{x}_r \)).

4.3. Discrete Empirical Interpolation Method (DEIM)

Purpose: To approximate the (potentially complex and high-dimensional) nonlinear function \( \\mathbf{f}(\\mathbf{x}) \) efficiently. If POD is applied directly, the nonlinear term becomes \( \\mathbf{V}_r^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) \), which still requires evaluating the full \( N \)-dimensional \( \\mathbf{f} \). DEIM avoids this.
Process:
    1.  Nonlinear Function Snapshot Collection: Collect snapshots of the nonlinear function \( \\mathbf{f} \) evaluated at the state snapshots: \( S_f = [\\mathbf{f}(\\mathbf{x}(t_1)), \\mathbf{f}(\\mathbf{x}(t_2)), ..., \\mathbf{f}(\\mathbf{x}(t_{N_s}))] \).
        *   In our project: This corresponds to why we are collecting snapshots of `rho_Nodes`, `cp_Nodes`, `viscosity_Nodes`, and `Q_loss_Pipe`, as these are the components that make up the nonlinearities in our thermal equations. The actual \( \mathbf{f} \) for DEIM would be the assembled vector from the system equations that these terms contribute to (e.g., the terms in `self.T_A` and `self.T_B` that depend nonlinearly on `T`).
    2.  SVD for Nonlinear Basis: Perform SVD on \( S_f \) to obtain a basis \( \\mathbf{U}_m \) (or \( \\Psi \) in the paper, denoted as \( V_f \) in their Fig 1) for the nonlinear term, of dimension \( m \).
    3.  Interpolation Index Selection (DEIM Algorithm):
        A greedy algorithm selects \( m \) interpolation indices (row indices) \( \\wp = [p_1, ..., p_m]^T \). These are the points where the original nonlinear function \( \\mathbf{f} \) will actually be evaluated.
        The selection process (Eq. 3.8 in the paper is the core DEIM point selection algorithm) aims to pick points that best capture the behavior of the function based on its projection onto the basis \( \mathbf{U}_m \).
Nonlinear Function Approximation: The nonlinear function \( \\mathbf{f}(\\mathbf{x}) \) is approximated as:
    \[ \\mathbf{f}(\\mathbf{x}) \\approx \\mathbf{U}_m ( (\\mathbf{P}^T \\mathbf{U}_m)^{-1} (\\mathbf{P}^T \\mathbf{f}(\\mathbf{x})) ) \]
    where \( \\mathbf{P} \) is a selection matrix ( \( \\mathbf{P}^T \\mathbf{g} = [g_{p_1}, ..., g_{p_m}]^T \) extracts the elements of \( \\mathbf{g} \) at the DEIM indices).
    (Paper's Eq. 3.10: \( f(\\cdot) \\approx \\Psi c_f \), and Eq. 3.11 for \( c_f \): \( c_f = (P^T \\Psi)^{-1} P^T f(\\cdot) \)).
    This means only \( m \) components of \( \\mathbf{f}(\\mathbf{x}) \) need to be computed.

4.4. Combined POD-DEIM Reduced Order Model (ROM)

Applying POD for the state and DEIM for the nonlinear term to the FOM:
\( \\dot{\\mathbf{x}} = \\mathbf{f}(\\mathbf{x}) + \\text{other linear terms} \)
Substitute \( \\mathbf{x} \\approx \\mathbf{V}_r \\mathbf{x}_r \):
\( \\mathbf{V}_r \\dot{\\mathbf{x}}_r = \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) + ... \)
Project onto the POD basis by multiplying with \( \\mathbf{V}_r^T \):
\( \\mathbf{V}_r^T \\mathbf{V}_r \\dot{\\mathbf{x}}_r = \\mathbf{V}_r^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) + ... \)
Since \( \\mathbf{V}_r^T \\mathbf{V}_r = \\mathbf{I}_r \) (if basis vectors are orthonormal):
\( \\dot{\\mathbf{x}}_r = \\mathbf{V}_r^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) + \\mathbf{V}_r^T (\text{other linear terms involving } \\mathbf{V}_r \\mathbf{x}_r ) \)

Now, apply DEIM to \( \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) \):
\( \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) \\approx \\mathbf{U}_m ( (\\mathbf{P}^T \\mathbf{U}_m)^{-1} (\\mathbf{P}^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r)) ) \)

The ROM becomes:
\( \\dot{\\mathbf{x}}_r = \\mathbf{V}_r^T \\mathbf{U}_m ( (\\mathbf{P}^T \\mathbf{U}_m)^{-1} (\\mathbf{P}^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r)) ) + \hat{A} \\mathbf{x}_r + \hat{B} \\mathbf{u} \)
where \( \hat{A} = \\mathbf{V}_r^T A \\mathbf{V}_r \) and \( \hat{B} = \\mathbf{V}_r^T B \) if the FOM was \( \\dot{x} = Ax + f(x) + Bu \).
The key is that \( \\mathbf{P}^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) \) only requires evaluating \( \\mathbf{f} \) at \( m \) specific entries corresponding to the original full state vector, but those entries are themselves computed based on the reduced state \( \mathbf{x}_r \).

The matrices \( \\mathbf{V}_r^T \\mathbf{U}_m \), \( (\\mathbf{P}^T \\mathbf{U}_m)^{-1} \\mathbf{P}^T \), and any projected linear system matrices can be precomputed offline.

5. Workflow (Simplified for ROM construction)

1.  Offline Stage:
     Generate state snapshots \( \\mathbf{x}(t_i) \) from the FOM.
     Construct POD basis \( \\mathbf{V}_r \) from state snapshots.
     Generate nonlinear function snapshots \( \\mathbf{f}(\\mathbf{x}(t_i)) \).
     Construct DEIM basis \( \\mathbf{U}_m \) from nonlinear function snapshots.
     Determine DEIM interpolation points \( \\wp \) and matrix \( \\mathbf{P} \).
      Precompute ROM matrices: \( \\hat{A} = \\mathbf{V}_r^T A \\mathbf{V}_r \), \( \\hat{B} = \\mathbf{V}_r^T B \), \( \\mathbf{V}_r^T \\mathbf{U}_m \), \( (\\mathbf{P}^T \\mathbf{U}_m)^{-1} \).
2.  Online Stage:
     Solve the much smaller ROM:
        \( \\dot{\\mathbf{x}}_r = \\mathbf{V}_r^T \\mathbf{U}_m ( (\\mathbf{P}^T \\mathbf{U}_m)^{-1} (\\mathbf{P}^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r)) ) + \hat{A} \\mathbf{x}_r + \hat{B} \\mathbf{u} \)
      To evaluate \( \\mathbf{P}^T \\mathbf{f}(\\mathbf{V}_r \\mathbf{x}_r) \):
        1.  Reconstruct the necessary full state components for the nonlinear function evaluation at DEIM points: \( \\mathbf{x}_{DEIM\_eval} = \\mathbf{V}_r \\mathbf{x}_r \). (Only elements of \( \mathbf{x}_{DEIM\_eval} \) that are inputs to the selected \( f_i \) are needed).
        2.  Evaluate only the \(m\) rows of \( \\mathbf{f} \) at these DEIM points.
     Reconstruct the full state if needed: \( \\mathbf{x}(t) \\approx \\mathbf{V}_r \\mathbf{x}_r(t) \).

(The paper's Fig. 1 shows a workflow for MPC using POD-DEIM, which includes the ROM construction as a core part).

 6. Benefits Highlighted

Significant computational speed-up: The ROM is much smaller than the FOM.
Preservation of nonlinear characteristics: DEIM helps maintain accuracy for nonlinear systems.
Suitability for real-time applications: Essential for MPC.
Systematic approach: Provides a structured way to derive reduced models.

7. Relevance to Our Project

State Reduction (POD):Our temperature vector `TimeSeries.T` is the primary state to be reduced. POD can find a compact basis for `T`.
   Nonlinear Term Approximation (DEIM): The key nonlinearities in `TimeSeries.timestepStationary` come from:
       `waterProperties()`: `rho_Nodes`, `cp_Nodes`, `viscosity_Nodes` depend on `T` (and `p`, though `p` is solved separately). These affect the hydraulic resistance `self.R` and the matrix `self.W` in the thermal equations.
       `buildTemperatureMatrix()`: `self.Q_loss_Pipe` depends on `T` and `T_out`. The advection terms (involving `self.W`) and heat loss terms contribute to `self.T_A` and `self.T_B`.
Applying POD-DEIM:
    1.  POD Basis for T: Use `T_state_snapshots.npy`.
    2.  DEIM for Nonlinear Functions:
           We need to define precisely which parts of the equations constitute \( \mathbf{f}(\mathbf{T}) \).
          The matrices `self.T_A` and `self.T_B` in `self.calcT()` (where `T_A * T = T_B` is solved) are constructed based on `self.W` (which includes `cp_Nodes`, `rho_Nodes` via `self.cp`, `self.rho`) and `self.Q_loss_Pipe`.
          Instead of directly forming a DEIM approximation for `T_A` and `T_B` (which are matrices), DEIM is typically applied to vector-valued functions.
          One approach could be to identify the assembled right-hand side vector or the vector of influential nonlinear components that contribute to the heat balance equations before they are put into matrix form for the linear solve `K_th * T = b_th`.
         Alternatively, as the paper handles general \( f(x) \), the "nonlinear part" of the equation system `K_th(T,p) * T = b_th(T,p, T_supply, T_out)` would be the terms in `K_th` and `b_th` that depend nonlinearly on `T`.
    3.  The snapshots `deim_rho_nodes_snapshots.npy`, `deim_cp_nodes_snapshots.npy`, `deim_viscosity_nodes_snapshots.npy`, and `deim_Q_loss_pipe_snapshots.npy` represent the *ingredients* of these nonlinearities. We will likely need to form snapshots of the actual assembled nonlinear vector term from the discretized heat equations.

