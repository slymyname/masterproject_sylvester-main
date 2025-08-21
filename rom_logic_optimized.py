import numba
import numpy as np

@numba.jit(nopython=True, cache=True)
def _average_properties_numba(prop_nodes, startnode, targetnode):
    """
    Numba-jitted function to average nodal properties onto edges (pipes, etc.).
    """
    num_edges = len(startnode)
    prop_edges = np.empty(num_edges, dtype=np.float64)
    for i in range(num_edges):
        prop_edges[i] = 0.5 * (prop_nodes[startnode[i]] + prop_nodes[targetnode[i]])
    return prop_edges

@numba.jit(nopython=True, cache=True)
def _calc_lambda_numba(Re, k, d, pipes_lambda_in):
    """
    Numba-jitted function to calculate the pipe friction factor lambda.
    Uses the Haaland approximation for turbulent flow.
    """
    pipes_lambda_out = pipes_lambda_in.copy()
    for i in range(len(Re)):
        if Re[i] <= 2320:
            if Re[i] > 0:
                pipes_lambda_out[i] = 64.0 / Re[i]
            else:
                pipes_lambda_out[i] = 1.0 # fallback for zero flow
        else:
            # Haaland approximation for turbulent flow
            inv_sqrt_lambda = -1.8 * np.log10((k[i] / (3.7 * d[i]))**1.11 + 6.9 / Re[i])
            pipes_lambda_out[i] = (1.0 / inv_sqrt_lambda)**2
    
    # Clip to handle potential division by zero for laminar flow Re=0
    for i in range(len(pipes_lambda_out)):
        if pipes_lambda_out[i] > 1.0:
             pipes_lambda_out[i] = 1.0
             
    return pipes_lambda_out

@numba.jit(nopython=True, cache=True)
def _build_and_project_thermic_system_numba(Sa_data, Sa_indices, Sa_indptr,
                                           Se_data, Se_indices, Se_indptr,
                                           V_r, W, C_diag, D):
    """
    Builds the full thermal system matrices and projects them to the reduced
    space using efficient, matrix-free operations.

    Args:
        Sa_data, Sa_indices, Sa_indptr: CSR components of the Sa matrix.
        Se_data, Se_indices, Se_indptr: CSR components of the Se matrix.
        V_r (N x r): POD basis matrix.
        W (M,): Advection weights vector.
        C_diag (M,): Diagonal of the C matrix.
        D (M,): Source/boundary vector.

    Returns:
        T_A_r (r x r): Reduced system matrix.
        T_B_r (r,): Reduced source vector.
    """
    M, N = Sa_indptr.shape[0] - 1, V_r.shape[0]
    r = V_r.shape[1]

    # VSa = Sa @ V_r
    # VSe = Se @ V_r
    VSa = np.zeros((M, r))
    VSe = np.zeros((M, r))

    for i in range(M):
        for j_ptr in range(Sa_indptr[i], Sa_indptr[i+1]):
            j = Sa_indices[j_ptr]
            for k in range(r):
                VSa[i, k] += V_r[j, k] # Sa_data is all 1s

    for i in range(M):
        for j_ptr in range(Se_indptr[i], Se_indptr[i+1]):
            j = Se_indices[j_ptr]
            for k in range(r):
                VSe[i, k] += V_r[j, k] # Se_data is all 1s

    # Term 1 of T_A_r: VSa.T @ diag(W * C_diag) @ VSe
    Term1_T_A_r = np.zeros((r, r))
    WC_diag = W * C_diag
    for i in range(r):
        for j in range(r):
            sum_val = 0.0
            for k in range(M):
                sum_val += VSa[k, i] * WC_diag[k] * VSe[k, j]
            Term1_T_A_r[i, j] = sum_val
            
    # Term 2 of T_A_r: VSe.T @ diag(W) @ VSe
    Term2_T_A_r = np.zeros((r, r))
    for i in range(r):
        for j in range(r):
            sum_val = 0.0
            for k in range(M):
                sum_val += VSe[k, i] * W[k] * VSe[k, j]
            Term2_T_A_r[i, j] = sum_val

    T_A_r = Term1_T_A_r - Term2_T_A_r

    # T_B_r = -VSa.T @ (W * D)
    WD_vec = W * D
    T_B_r = np.zeros(r)
    for i in range(r):
        sum_val = 0.0
        for k in range(M):
            sum_val += VSa[k, i] * WD_vec[k]
        T_B_r[i] = -sum_val

    return T_A_r, T_B_r
