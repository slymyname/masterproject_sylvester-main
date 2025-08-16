import numpy as np
import os

def select_deim_points(num_interpolation_points, basis_matrix_path, output_indices_path, output_proj_matrix_path, data_label="Data"):
    """
    Selects DEIM interpolation points using the greedy algorithm.

    Args:
        num_interpolation_points (int): The number of DEIM points to select (m).
        basis_matrix_path (str): Path to the .npy file containing the DEIM basis matrix U (num_elements x num_modes).
        output_indices_path (str): Path to save the selected DEIM indices.
        output_proj_matrix_path (str): Path to save the DEIM projection matrix.
        data_label (str): Label for the data being processed, for print statements.
    """
    print(f"--- Selecting DEIM points for {data_label} (m = {num_interpolation_points}) ---")
    if not os.path.exists(basis_matrix_path):
        print(f"ERROR: DEIM basis matrix file not found: {basis_matrix_path}")
        return False

    try:
        U_full = np.load(basis_matrix_path)
        print(f"Successfully loaded DEIM basis U_full from {basis_matrix_path}")
        print(f"Shape of U_full matrix: {U_full.shape}")

        num_elements, num_available_modes = U_full.shape

        if num_interpolation_points > num_available_modes:
            print(f"ERROR: num_interpolation_points ({num_interpolation_points}) for {data_label} is greater than the number of available modes ({num_available_modes}).")
            return False
        if num_interpolation_points <= 0:
            print(f"ERROR: num_interpolation_points for {data_label} must be positive.")
            return False

        U_m = U_full[:, :num_interpolation_points]
        print(f"Using first {num_interpolation_points} modes for DEIM basis U_m for {data_label}. Shape: {U_m.shape}")

        deim_indices = []
        idx = np.argmax(np.abs(U_m[:, 0]))
        deim_indices.append(idx)

        for j in range(1, num_interpolation_points):
            P_T_U_cols_1_to_j = U_m[deim_indices, :j]
            P_T_U_col_j_plus_1 = U_m[deim_indices, j]
            try:
                c, _, _, _ = np.linalg.lstsq(P_T_U_cols_1_to_j, P_T_U_col_j_plus_1, rcond=None)
            except np.linalg.LinAlgError as e:
                print(f"Linear algebra error for {data_label} at point {j+1}: {e}")
                return False
            residual_vector = U_m[:, j] - (U_m[:, :j] @ c)
            new_idx = np.argmax(np.abs(residual_vector))
            if new_idx in deim_indices:
                 # If the argmax is already chosen, find the next best unique one
                 sorted_residual_indices = np.argsort(np.abs(residual_vector))[::-1]
                 for sorted_idx in sorted_residual_indices:
                     if sorted_idx not in deim_indices:
                         new_idx = sorted_idx
                         break
                 else: # All remaining points have zero residual or are already picked
                     print(f"Warning: Could not find a unique new DEIM index for {data_label} at point {j+1}. Residuals might be zero or all high-value points taken.")
                     # Fallback: pick the first available unique index if any, otherwise break or handle error
                     available_indices = [i for i in range(num_elements) if i not in deim_indices]
                     if available_indices:
                         new_idx = available_indices[0]
                         print(f"Fallback: picking first available unique index {new_idx}")
                     else:
                         print(f"Error: No more unique indices to pick for {data_label}.")
                         return False 
            deim_indices.append(new_idx)

        deim_indices = np.array(deim_indices)
        print(f"\nSelected DEIM interpolation indices for {data_label} (0-indexed):\n{deim_indices}")
        np.save(output_indices_path, deim_indices)
        print(f"Saved DEIM indices for {data_label} to: {output_indices_path}")

        P_T_U_m_matrix = U_m[deim_indices, :] 
        if P_T_U_m_matrix.shape[0] != P_T_U_m_matrix.shape[1]:
            print(f"ERROR: Matrix U_m[deim_indices,:] (P^T U_m) for {data_label} is not square ({P_T_U_m_matrix.shape}).")
            return False
        try:
            inv_P_T_U_m_matrix = np.linalg.inv(P_T_U_m_matrix)
            deim_projection_matrix = U_m @ inv_P_T_U_m_matrix
            np.save(output_proj_matrix_path, deim_projection_matrix)
            print(f"Saved DEIM projection matrix for {data_label} ({deim_projection_matrix.shape}) to: {output_proj_matrix_path}")
        except np.linalg.LinAlgError as e:
            print(f"Error inverting P_T_U_m_matrix or forming projection matrix for {data_label}: {e}")
            return False
        print(f"--- DEIM point selection for {data_label} complete. ---")
        return True

    except Exception as e:
        print(f"An error occurred in select_deim_points for {data_label}: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    fom_snapshot_dir = 'fom_snapshots'
    
    # --- Parameters for Q_loss_Pipe DEIM --- #
    m_q_loss = 2
    q_loss_basis_matrix_path = os.path.join(fom_snapshot_dir, 'deim_q_loss_U_matrix.npy')
    q_loss_indices_path = os.path.join(fom_snapshot_dir, f'deim_q_loss_indices_m{m_q_loss}.npy')
    q_loss_proj_matrix_path = os.path.join(fom_snapshot_dir, f'deim_q_loss_proj_matrix_m{m_q_loss}.npy')
    print(f"\n--- Configuring DEIM point selection for Q_loss_Pipe with m = {m_q_loss} ---")
    select_deim_points(
        num_interpolation_points=m_q_loss,
        basis_matrix_path=q_loss_basis_matrix_path,
        output_indices_path=q_loss_indices_path,
        output_proj_matrix_path=q_loss_proj_matrix_path,
        data_label="Q_loss_Pipe"
    )
    print("-----------------------------------------------------\n")

    # --- Parameters for rho_Nodes DEIM --- #
    m_rho = 2
    rho_basis_matrix_path = os.path.join(fom_snapshot_dir, 'deim_rho_U_matrix.npy')
    rho_indices_path = os.path.join(fom_snapshot_dir, f'deim_rho_indices_m{m_rho}.npy')
    rho_proj_matrix_path = os.path.join(fom_snapshot_dir, f'deim_rho_proj_matrix_m{m_rho}.npy')
    print(f"\n--- Configuring DEIM point selection for rho_Nodes with m = {m_rho} ---")
    select_deim_points(
        num_interpolation_points=m_rho,
        basis_matrix_path=rho_basis_matrix_path,
        output_indices_path=rho_indices_path,
        output_proj_matrix_path=rho_proj_matrix_path,
        data_label="rho_Nodes"
    )
    print("-----------------------------------------------------\n")

    # --- Parameters for cp_Nodes DEIM --- #
    m_cp = 2
    cp_basis_matrix_path = os.path.join(fom_snapshot_dir, 'deim_cp_U_matrix.npy')
    cp_indices_path = os.path.join(fom_snapshot_dir, f'deim_cp_indices_m{m_cp}.npy')
    cp_proj_matrix_path = os.path.join(fom_snapshot_dir, f'deim_cp_proj_matrix_m{m_cp}.npy')
    print(f"\n--- Configuring DEIM point selection for cp_Nodes with m = {m_cp} ---")
    select_deim_points(
        num_interpolation_points=m_cp,
        basis_matrix_path=cp_basis_matrix_path,
        output_indices_path=cp_indices_path,
        output_proj_matrix_path=cp_proj_matrix_path,
        data_label="cp_Nodes"
    )
    print("-----------------------------------------------------\n")

    # --- Parameters for viscosity_Nodes DEIM --- #
    m_viscosity = 2
    viscosity_basis_matrix_path = os.path.join(fom_snapshot_dir, 'deim_viscosity_U_matrix.npy')
    viscosity_indices_path = os.path.join(fom_snapshot_dir, f'deim_viscosity_indices_m{m_viscosity}.npy')
    viscosity_proj_matrix_path = os.path.join(fom_snapshot_dir, f'deim_viscosity_proj_matrix_m{m_viscosity}.npy')
    print(f"\n--- Configuring DEIM point selection for viscosity_Nodes with m = {m_viscosity} ---")
    select_deim_points(
        num_interpolation_points=m_viscosity,
        basis_matrix_path=viscosity_basis_matrix_path,
        output_indices_path=viscosity_indices_path,
        output_proj_matrix_path=viscosity_proj_matrix_path,
        data_label="viscosity_Nodes"
    )
    print("-----------------------------------------------------\n") 