import numpy as np
import os

def create_pod_basis(r_value):
    fom_snapshot_dir = 'fom_snapshots'
    u_matrix_path = os.path.join(fom_snapshot_dir, 'pod_U_matrix.npy')
    vr_basis_path = os.path.join(fom_snapshot_dir, f'pod_Vr_basis_r{r_value}.npy')

    if not os.path.exists(u_matrix_path):
        print(f"ERROR: U matrix file not found: {u_matrix_path}")
        print("Please ensure you have run perform_pod_analysis.py first.")
        return

    try:
        U = np.load(u_matrix_path)
        print(f"Successfully loaded {u_matrix_path}")
        print(f"Shape of full U matrix: {U.shape}")

        if r_value > U.shape[1]:
            print(f"ERROR: r_value ({r_value}) is greater than the number of available modes ({U.shape[1]}).")
            print(f"Please choose an r_value <= {U.shape[1]}.")
            return
        
        # Extract the first r_value columns to form the POD basis V_r
        V_r = U[:, :r_value]
        print(f"\nSelected r = {r_value} modes.")
        print(f"Shape of reduced POD basis V_r: {V_r.shape}")

        np.save(vr_basis_path, V_r)
        print(f"Saved reduced POD basis V_r to: {vr_basis_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    # choose the number of modes (r)
    r_modes_to_keep = 2 
    
    
    print(f"Creating POD basis with r = {r_modes_to_keep} modes...")
    create_pod_basis(r_modes_to_keep) 