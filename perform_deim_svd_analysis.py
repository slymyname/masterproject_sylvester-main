import numpy as np
import os
import matplotlib.pyplot as plt

def perform_deim_svd_analysis(snapshot_file_name, output_prefix, data_label):
    fom_snapshot_dir = 'fom_snapshots'
    snapshots_path = os.path.join(fom_snapshot_dir, snapshot_file_name)

    if not os.path.exists(snapshots_path):
        print(f"ERROR: Snapshot file not found: {snapshots_path}")
        return

    try:
        snapshots = np.load(snapshots_path)
        print(f"Successfully loaded {snapshots_path}")
        print(f"Shape of {data_label} snapshots: {snapshots.shape}")

        if snapshots.size == 0:
            print(f"{data_label} snapshots array is empty. Cannot perform SVD analysis.")
            return
        
        num_elements, num_s = snapshots.shape
        if num_s < 1:
            print("Not enough snapshots to perform SVD.")
            return

        print(f"\nPerforming SVD on {data_label} snapshots...")
        U, s, Vt = np.linalg.svd(snapshots, full_matrices=False)
        print("SVD complete.")

        print(f"Shape of U (basis for {data_label}): {U.shape}")
        print(f"Shape of s (singular values for {data_label}): {s.shape}")

        print(f"\nSingular values for {data_label}:")
        print(s)

        squared_singular_values = s**2
        total_energy = np.sum(squared_singular_values)
        print(f"\nTotal energy for {data_label}: {total_energy:.4e}")
        
        cumulative_energy = np.cumsum(squared_singular_values) / total_energy
        print(f"\nCumulative energy captured for {data_label}:")
        for i, energy_fraction in enumerate(cumulative_energy):
            print(f"  {i+1} mode(s): {energy_fraction*100:.4f}%")

        u_matrix_path = os.path.join(fom_snapshot_dir, f'{output_prefix}_U_matrix.npy')
        singular_values_path = os.path.join(fom_snapshot_dir, f'{output_prefix}_s_values.npy')
        
        np.save(u_matrix_path, U)
        print(f"\nSaved U matrix for {data_label} to: {u_matrix_path}")
        np.save(singular_values_path, s)
        print(f"Saved singular values for {data_label} to: {singular_values_path}")

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.semilogy(np.arange(1, len(s) + 1), s, 'o-', markersize=5)
        plt.title(f'Singular Values for {data_label} (Log Scale)')
        plt.xlabel('Mode Number'); plt.ylabel('Singular Value')
        plt.grid(True, which="both", ls="--"); plt.xticks(np.arange(1, len(s) + 1))

        plt.subplot(1, 2, 2)
        plt.plot(np.arange(1, len(cumulative_energy) + 1), cumulative_energy * 100, 'o-', markersize=5)
        plt.title(f'Cumulative Energy for {data_label}')
        plt.xlabel('Number of Modes Retained'); plt.ylabel('Cumulative Energy (%)')
        plt.grid(True, which="both", ls="--"); plt.xticks(np.arange(1, len(cumulative_energy) + 1))
        plt.yticks(np.arange(0, 101, 10))
        
        plt.tight_layout()
        plot_path = os.path.join(fom_snapshot_dir, f'{output_prefix}_svalue_analysis.png')
        plt.savefig(plot_path)
        print(f"Saved {data_label} singular value plots to: {plot_path}")

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    # Analysis for rho_Nodes 
    print("--- Starting DEIM SVD Analysis for rho_Nodes ---")
    perform_deim_svd_analysis(
        snapshot_file_name='fom_rho_nodes_snapshots.npy',
        output_prefix='deim_rho',
        data_label='rho_Nodes'
    )
    print("--------------------------------------------------\n")

    # Analysis for cp_Nodes
    print("--- Starting DEIM SVD Analysis for cp_Nodes ---")
    perform_deim_svd_analysis(
        snapshot_file_name='fom_cp_nodes_snapshots.npy',
        output_prefix='deim_cp',
        data_label='cp_Nodes'
    )
    print("--------------------------------------------------\n")

    # Analysis for viscosity_Nodes
    print("--- Starting DEIM SVD Analysis for viscosity_Nodes ---")
    perform_deim_svd_analysis(
        snapshot_file_name='fom_viscosity_nodes_snapshots.npy',
        output_prefix='deim_viscosity',
        data_label='viscosity_Nodes'
    )
    print("--------------------------------------------------\n") 