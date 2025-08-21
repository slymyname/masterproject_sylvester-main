import numpy as np
import os
import matplotlib.pyplot as plt

def perform_pod_analysis():
    fom_snapshot_dir = 'fom_snapshots'
    t_snapshots_path = os.path.join(fom_snapshot_dir, 'fom_T_snapshots.npy')

    if not os.path.exists(t_snapshots_path):
        print(f"ERROR: Snapshot file not found: {t_snapshots_path}")
        print("Please ensure you have run Main.py to generate snapshots.")
        return

    try:
        # Load the temperature snapshots (N x N_s)
        # Each column is a snapshot, each row is a node
        T_snapshots = np.load(t_snapshots_path)
        print(f"Successfully loaded {t_snapshots_path}")
        print(f"Shape of T_snapshots: {T_snapshots.shape}")

        if T_snapshots.size == 0:
            print("T_snapshots array is empty. Cannot perform POD.")
            return
        
        num_nodes, num_snapshots = T_snapshots.shape
        if num_snapshots < 1:
            print("Not enough snapshots to perform POD.")
            return

        # Perform Singular Value Decomposition (SVD)
        # S = U * Sigma * V_transpose
        # U contains the POD modes (spatial basis vectors)
        # Sigma contains the singular values
        print("\nPerforming SVD...")
        U, s, Vt = np.linalg.svd(T_snapshots, full_matrices=False) # full_matrices=False for thin SVD
        print("SVD complete.")

        print(f"Shape of U (POD modes): {U.shape}")
        print(f"Shape of s (singular values): {s.shape}")
        print(f"Shape of Vt (V transpose): {Vt.shape}")

        # Singular values indicate the energy of each mode
        print("\nSingular values (s):")
        print(s)

        # Calculate the energy captured by the modes
        squared_singular_values = s**2
        total_energy = np.sum(squared_singular_values)

        print(f"\nTotal energy (sum of squared singular values): {total_energy:.4e}")
        
        cumulative_energy = np.cumsum(squared_singular_values) / total_energy

        print("\nCumulative energy captured by retaining k modes:")
        for i, energy_fraction in enumerate(cumulative_energy):
            print(f"  {i+1} mode(s): {energy_fraction*100:.4f}%")

        # Save POD results for later use
        u_matrix_path = os.path.join(fom_snapshot_dir, 'pod_U_matrix.npy')
        singular_values_path = os.path.join(fom_snapshot_dir, 'pod_singular_values.npy')
        
        np.save(u_matrix_path, U)
        print(f"\nSaved U matrix (all potential POD modes) to: {u_matrix_path}")
        np.save(singular_values_path, s)
        print(f"Saved singular values to: {singular_values_path}")

        #Plotting Singular Values
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.semilogy(np.arange(1, len(s) + 1), s, 'o-', markersize=5)
        plt.title('Singular Values (Log Scale)')
        plt.xlabel('Mode Number')
        plt.ylabel('Singular Value')
        plt.grid(True, which="both", ls="--")
        plt.xticks(np.arange(1, len(s) + 1))

        plt.subplot(1, 2, 2)
        plt.plot(np.arange(1, len(cumulative_energy) + 1), cumulative_energy * 100, 'o-', markersize=5)
        plt.title('Cumulative Energy Spectrum')
        plt.xlabel('Number of Modes Retained')
        plt.ylabel('Cumulative Energy (%)')
        plt.grid(True, which="both", ls="--")
        plt.xticks(np.arange(1, len(cumulative_energy) + 1))
        plt.yticks(np.arange(0, 101, 10))
        
        plt.tight_layout()
        plot_path = os.path.join(fom_snapshot_dir, 'pod_singular_value_analysis.png')
        plt.savefig(plot_path)
        print(f"Saved singular value plots to: {plot_path}")
        

    except FileNotFoundError:
        print(f"ERROR: Snapshot file not found at {t_snapshots_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    perform_pod_analysis() 