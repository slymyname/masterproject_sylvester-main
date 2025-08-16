import numpy as np
import os
import matplotlib.pyplot as plt

def perform_deim_analysis_q_loss():
    fom_snapshot_dir = 'fom_snapshots'
    snapshots_path = os.path.join(fom_snapshot_dir, 'fom_q_loss_pipe_snapshots.npy')
    output_prefix = 'deim_q_loss' # Prefix for output files

    if not os.path.exists(snapshots_path):
        print(f"ERROR: Snapshot file not found: {snapshots_path}")
        print("Please ensure you have run Main.py to generate snapshots.")
        return

    try:
        # Load the Q_loss_Pipe snapshots (num_pipes x num_snapshots)
        # Each column is a snapshot, each row corresponds to a pipe
        q_loss_snapshots = np.load(snapshots_path)
        print(f"Successfully loaded {snapshots_path}")
        print(f"Shape of Q_loss_Pipe snapshots: {q_loss_snapshots.shape}")

        if q_loss_snapshots.size == 0:
            print("Q_loss_Pipe snapshots array is empty. Cannot perform DEIM analysis.")
            return
        
        num_elements, num_s = q_loss_snapshots.shape # num_elements is num_pipes
        if num_s < 1:
            print("Not enough snapshots to perform SVD.")
            return

        # Perform Singular Value Decomposition (SVD)
        # Snapshots = U * Sigma * V_transpose
        # U contains the basis vectors for Q_loss_Pipe
        # Sigma contains the singular values
        print("\nPerforming SVD on Q_loss_Pipe snapshots...")
        U, s, Vt = np.linalg.svd(q_loss_snapshots, full_matrices=False)
        print("SVD complete.")

        print(f"Shape of U (basis for Q_loss_Pipe): {U.shape}")
        print(f"Shape of s (singular values for Q_loss_Pipe): {s.shape}")
        print(f"Shape of Vt (V transpose): {Vt.shape}")

        # Singular values indicate the energy of each mode
        print("\nSingular values for Q_loss_Pipe (s_q_loss):")
        print(s)

        # Calculate the energy captured by the modes
        squared_singular_values = s**2
        total_energy = np.sum(squared_singular_values)

        print(f"\nTotal energy (sum of squared s_q_loss values): {total_energy:.4e}")
        
        cumulative_energy = np.cumsum(squared_singular_values) / total_energy

        print("\nCumulative energy captured by retaining k modes for Q_loss_Pipe:")
        for i, energy_fraction in enumerate(cumulative_energy):
            print(f"  {i+1} mode(s): {energy_fraction*100:.4f}%")

        # Save DEIM analysis results for later use
        u_matrix_path = os.path.join(fom_snapshot_dir, f'{output_prefix}_U_matrix.npy')
        singular_values_path = os.path.join(fom_snapshot_dir, f'{output_prefix}_s_values.npy')
        
        np.save(u_matrix_path, U)
        print(f"\nSaved U matrix for Q_loss_Pipe to: {u_matrix_path}")
        np.save(singular_values_path, s)
        print(f"Saved singular values for Q_loss_Pipe to: {singular_values_path}")

        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.semilogy(np.arange(1, len(s) + 1), s, 'o-', markersize=5)
        plt.title(f'Singular Values for Q_loss_Pipe (Log Scale)')
        plt.xlabel('Mode Number')
        plt.ylabel('Singular Value')
        plt.grid(True, which="both", ls="--")
        plt.xticks(np.arange(1, len(s) + 1))

        plt.subplot(1, 2, 2)
        plt.plot(np.arange(1, len(cumulative_energy) + 1), cumulative_energy * 100, 'o-', markersize=5)
        plt.title(f'Cumulative Energy for Q_loss_Pipe')
        plt.xlabel('Number of Modes Retained')
        plt.ylabel('Cumulative Energy (%)')
        plt.grid(True, which="both", ls="--")
        plt.xticks(np.arange(1, len(cumulative_energy) + 1))
        plt.yticks(np.arange(0, 101, 10))
        
        plt.tight_layout()
        plot_path = os.path.join(fom_snapshot_dir, f'{output_prefix}_svalue_analysis.png')
        plt.savefig(plot_path)
        print(f"Saved Q_loss_Pipe singular value plots to: {plot_path}")

    except FileNotFoundError:
        print(f"ERROR: Snapshot file not found at {snapshots_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    perform_deim_analysis_q_loss() 