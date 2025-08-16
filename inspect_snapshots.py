import numpy as np
import os

fom_snapshot_dir = 'fom_snapshots'

t_snapshots_path = os.path.join(fom_snapshot_dir, 'fom_T_snapshots.npy')
print(f"--- Inspecting Temperature Snapshots ({t_snapshots_path}) ---")
if os.path.exists(t_snapshots_path):
    try:
        t_snapshots = np.load(t_snapshots_path)
        print(f"Successfully loaded {t_snapshots_path}")
        print(f"Shape of T_snapshots: {t_snapshots.shape}")
        
        if t_snapshots.size > 0:
            num_nodes, num_snapshots = t_snapshots.shape
            print(f"Number of nodes: {num_nodes}")
            print(f"Number of snapshots: {num_snapshots}")
            
            print("\nFirst T_snapshot (first 5 and last 5 temps):")
            # if num_nodes > 10:
            #     print(f"  First 5: {t_snapshots[:5, 0]}")
            #     print(f"  Last 5:  {t_snapshots[-5:, 0]}")
            # else:
            print(t_snapshots[:, 0])
            
            if num_snapshots > 1:
                print("\nLast T_snapshot (first 5 and last 5 temps):")
                if num_nodes > 10:
                    print(f"  First 5: {t_snapshots[:5, -1]}")
                    print(f"  Last 5:  {t_snapshots[-5:, -1]}")
                else:
                    print(t_snapshots[:, -1])
            
            print(f"\nOverall min temperature: {np.min(t_snapshots):.2f}°C")
            print(f"Overall max temperature: {np.max(t_snapshots):.2f}°C")
            print(f"Mean temperature of first snapshot: {np.mean(t_snapshots[:, 0]):.2f}°C")
            if num_snapshots > 1:
                print(f"Mean temperature of last snapshot: {np.mean(t_snapshots[:, -1]):.2f}°C")
        else:
            print("T_snapshots array is empty.")
            
    except Exception as e:
        print(f"Error loading or processing {t_snapshots_path}: {e}")
else:
    print(f"File not found: {t_snapshots_path}")

print("\n" + "-"*50 + "\n") # Separator

tsupply_snapshots_path = os.path.join(fom_snapshot_dir, 'fom_tsupply_snapshots.npy')
print(f"--- Inspecting Tsupply Snapshots ({tsupply_snapshots_path}) ---")
if os.path.exists(tsupply_snapshots_path):
    try:
        tsupply_snapshots = np.load(tsupply_snapshots_path)
        print(f"Successfully loaded {tsupply_snapshots_path}")
        print(f"Shape of Tsupply_snapshots: {tsupply_snapshots.shape}") 

    
        if tsupply_snapshots.ndim > 1:
            tsupply_snapshots = tsupply_snapshots.flatten()
            print(f"Flattened Tsupply_snapshots to shape: {tsupply_snapshots.shape}")

        if tsupply_snapshots.size > 0:
            num_tsupply_values = tsupply_snapshots.shape[0]
            print(f"Number of Tsupply values recorded: {num_tsupply_values}")
            
            print("\nFirst 5 Tsupply values:")
            print(tsupply_snapshots[:5])
            
            if num_tsupply_values > 5:
                print("\nLast 5 Tsupply values:")
                print(tsupply_snapshots[-5:])
               
            print(f"\nOverall min Tsupply: {np.min(tsupply_snapshots):.2f}°C")
            print(f"Overall max Tsupply: {np.max(tsupply_snapshots):.2f}°C")
            print(f"Mean Tsupply: {np.mean(tsupply_snapshots):.2f}°C")
        else:
            print("Tsupply_snapshots array is empty.")
            
    except Exception as e:
        print(f"Error loading or processing {tsupply_snapshots_path}: {e}")
else:
    print(f"File not found: {tsupply_snapshots_path}") 

print("\n" + "-"*50 + "\n") # Separator

#Inspect Q_loss_Pipe_snapshots
q_loss_pipe_snapshots_path = os.path.join(fom_snapshot_dir, 'fom_q_loss_pipe_snapshots.npy')
print(f"--- Inspecting Q_loss_Pipe Snapshots ({q_loss_pipe_snapshots_path}) ---")
if os.path.exists(q_loss_pipe_snapshots_path):
    try:
        q_loss_pipe_snapshots = np.load(q_loss_pipe_snapshots_path)
        print(f"Successfully loaded {q_loss_pipe_snapshots_path}")
        print(f"Shape of Q_loss_Pipe_snapshots: {q_loss_pipe_snapshots.shape}")
        
        if q_loss_pipe_snapshots.size > 0:
            num_pipes, num_snapshots_qloss = q_loss_pipe_snapshots.shape
            print(f"Number of pipes: {num_pipes}")
            print(f"Number of snapshots: {num_snapshots_qloss}")
            
            print("\nFirst Q_loss_Pipe snapshot (first 5 and last 5 values if available):")
            if num_pipes > 10:
                print(f"  First 5: {q_loss_pipe_snapshots[:5, 0]}")
                print(f"  Last 5:  {q_loss_pipe_snapshots[-5:, 0]}")
            else:
                print(q_loss_pipe_snapshots[:, 0])
            
            if num_snapshots_qloss > 1:
                print("\nLast Q_loss_Pipe snapshot (first 5 and last 5 values if available):")
                if num_pipes > 10:
                    print(f"  First 5: {q_loss_pipe_snapshots[:5, -1]}")
                    print(f"  Last 5:  {q_loss_pipe_snapshots[-5:, -1]}")
                else:
                    print(q_loss_pipe_snapshots[:, -1])
            
            print(f"\nOverall min Q_loss_Pipe value: {np.min(q_loss_pipe_snapshots):.2f} W")
            print(f"Overall max Q_loss_Pipe value: {np.max(q_loss_pipe_snapshots):.2f} W")
            print(f"Mean Q_loss_Pipe of first snapshot: {np.mean(q_loss_pipe_snapshots[:, 0]):.2f} W")
            if num_snapshots_qloss > 1:
                print(f"Mean Q_loss_Pipe of last snapshot: {np.mean(q_loss_pipe_snapshots[:, -1]):.2f} W")
            # Sum of heat losses for the first and last snapshot
            print(f"Sum of Q_loss_Pipe for first snapshot: {np.sum(q_loss_pipe_snapshots[:, 0]):.2f} W")
            if num_snapshots_qloss > 1:
                print(f"Sum of Q_loss_Pipe for last snapshot: {np.sum(q_loss_pipe_snapshots[:, -1]):.2f} W")
        else:
            print("Q_loss_Pipe_snapshots array is empty.")
            
    except Exception as e:
        print(f"Error loading or processing {q_loss_pipe_snapshots_path}: {e}")
else:
    print(f"File not found: {q_loss_pipe_snapshots_path}") 