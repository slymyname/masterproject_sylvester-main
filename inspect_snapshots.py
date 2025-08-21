import numpy as np
import os

<<<<<<< HEAD
def inspect_snapshot(file_path):
    """
    Loads a .npy snapshot file and prints detailed information about its contents.

    Args:
        file_path (str): The full path to the .npy file to inspect.
    """
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}\n")
        return

    try:
        data = np.load(file_path, allow_pickle=True)
        print(f"--- Inspection Report for: {os.path.basename(file_path)} ---")
        
        # Check if the array is empty
        if data.size == 0:
            print("  Status: Array is empty.")
            print("-" * (28 + len(os.path.basename(file_path))) + "\n")
            return

        # Print basic metadata
        print(f"  Shape: {data.shape}")
        print(f"  Data Type: {data.dtype}")
        
        # For numeric data, calculate and print statistics
        if np.issubdtype(data.dtype, np.number):
            print("\n  Statistics:")
            print(f"    Mean:   {np.mean(data):.4f}")
            print(f"    Min:    {np.min(data):.4f}")
            print(f"    Max:    {np.max(data):.4f}")
            print(f"    Std Dev:{np.std(data):.4f}")
        
        # Print a sample of the data
        print("\n  Data Sample:")
        if data.ndim == 0: # Scalar
            print(f"    {data}")
        elif data.ndim == 1: # 1D Array
            print(f"    First 5 elements: {data[:5]}")
        elif data.ndim == 2: # 2D Array
            print("    Top-left 5x5 corner:")
            print(data[:5, :5])
        else: # Higher dimension
            print("    (Data has more than 2 dimensions, showing slice from the start)")
            print(data.flatten()[:10]) # Show first 10 flattened elements as a sample

        print("-" * (28 + len(os.path.basename(file_path))) + "\n")

    except Exception as e:
        print(f"--- Error inspecting file: {os.path.basename(file_path)} ---")
        print(f"  An error occurred: {e}\n")

def main():
    """
    Main function to find and inspect all ROM snapshots in the target directory.
    """
    rom_snapshot_dir = 'rom_snapshots'
    print(f"--- Starting Inspection of ROM Snapshots in '{rom_snapshot_dir}/' ---\n")

    if not os.path.isdir(rom_snapshot_dir):
        print(f"ERROR: Directory not found: '{rom_snapshot_dir}'")
        print("Please ensure you are running this script from the project's root directory.")
        return

    snapshot_files = sorted([f for f in os.listdir(rom_snapshot_dir) if f.endswith('.npy')])

    if not snapshot_files:
        print("No .npy snapshot files found in the directory.")
        return

    print(f"Found {len(snapshot_files)} snapshot files to inspect.\n")

    for filename in snapshot_files:
        full_path = os.path.join(rom_snapshot_dir, filename)
        inspect_snapshot(full_path)

    print("--- Inspection Complete ---")

if __name__ == '__main__':
    main() 
=======
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
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b
