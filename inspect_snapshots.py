import numpy as np
import os

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