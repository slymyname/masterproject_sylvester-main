'''
Compares FOM and ROM snapshots and calculates accuracy metrics.
'''
import numpy as np
import os
import datetime

FOM_DIR = 'fom_snapshots'
ROM_DIR = 'rom_snapshots'

# Variables to compare: 
# - key: internal reference name
# - name: display name for reports
# - fom_file_pattern: pattern for FOM snapshot filename 
# - rom_file_pattern: pattern for ROM snapshot filename 
# - dim_type: 'vector_x_time': nodes/pipes vs time or 'scalar_over_time': sum vs time
VARIABLES_TO_COMPARE = {
    'T': {'name': 'Nodal Temperatures', 'file_suffix': '_T_snapshots.npy', 'dim_type': 'vector_x_time', 'element_name': 'Node'},
    'q_loss_pipe': {'name': 'Pipe Heat Losses', 'file_suffix': '_q_loss_pipe_snapshots.npy', 'dim_type': 'vector_x_time', 'element_name': 'Pipe'},
    'heat_losses_sum': {'name': 'Total Heat Loss Sum (per step)', 'file_suffix': '_heat_losses_sum_snapshots.npy', 'dim_type': 'scalar_over_time'},
    'tsupply': {'name': 'Supply Temperature Used', 'file_suffix': '_tsupply_snapshots.npy', 'dim_type': 'scalar_over_time'},
    't_out': {'name': 'Outside Temperature Used', 'file_suffix': '_t_out_snapshots.npy', 'dim_type': 'scalar_over_time'},
    'rho_nodes': {'name': 'Nodal Densities', 'file_suffix': '_rho_nodes_snapshots.npy', 'dim_type': 'vector_x_time', 'element_name': 'Node'},
    'cp_nodes': {'name': 'Nodal Specific Heats', 'file_suffix': '_cp_nodes_snapshots.npy', 'dim_type': 'vector_x_time', 'element_name': 'Node'},
}

DATETIME_FILE_SUFFIX = '_datetimes.npy'

def load_data(directory, engine_prefix, file_suffix):
    '''Loads .npy file from the specified directory.'''
    file_path = os.path.join(directory, f"{engine_prefix}{file_suffix}")
    if not os.path.exists(file_path):
        print(f"  WARNING: File not found: {file_path}")
        return None
    try:
        data = np.load(file_path, allow_pickle=True)
        # Ensure tsupply and t_out are 1D if they were saved as (Ns, 1)
        if file_suffix == '_tsupply_snapshots.npy' or file_suffix == '_t_out_snapshots.npy':
            if data.ndim > 1:
                data = data.flatten()
        return data
    except Exception as e:
        print(f"  ERROR: Could not load file {file_path}: {e}")
        return None

def calculate_metrics(fom_data, rom_data, var_config):
    '''Calculates and prints accuracy metrics.'''
    data_name = var_config['name']
    dim_type = var_config['dim_type']
    element_name = var_config.get('element_name', 'Element')

    print(f"--- Metrics for: {data_name} ---")

    if fom_data.shape != rom_data.shape:
        print(f"  ERROR: Shape mismatch! FOM: {fom_data.shape}, ROM: {rom_data.shape}. Cannot compare.")
        return

    abs_error = np.abs(fom_data - rom_data)
    squared_error = (fom_data - rom_data)**2

    global_mae = np.mean(abs_error)
    global_rmse = np.sqrt(np.mean(squared_error))
    
    # Calculate RRMSE, handle potential division by zero if FOM data is all zeros
    fom_norm_sq = np.mean(fom_data**2)
    if fom_norm_sq < 1e-12: # Threshold to avoid division by near-zero
        global_rrmse = np.inf if global_rmse > 1e-9 else 0.0
        print(f"  WARNING: FOM data norm is near zero for {data_name}. RRMSE might be inf or 0.")
    else:
        global_rrmse = global_rmse / np.sqrt(fom_norm_sq)

    print(f"  Global MAE: {global_mae:.4e}")
    print(f"  Global RMSE: {global_rmse:.4e}")
    print(f"  Global RRMSE: {global_rrmse:.4%}")

    if dim_type == 'vector_x_time':
        # Metrics per element, averaged over time
        mae_per_element = np.mean(abs_error, axis=1)
        rmse_per_element = np.sqrt(np.mean(squared_error, axis=1))
        
        # RRMSE per element
        fom_norm_sq_per_element = np.mean(fom_data**2, axis=1)
        rrmse_per_element = np.full_like(rmse_per_element, np.inf)
        valid_fom_norm_mask = fom_norm_sq_per_element >= 1e-12
        
        rrmse_per_element[valid_fom_norm_mask] = rmse_per_element[valid_fom_norm_mask] / np.sqrt(fom_norm_sq_per_element[valid_fom_norm_mask])
        # For elements where FOM norm is zero, if RMSE is also zero, RRMSE is 0, else Inf.
        rrmse_per_element[~valid_fom_norm_mask & (rmse_per_element < 1e-9)] = 0.0
        
        print(f"  Average MAE across all {element_name}s (time-averaged): {np.mean(mae_per_element):.4e}")
        print(f"  Average RMSE across all {element_name}s (time-averaged): {np.mean(rmse_per_element):.4e}")
        print(f"  Average RRMSE across all {element_name}s (time-averaged): {np.mean(rrmse_per_element[np.isfinite(rrmse_per_element)]):.4%}") # Avg of finite RRMSEs
        
        

    elif dim_type == 'scalar_over_time':
       
        safe_fom_data = np.where(np.abs(fom_data) < 1e-9, 1e-9, fom_data) # Avoid division by zero
        mape = np.mean(np.abs((fom_data - rom_data) / safe_fom_data))
        print(f"  MAPE: {mape:.4%}")
    print("-" * 25)

def main():
    print("Starting FOM vs ROM Snapshot Comparison...")

    fom_datetimes = load_data(FOM_DIR, 'fom', DATETIME_FILE_SUFFIX)
    rom_datetimes = load_data(ROM_DIR, 'rom', DATETIME_FILE_SUFFIX)

    if fom_datetimes is None or rom_datetimes is None:
        print("ERROR: Could not load datetime arrays. Cannot proceed.")
        return
    
    if fom_datetimes.size == 0 or rom_datetimes.size == 0:
        print("ERROR: Datetime arrays are empty. Cannot proceed.")
        return
    
    common_datetimes_list = []
    fom_common_indices = []
    rom_common_indices = []

   
    try:
        set_rom_datetimes = set(rom_datetimes)
        for i, dt_fom in enumerate(fom_datetimes):
            if dt_fom in set_rom_datetimes:
                # Find the corresponding index in rom_datetimes (can be slow for large non-unique lists)
                # This assumes unique datetimes in rom_datetimes for simplicity here.
               
                try:
                    j = np.where(rom_datetimes == dt_fom)[0][0]
                    common_datetimes_list.append(dt_fom)
                    fom_common_indices.append(i)
                    rom_common_indices.append(j)
                except IndexError:
                    pass # Should not happen if dt_fom was in set_rom_datetimes and unique
    except TypeError: # If datetime objects are not directly hashable in a set as expected
        print("WARNING: Datetime objects might not be directly hashable for set optimization. Using slower list iteration for alignment.")
        for i, dt_fom in enumerate(fom_datetimes):
            for j, dt_rom in enumerate(rom_datetimes):
                if dt_fom == dt_rom:
                    common_datetimes_list.append(dt_fom)
                    fom_common_indices.append(i)
                    rom_common_indices.append(j)
                    break # Assuming first match is the one we want (unique datetimes)

    if not common_datetimes_list:
        print("ERROR: No common timesteps found between FOM and ROM snapshots. Cannot compare.")
        return

    print(f"Found {len(common_datetimes_list)} common timesteps for comparison.")
    print(f"Comparison period: {common_datetimes_list[0]} to {common_datetimes_list[-1]}")

    fom_common_indices = np.array(fom_common_indices, dtype=int)
    rom_common_indices = np.array(rom_common_indices, dtype=int)

    for var_key, var_config in VARIABLES_TO_COMPARE.items():
        print(f"\nComparing: {var_config['name']} ({var_key})")
        fom_data_full = load_data(FOM_DIR, 'fom', var_config['file_suffix'])
        rom_data_full = load_data(ROM_DIR, 'rom', var_config['file_suffix'])

        if fom_data_full is None or rom_data_full is None:
            print(f"  Skipping {var_config['name']} due to missing data.")
            continue

        # Slice data to common timesteps
        # For data like (num_elements, num_snapshots), slice along axis 1 (time)
        # For data like (num_snapshots,), slice along axis 0 (time)
        if var_config['dim_type'] == 'vector_x_time':
            if fom_data_full.ndim < 2 or rom_data_full.ndim < 2:
                print(f"  ERROR: Expected 2D data for {var_config['name']}, but got FOM shape {fom_data_full.shape}, ROM shape {rom_data_full.shape}")
                continue
            fom_data_common = fom_data_full[:, fom_common_indices]
            rom_data_common = rom_data_full[:, rom_common_indices]
        elif var_config['dim_type'] == 'scalar_over_time':
            if fom_data_full.ndim != 1 or rom_data_full.ndim != 1:
                 # Handle cases where it might be (Ns, 1) due to np.array([val]) in loop
                if fom_data_full.ndim == 2 and fom_data_full.shape[1] == 1: fom_data_full = fom_data_full.flatten()
                if rom_data_full.ndim == 2 and rom_data_full.shape[1] == 1: rom_data_full = rom_data_full.flatten()
                if fom_data_full.ndim != 1 or rom_data_full.ndim != 1:
                    print(f"  ERROR: Expected 1D data for {var_config['name']}, but got FOM shape {fom_data_full.shape}, ROM shape {rom_data_full.shape} after flatten attempt.")
                    continue
            fom_data_common = fom_data_full[fom_common_indices]
            rom_data_common = rom_data_full[rom_common_indices]
        else:
            print(f"  ERROR: Unknown dim_type '{var_config['dim_type']}' for {var_config['name']}")
            continue
            
        if fom_data_common.size == 0 or rom_data_common.size == 0:
            print(f"  WARNING: Data for {var_config['name']} is empty after slicing to common timesteps.")
            continue

        calculate_metrics(fom_data_common, rom_data_common, var_config)

    print("\nComparison complete.")

if __name__ == '__main__':
    main() 