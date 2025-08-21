# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:57:46 2024

@author: sander
"""
import time
import datetime
import os # Ensure os is imported
import cProfile 
import pstats   
import io       

import numpy as np
import matplotlib.pyplot as plt
from StaticData import StaticData

from ProgramModes import ProgramModes
from hydraulicMesh import hydraulicMesh
from TimeSeries import TimeSeries
from ROMtimeseries import ROMTimeSeries # Import ROMTimeSeries
from LoadMeasValues import LoadMeasValues
from optimize import OptimizeStationary

from PlotResults import PlotResult

# Create a profiler object
profiler = cProfile.Profile()

start = time.time()

ProgramModes=ProgramModes()
ProgramModes.SystemMode='TimeSeries'
ProgramModes.MeasurementMode='ideal'
ProgramModes.TempratureCalculationMode='stationary' # 'dynamic' or 'stationary'
ProgramModes.TemperatureDynamicMode='ImplicitUpwind2'   # 'ImplicitUpwind' 'ImplicitUpwind2' and 'Icking' working
ProgramModes.HeatLossMode='Norm' # 'simple' , 'Norm' and 'EMCM'
ProgramModes.HeatLossNextPipe='No' #'Yes' or 'No' 
ProgramModes.randomConsumerMode='default' # 'Gauss'
ProgramModes.PowerInput='None' #for Measurement Values 'None' 'Q' and 'm' 
ProgramModes.TsupplyMode = 'Custom' # 'Default', 'BenchmarkReplay', or 'Custom'

# FOM snapshot generation
ProgramModes.TInput='None' # Optimizer determines Tsupply
ProgramModes.Optimize='No'  # Ensure optimizer is OFF for benchmark replay
ProgramModes.SimulationEngine = 'ROM' # Set to FOM for the first benchmark run

ProgramModes.SnapshotGeneration='No' # This flag is not used by current snapshot logic
ProgramModes.EnableSnapshotSaving = True #snapshots are saved for comparison
ProgramModes.startTime=datetime.datetime(2022, 2, 1, 0, 0, 0, 0) #change startdate
ProgramModes.endTime=datetime.datetime(2022,2, 2, 0, 0,0, 0) 


StaticData=StaticData()
StaticData.readData(ProgramModes)
StaticData.convertAllData()

hydraulicMesh=hydraulicMesh()
hydraulicMesh.MeshPrep(StaticData)

# Conditional instantiation based on SimulationEngine
if ProgramModes.SimulationEngine == 'ROM':
    print("\n[[ INFO ]] Initializing ROMTimeSeries Engine...")
    SimEngine = ROMTimeSeries(StaticData, hydraulicMesh, ProgramModes)
elif ProgramModes.SimulationEngine == 'FOM':
    print("\n[[ INFO ]] Initializing TimeSeries Engine (FOM)...")
    SimEngine = TimeSeries(StaticData, hydraulicMesh, ProgramModes)
else:
    raise ValueError(f"Unknown SimulationEngine: {ProgramModes.SimulationEngine}")

SimEngine.calcThermalResistanceVertical(StaticData)
SimEngine.datestart(ProgramModes)

LoadMeasValues=LoadMeasValues()

# This section is for ROM, so it won't load Tsupply snapshots when Engine is FOM
loaded_fom_tsupply_values = None
fom_tsupply_idx = 0
using_fom_tsupply_snapshots = False

if ProgramModes.SimulationEngine == 'ROM': 
    tsupply_snapshot_path = os.path.join('fom_snapshots', 'fom_tsupply_snapshots.npy')
    try:
        loaded_fom_tsupply_values = np.load(tsupply_snapshot_path)
        if loaded_fom_tsupply_values.ndim > 1:
            loaded_fom_tsupply_values = loaded_fom_tsupply_values.flatten()
        if loaded_fom_tsupply_values.size > 0:
            using_fom_tsupply_snapshots = True
            print(f"[{datetime.datetime.now()}] INFO: Successfully loaded FOM Tsupply snapshots from {tsupply_snapshot_path} for ROM run. Shape: {loaded_fom_tsupply_values.shape}")
        else:
            print(f"[{datetime.datetime.now()}] WARNING: FOM Tsupply snapshot file {tsupply_snapshot_path} is empty. ROM will not use Tsupply snapshots.")
    except FileNotFoundError:
        print(f"[{datetime.datetime.now()}] WARNING: FOM Tsupply snapshot file not found at {tsupply_snapshot_path}. ROM will not use Tsupply snapshots.")
    except Exception as e:
        print(f"[{datetime.datetime.now()}] WARNING: Error loading FOM Tsupply snapshot file {tsupply_snapshot_path}: {e}. ROM will not use Tsupply snapshots.")



#Custom Tsupply Profile Definition

#Define your custom Tsupply values here. This list will be used if ProgramModes.TsupplyMode is 'Custom'.
custom_tsupply_values = [110.50538848, 110.30381051, 111.06077353, 110.59445023, 110.67591816,
                         110.45126529, 109.27746027, 108.51502387, 108.5750523, 108.95780088,
                         108.63269546, 108.00233205, 107.99820189, 107.19264065, 106.40765986,
                         106.56624785, 105.54677866, 105.89819487, 106.7497089, 107.12062118,
                         107.8229505, 107.40800564, 108.67829473, 109.02263402, 107.93369876,
                         107.93410031, 108.00846395, 106.7354604, 108.58462167, 107.882428,
                         107.39443598, 106.55677488, 105.28379949, 106.30298588, 105.50441555,
                         104.52176034, 104.0195072, 104.97679481, 104.72669458, 106.3239832,
                         105.74753909, 106.35919123, 105.90101016, 105.60284368, 105.32943348,
                         105.74192032, 105.05745241, 105.58311348, 105.63719044, 106.30804619,
                         105.58606245, 106.30616779, 105.09004877, 104.27003531, 105.50430354,
                         103.71533864, 104.63805578, 104.35578993, 102.96814555, 103.16421547,
                         103.07314541, 100.96778802, 101.81307898, 102.05518588, 102.59435524,
                         102.12825618, 103.52729421, 103.03389754, 103.78004194, 104.2672075,
                         105.52277801, 106.10188857, 106.88698283, 106.61876888, 105.35341783,
                         104.93315822, 105.69336581, 106.26036153, 106.62847915, 107.63611924,
                         107.29263901, 106.889335, 107.45154014, 105.05471139, 105.41793582,
                         105.57988551, 105.91851627, 105.22897674, 104.37049323, 106.08168836,
                         106.51862453, 106.03608488, 105.78135172, 104.63464744, 104.19230685,
                         103.86042585, 103.52817392]
#Custom Tsupply Profile Definition 

#Benchmark Tsupply Profile Loading 
benchmark_tsupply_values = None
benchmark_idx = 0
using_benchmark_tsupply = False
if ProgramModes.TsupplyMode == 'BenchmarkReplay':
    benchmark_profile_path = os.path.join('fom_snapshots', 'benchmark_tsupply_profile.npy')
    try:
        benchmark_tsupply_values = np.load(benchmark_profile_path)
        if benchmark_tsupply_values.size > 0:
            using_benchmark_tsupply = True
            print(f"[{datetime.datetime.now()}] INFO: Successfully loaded Benchmark Tsupply profile from {benchmark_profile_path}. Shape: {benchmark_tsupply_values.shape}")
        else:
            print(f"[{datetime.datetime.now()}] CRITICAL: Benchmark Tsupply profile file {benchmark_profile_path} is empty. Exiting.")
            exit()
    except FileNotFoundError:
        print(f"[{datetime.datetime.now()}] CRITICAL: Benchmark Tsupply profile file not found at {benchmark_profile_path}. Please run create_benchmark_tsupply.py first. Exiting.")
        exit()
elif ProgramModes.TsupplyMode == 'Custom':
    if custom_tsupply_values:
        benchmark_tsupply_values = np.array(custom_tsupply_values)
        using_benchmark_tsupply = True
        print(f"[{datetime.datetime.now()}] INFO: Using custom Tsupply profile defined in Main.py. Shape: {benchmark_tsupply_values.shape}")
    else:
        print(f"[{datetime.datetime.now()}] WARNING: Custom Tsupply profile is empty. Continuing with default Tsupply behavior.")

# END: Benchmark Tsupply Profile Loading 
elif ProgramModes.TsupplyMode == 'Default':
    Tsupply_current_step_value = np.array([120.0])
    print(f"[{datetime.datetime.now()}] INFO: Using default Tsupply value: {Tsupply_current_step_value}")

if not(ProgramModes.PowerInput=='None' and ProgramModes.TInput=='None'):
    LoadMeasValues.loadmat(SimEngine.datetimeStationary)

# Instantiate Optimizer if ProgramModes.Optimize is 'Yes'
if ProgramModes.Optimize=='Yes': 
    Optimize=OptimizeStationary()
    print(f"[{datetime.datetime.now()}] INFO: Optimizer is ENABLED.")
else:
    Optimize = None # Explicitly set to None if not used
    print(f"[{datetime.datetime.now()}] INFO: Optimizer is DISABLED.")

PlotResult=PlotResult(ProgramModes)

LoadMeasValues.TemperatureOutside(r'prm/produkt_tu_stunde_19970701_20231231_03379.txt')

#Snapshot Collection Setup  
# THIS IS THE CORRECTED LOGIC
if ProgramModes.SimulationEngine == 'ROM':
    snapshot_dir = 'rom_snapshots'
    if ProgramModes.EnableSnapshotSaving:
        print(f"[[ INFO ]] ROM snapshots will be saved to: {snapshot_dir}")
    else:
        print(f"[[ INFO ]] ROM run: Snapshot saving is DISABLED for this run.")
elif ProgramModes.SimulationEngine == 'FOM':
    # Special case for benchmark runs if needed, otherwise default to fom_snapshots
    if ProgramModes.TsupplyMode == 'BenchmarkReplay':
        snapshot_dir = 'benchmark_fom_output'
        if ProgramModes.EnableSnapshotSaving:
            print(f"[[ INFO ]] BENCHMARK MODE: FOM snapshots will be saved to: {snapshot_dir}")
    else:
        snapshot_dir = 'fom_snapshots' # Default FOM snapshots go here
        if ProgramModes.EnableSnapshotSaving:
            print(f"[[ INFO ]] FOM snapshots will be saved to: {snapshot_dir} (will overwrite existing if any)")
        else:
            print(f"[[ INFO ]] FOM run: Snapshot saving is DISABLED. Existing fom_snapshots will NOT be overwritten.")
else:
    snapshot_dir = 'unknown_engine_snapshots' 
    print(f"[[ WARNING ]] Unknown simulation engine for snapshot directory: {snapshot_dir}")

if not os.path.exists(snapshot_dir):
    os.makedirs(snapshot_dir)
    print(f"Created directory: {snapshot_dir}")

current_T_snapshots = []
current_rho_nodes_snapshots = []
current_cp_nodes_snapshots = []
current_viscosity_nodes_snapshots = []
current_heat_losses_snapshots = [] 
current_q_loss_pipe_snapshots = [] 
current_t_out_snapshots = []
current_tsupply_snapshots = [] # This will capture the optimized Tsupply
current_datetime_snapshots = []
# End Snapshot Collection Setup 

# Initial Tsupply for the very first step (or if optimizer is off)
# The optimizer will refine this value at each step if active.
Tsupply_current_step_value = np.array([120.0]) 
print(f"[{datetime.datetime.now()}] Initial Tsupply guess for first step: {Tsupply_current_step_value}")
Tsupply_default_fallback=120.0
if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
    snapshot_status_message = "ENABLED" if ProgramModes.EnableSnapshotSaving else "DISABLED"
    print(f"[{datetime.datetime.now()}] Starting stationary simulation loop for {ProgramModes.SimulationEngine} (snapshot collection {snapshot_status_message})...")
    
    profiler.enable() 

    while SimEngine.datetimeStationary[-1]<ProgramModes.endTime:
        print(f"\n[{datetime.datetime.now()}] --- Main Timestep Start: {SimEngine.datetimeStationary[-1]} ({ProgramModes.SimulationEngine}) ---")
        
        LoadMeasValues.updateTout(SimEngine, ProgramModes)
        print(f"[{datetime.datetime.now()}] Updated T_out: {SimEngine.T_out:.2f}°C")
        
        # Determine Tsupply for the current step
        if using_benchmark_tsupply:
            if benchmark_idx < len(benchmark_tsupply_values):
                Tsupply_current_step_value = np.array([benchmark_tsupply_values[benchmark_idx]])
                benchmark_idx += 1
            else:
                # FOM-compatible behavior: default to 120.0 when Tsupply values are exhausted
                Tsupply_current_step_value = np.array([120.0])
                print(f"[{datetime.datetime.now()}] BENCHMARK Tsupply profile exhausted, defaulting to 120.0°C")
            print(f"[{datetime.datetime.now()}] BENCHMARK using Tsupply from profile: {Tsupply_current_step_value}")

        elif ProgramModes.SimulationEngine == 'FOM':
            if ProgramModes.Optimize == 'Yes' and Optimize is not None:
                # The Tsupply_current_step_value from the PREVIOUS step (or initial guess) is used as the starting point for the optimizer
                print(f"[{datetime.datetime.now()}] FOM Optimizer input Tsupply: {Tsupply_current_step_value}")
                Tsupply_optimized = Optimize.optimize(Tsupply_current_step_value, StaticData, SimEngine, hydraulicMesh, ProgramModes, LoadMeasValues)
                Tsupply_current_step_value = Tsupply_optimized # Update for this step and for next step's initial guess
                print(f"[{datetime.datetime.now()}] FOM Optimizer output Tsupply: {Tsupply_current_step_value}")
            elif ProgramModes.TInput == 'Yes': # This case should not be active due to TInput='None'
                Tsupply_from_meas = np.array([LoadMeasValues.SupplyTemperature(SimEngine.datetimeStationary[-1])])
                if Tsupply_from_meas.ndim > 1: Tsupply_from_meas = Tsupply_from_meas.flatten()
                Tsupply_current_step_value = Tsupply_from_meas
                print(f"[{datetime.datetime.now()}] FOM Tsupply from measurements: {Tsupply_current_step_value}")
            # else: Tsupply_current_step_value remains the default/previous optimized value
        
        elif ProgramModes.SimulationEngine == 'ROM':
            # ROM Tsupply logic (using loaded_fom_tsupply_values)
            if ProgramModes.Optimize == 'Yes' and Optimize is not None:
                # The Tsupply_current_step_value from the PREVIOUS step (or initial guess) is used as the starting point for the optimizer
                print(f"[{datetime.datetime.now()}] FOM Optimizer input Tsupply: {Tsupply_current_step_value}")
                Tsupply_optimized = Optimize.optimize(Tsupply_current_step_value, StaticData, SimEngine, hydraulicMesh, ProgramModes, LoadMeasValues)
                Tsupply_current_step_value = Tsupply_optimized # Update for this step and for next step's initial guess
                print(f"[{datetime.datetime.now()}] FOM Optimizer output Tsupply: {Tsupply_current_step_value}")
            
            elif using_fom_tsupply_snapshots:
                if fom_tsupply_idx < len(loaded_fom_tsupply_values):
                    Tsupply_current_step_value = np.array([loaded_fom_tsupply_values[fom_tsupply_idx]])
                    if Tsupply_current_step_value.ndim > 1: Tsupply_current_step_value = Tsupply_current_step_value.flatten()
                    fom_tsupply_idx += 1
                else:
                    # FOM-compatible behavior: default to 120.0 when FOM Tsupply snapshots are exhausted
                    Tsupply_current_step_value = np.array([120.0])
                    if Tsupply_current_step_value.ndim > 1: Tsupply_current_step_value = Tsupply_current_step_value.flatten()
                    print(f"[{datetime.datetime.now()}] ROM FOM Tsupply snapshots exhausted, defaulting to 120.0°C")
                print(f"[{datetime.datetime.now()}] ROM using Tsupply from snapshot: {Tsupply_current_step_value}")
            else:
                Tsupply_current_step_value = Tsupply_default_fallback # Defined earlier
                print(f"[{datetime.datetime.now()}] ROM Tsupply: Not using snapshots, using default: {Tsupply_current_step_value}")

        print(f"[{datetime.datetime.now()}] Using Tsupply for timestepStationary: {Tsupply_current_step_value}")
        SimEngine.timestepStationary(StaticData,hydraulicMesh,ProgramModes,LoadMeasValues, Tsupply_current_step_value,True)
        
        PlotResult.T.append(SimEngine.T)
        PlotResult.m_i.append(SimEngine.m_i)
        PlotResult.p.append(SimEngine.p)
        PlotResult.m_consumer.append(SimEngine.m_consumer)
        PlotResult.HL1.append(SimEngine.HeatLoses) 
        PlotResult.PumpPower.append(SimEngine.hydraulicPower)
        
        # --- MPC State Collection ---
        # Store the Tsupply that was used for this step
        PlotResult.T_supply.append(Tsupply_current_step_value)

        # Store plant pressures (Outlet from Supplies, Inlet from Pumps)
        # Assuming the first supply is the main plant outlet and first pump is main plant inlet
        plant_outlet_pressure = SimEngine.p[StaticData.vorlauf_SupplieNodes[0]]
        plant_inlet_pressure = SimEngine.p[StaticData.Pumps_startnode[0]]
        PlotResult.p_plant.append([plant_outlet_pressure, plant_inlet_pressure])

        # Calculate and store the minimum pressure difference across all consumers
        consumer_p_diffs = SimEngine.p[StaticData.vorlauf_consumptionnode] - SimEngine.p[StaticData.ruecklauf_consumptionnode]
        min_consumer_p_diff = np.min(consumer_p_diffs)
        PlotResult.p_diff_consumer_min.append(min_consumer_p_diff)
        
        # --- Record Data for NN Training ---
        current_ambient_temp = SimEngine.T_out # Corrected: Get temp from SimEngine
        # Calculate total heat delivered to consumers (a proxy for demand)
        cp_water = 4186 # J/kgK, approximate
        heat_delivered = np.sum(SimEngine.m_consumer * cp_water * 
                                (SimEngine.T[StaticData.vorlauf_consumptionnode] - SimEngine.T[StaticData.ruecklauf_consumptionnode]))
        
        snapshot = {
            # Inputs
            'T_supply': Tsupply_current_step_value,
            'T_ambient': current_ambient_temp,
            'Heat_Demand': heat_delivered,
            # Outputs
            'Heat_Loss': SimEngine.HeatLoses,
            'Pump_Power': SimEngine.hydraulicPower,
            'P_outlet': plant_outlet_pressure,
            'P_inlet': plant_inlet_pressure,
            'P_diff_min': min_consumer_p_diff
        }
        mpc_data_recorder.append(snapshot)

        
        if ProgramModes.EnableSnapshotSaving: 
            current_datetime_snapshots.append(SimEngine.datetimeStationary[-1])
            current_T_snapshots.append(SimEngine.T.copy())
            current_rho_nodes_snapshots.append(SimEngine.rho_Nodes.copy())
            current_cp_nodes_snapshots.append(SimEngine.cp_Nodes.copy())
            current_viscosity_nodes_snapshots.append(SimEngine.viscosity_Nodes.copy())
            current_heat_losses_snapshots.append(np.sum(SimEngine.Q_loss_Pipe)) 
            current_q_loss_pipe_snapshots.append(SimEngine.Q_loss_Pipe.copy()) 
            current_t_out_snapshots.append(SimEngine.T_out)
            current_tsupply_snapshots.append(Tsupply_current_step_value.copy()) # Save the Tsupply used for this step
        
        print(f"[{datetime.datetime.now()}] --- Main Timestep End: {SimEngine.datetimeStationary[-1]} ({ProgramModes.SimulationEngine}) ---")

    profiler.disable() 

elif ProgramModes.TempratureCalculationMode=='dynamic':
    print(f"[{datetime.datetime.now()}] Starting dynamic simulation setup ({ProgramModes.SimulationEngine})...")
    LoadMeasValues.updateTout(SimEngine, ProgramModes)
    if ProgramModes.TInput=='Yes':
        Tsupply=np.array([LoadMeasValues.SupplyTemperature(SimEngine.datetimeStationary[-1])]) 
    SimEngine.timestepStationary(StaticData,hydraulicMesh,ProgramModes,LoadMeasValues,Tsupply,False)
    PlotResult.m_i.append(SimEngine.m_i.copy())
    PlotResult.p.append(SimEngine.p.copy())
    PlotResult.T_s.append(SimEngine.T.copy())
    PlotResult.HL1.append(SimEngine.HeatLoses)
    PlotResult.PumpPower.append(SimEngine.hydraulicPower)
    while SimEngine.datetimeStationary[-1]<ProgramModes.endTime:


        while SimEngine.datetimeDynamic[-1]<SimEngine.datetimeStationary[-1]+SimEngine.time_change_stationary:
            if ProgramModes.TInput=='Yes':
                Tsupply=np.array([LoadMeasValues.SupplyTemperature(SimEngine.datetimeDynamic[-1])]) 
            SimEngine.timestepDynamic(StaticData,ProgramModes,Tsupply)
            PlotResult.T.append(SimEngine.T.copy())
            
        LoadMeasValues.updateTout(SimEngine, ProgramModes)
        if not(ProgramModes.PowerInput=='None' and ProgramModes.TInput=='None'):
            LoadMeasValues.loadmat(SimEngine.datetimeStationary)           
        if len(SimEngine.datetimeStationary)>0:
            if ProgramModes.TInput=='Yes':
                Tsupply=np.array([LoadMeasValues.SupplyTemperature(SimEngine.datetimeDynamic[-1])]) 
            SimEngine.timestepDynamicHydraulic(StaticData, hydraulicMesh,ProgramModes,LoadMeasValues,Tsupply)
            PlotResult.m_i.append(SimEngine.m_i.copy())
            PlotResult.p.append(SimEngine.p.copy())
            PlotResult.m_consumer.append(SimEngine.m_consumer)
            PlotResult.T_s.append(SimEngine.T.copy())
            PlotResult.HL1.append(SimEngine.HeatLoses)
            PlotResult.PumpPower.append(SimEngine.hydraulicPower)
            
        
def export_to_csv(plot_result, num_nodes, output_dir='.', filename_prefix='simulation_results'):
    """Exports simulation data to a CSV file.

    Args:
        plot_result (PlotResult): The PlotResult object containing the simulation data.
        num_nodes (int): The number of nodes in the simulation.
        output_dir (str): The directory to save the CSV file in.
        filename_prefix (str): The prefix for the output CSV filename.
    """
    import pandas as pd
    import os

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Prepare data for each time step
    all_data = []
    for i in range(len(plot_result.T)):
        for j in range(num_nodes):
            all_data.append({
                'timestep': i,
                'node_number': j,
                'temperature': plot_result.T[i][j] if i < len(plot_result.T) and j < len(plot_result.T[i]) else 'N/A',
                'pressure': plot_result.p[i][j] if i < len(plot_result.p) and j < len(plot_result.p[i]) else 'N/A',
                'mass_flow': plot_result.m_i[i][j] if i < len(plot_result.m_i) and j < len(plot_result.m_i[i]) else 'N/A'
            })

    # Create a DataFrame and save to CSV
    df = pd.DataFrame(all_data)
    filename = f"{filename_prefix}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    filepath = os.path.join(output_dir, filename)
    df.to_csv(filepath, index=False)
    print(f"Data exported to {filepath}")
    
    # Create a separate summary CSV with network-wide values
    summary_data = []
    for i in range(len(plot_result.T)):
        summary_data.append({
            'timestep': i,
            'total_heat_loss_network': plot_result.HL1[i] if i < len(plot_result.HL1) else 'N/A',
            'total_pump_power_network': plot_result.PumpPower[i] if i < len(plot_result.PumpPower) else 'N/A'
        })
    
    # Save summary CSV
    summary_df = pd.DataFrame(summary_data)
    summary_filename = f"{filename_prefix}_network_summary_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    summary_filepath = os.path.join(output_dir, summary_filename)
    summary_df.to_csv(summary_filepath, index=False)
    print(f"Network summary exported to {summary_filepath}")

# record end time
end = time.time()
 
# print the difference between start 

print("The time of execution of above program is :",
      (end-start), "s")

# Snapshot Saving (Adjusted for dynamic snapshot_dir) 
if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
    if current_datetime_snapshots: # Check if snapshots were collected
        if ProgramModes.EnableSnapshotSaving:
            print(f"\n[{datetime.datetime.now()}] Saving {ProgramModes.SimulationEngine} snapshots to {snapshot_dir}...")
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_datetimes.npy'), np.array(current_datetime_snapshots))
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_T_snapshots.npy'), np.array(current_T_snapshots).T)
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_rho_nodes_snapshots.npy'), np.array(current_rho_nodes_snapshots).T)
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_cp_nodes_snapshots.npy'), np.array(current_cp_nodes_snapshots).T)
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_viscosity_nodes_snapshots.npy'), np.array(current_viscosity_nodes_snapshots).T)
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_heat_losses_sum_snapshots.npy'), np.array(current_heat_losses_snapshots))
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_q_loss_pipe_snapshots.npy'), np.array(current_q_loss_pipe_snapshots).T)
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_t_out_snapshots.npy'), np.array(current_t_out_snapshots))
            np.save(os.path.join(snapshot_dir, f'{ProgramModes.SimulationEngine.lower()}_tsupply_snapshots.npy'), np.array(current_tsupply_snapshots).T if current_tsupply_snapshots and isinstance(current_tsupply_snapshots[0], np.ndarray) else np.array(current_tsupply_snapshots))
            print(f"[{datetime.datetime.now()}] {ProgramModes.SimulationEngine} snapshots saved to '{snapshot_dir}/' directory.")
        else:
            print(f"\n[{datetime.datetime.now()}] Snapshot saving is DISABLED for {ProgramModes.SimulationEngine} run (EnableSnapshotSaving=False).")
    else:
        print(f"[{datetime.datetime.now()}] No {ProgramModes.SimulationEngine} snapshots collected to save.")
#  End Snapshot Saving 

# Export data to CSV
export_to_csv(PlotResult, StaticData.nbNodes, output_dir='simulation_output')
    
PlotResult.plot_combined(SimEngine,ProgramModes) 
PlotResult.plot_mpc_states(SimEngine, StaticData)
# print(sum(PlotResult.HL1)*ProgramModes.stationaryTimeStep/60*1e-6,'MWh') # HL1 is already sums
total_heat_loss_joules = sum(h_sum * ProgramModes.stationaryTimeStep * 60 for h_sum in PlotResult.HL1 if isinstance(h_sum, (int, float, np.number)))
print(f"Total Heat Loss from PlotResult.HL1: {total_heat_loss_joules * 1e-6 / 3600:.4f} MWh ({ProgramModes.SimulationEngine} run)") # Convert Joules to MWh



print("\n\n--- cProfile Stats ---") 
s = io.StringIO()                  
sortby = pstats.SortKey.CUMULATIVE 
ps = pstats.Stats(profiler, stream=s).sort_stats(sortby) 
ps.print_stats(50) 
print(s.getvalue())                
# profiler.dump_stats(os.path.join(snapshot_dir, 'rom_profile.prof')) 
# print(f"Profiling data saved to {os.path.join(snapshot_dir, 'rom_profile.prof')}") 

if __name__ == "__main__":
    
    # --- Initialize Data Structures ---
    PlotResult = PlotResult(ProgramModes)
    mpc_data_recorder = [] # To store data for NN training
    
    if ProgramModes.SystemMode=='TimeSeries' and ProgramModes.TempratureCalculationMode=='stationary':
        
        # --- Save Recorded Data for NN Training ---
        if mpc_data_recorder:
            # Convert list of dicts to a structured NumPy array for easier handling later
            # We define the data type for each column
            dtype = [('T_supply', 'f8'), ('T_ambient', 'f8'), ('Heat_Demand', 'f8'),
                     ('Heat_Loss', 'f8'), ('Pump_Power', 'f8'), ('P_outlet', 'f8'), 
                     ('P_inlet', 'f8'), ('P_diff_min', 'f8')]
            
            # Create an empty array with the correct shape and dtype
            data_array = np.empty(len(mpc_data_recorder), dtype=dtype)
            
            # Fill the array
            for i, record in enumerate(mpc_data_recorder):
                data_array[i] = (record['T_supply'], record['T_ambient'], record['Heat_Demand'],
                                 record['Heat_Loss'], record['Pump_Power'], record['P_outlet'],
                                 record['P_inlet'], record['P_diff_min'])
            
            output_filename = 'mpc_training_data.npy'
            np.save(output_filename, data_array)
            print(f"\nSuccessfully saved MPC training data to '{output_filename}'")
            print(f"Recorded {len(data_array)} timesteps.")

    
    elif ProgramModes.SystemMode=='coupling' and ProgramModes.TempratureCalculationMode=='stationary':
        
        # --- Save Recorded Data for NN Training ---
        if mpc_data_recorder:
            # Convert list of dicts to a structured NumPy array for easier handling later
            # We define the data type for each column
            dtype = [('T_supply', 'f8'), ('T_ambient', 'f8'), ('Heat_Demand', 'f8'),
                     ('Heat_Loss', 'f8'), ('Pump_Power', 'f8'), ('P_outlet', 'f8'), 
                     ('P_inlet', 'f8'), ('P_diff_min', 'f8')]
            
            # Create an empty array with the correct shape and dtype
            data_array = np.empty(len(mpc_data_recorder), dtype=dtype)
            
            # Fill the array
            for i, record in enumerate(mpc_data_recorder):
                data_array[i] = (record['T_supply'], record['T_ambient'], record['Heat_Demand'],
                                 record['Heat_Loss'], record['Pump_Power'], record['P_outlet'],
                                 record['P_inlet'], record['P_diff_min'])
            
            output_filename = 'mpc_training_data.npy'
            np.save(output_filename, data_array)
            print(f"\nSuccessfully saved MPC training data to '{output_filename}'")
            print(f"Recorded {len(data_array)} timesteps.") 
