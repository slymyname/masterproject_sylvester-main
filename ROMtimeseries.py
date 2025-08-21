# -*- coding: utf-8 -*-
import numpy as np
import math
from pyXSteam.XSteam import XSteam
import datetime
from StandardConsumptionProfile import StandardConsumptionProfile
from DynamicTemperature import DynamicTemperature
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve,inv, lsqr, splu
from scipy.integrate import quad
import random
import copy
import os # Added for path joining
from rom_logic_optimized import _average_properties_numba, _calc_lambda_numba, _build_and_project_thermic_system_numba



class ROMTimeSeries:
    """Reduced Order Model (ROM) for the network.

    This class runs a fast surrogate of the full thermo‑hydraulic model.
    - Temperatures are reduced with Proper Orthogonal Decomposition (POD)
      using the basis matrix `V_r`.
    - Nonlinear material properties (rho, cp, viscosity) and heat‑loss terms
      are approximated via the Discrete Empirical Interpolation Method (DEIM).
    - Hydraulics (flows/pressures) are purposely kept full‑order to preserve
      robustness; this is a hybrid ROM that trades some speed for reliability.

    The intent is to retain the important physics while cutting the majority
    of the computational cost. Accuracy mainly depends on how representative
    the POD/DEIM snapshots are and on the chosen reduced dimensions.
    """
    
    def __init__(self,StaticData,hydraulicMesh,ProgramModes):
        """Set up ROM state, load POD/DEIM assets and precompute operators.

        Loads the precomputed POD basis and DEIM objects from `fom_snapshots/`,
        prepares constant matrices, and builds small helper structures that are
        reused each timestep (e.g., factorized linear systems for hydraulics).
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        # ROM assets and sizes (chosen offline based on snapshot analyses)
        self.r_dim = 2
        self.m_rho = 2
<<<<<<< HEAD
        self.m_cp = 1
        self.m_q_loss = 5
=======
        self.m_cp = 2
        self.m_q_loss = 2
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b

        fom_snapshot_dir = 'fom_snapshots' # Base directory for ROM data

        # Load POD Basis V_r
        vr_basis_path = os.path.join(fom_snapshot_dir, f'pod_Vr_basis_r{self.r_dim}.npy')
        try:
            self.V_r = np.load(vr_basis_path)
            print(f"Successfully loaded POD basis V_r ({self.V_r.shape}) from {vr_basis_path}")
        except FileNotFoundError:
            print(f"ERROR: POD basis file not found at {vr_basis_path}. Exiting.")
            raise

        # Load DEIM components for rho_Nodes
        rho_indices_path = os.path.join(fom_snapshot_dir, f'deim_rho_indices_m{self.m_rho}.npy')
        rho_proj_mat_path = os.path.join(fom_snapshot_dir, f'deim_rho_proj_matrix_m{self.m_rho}.npy')
        try:
            self.deim_rho_indices = np.load(rho_indices_path)
            self.deim_rho_proj_matrix = np.load(rho_proj_mat_path)
            print(f"Loaded DEIM rho indices ({self.deim_rho_indices.shape}) and proj matrix ({self.deim_rho_proj_matrix.shape})")
        except FileNotFoundError:
            print(f"ERROR: DEIM rho files not found. Searched: {rho_indices_path}, {rho_proj_mat_path}. Exiting.")
            raise

        # Load DEIM components for cp_Nodes
        cp_indices_path = os.path.join(fom_snapshot_dir, f'deim_cp_indices_m{self.m_cp}.npy')
        cp_proj_mat_path = os.path.join(fom_snapshot_dir, f'deim_cp_proj_matrix_m{self.m_cp}.npy')
        try:
            self.deim_cp_indices = np.load(cp_indices_path)
            self.deim_cp_proj_matrix = np.load(cp_proj_mat_path)
            print(f"Loaded DEIM cp indices ({self.deim_cp_indices.shape}) and proj matrix ({self.deim_cp_proj_matrix.shape})")
        except FileNotFoundError:
            print(f"ERROR: DEIM cp files not found. Searched: {cp_indices_path}, {cp_proj_mat_path}. Exiting.")
            raise

        # Load DEIM components for Q_loss_Pipe
        q_loss_indices_path = os.path.join(fom_snapshot_dir, f'deim_q_loss_indices_m{self.m_q_loss}.npy')
        q_loss_proj_mat_path = os.path.join(fom_snapshot_dir, f'deim_q_loss_proj_matrix_m{self.m_q_loss}.npy')
        try:
            self.deim_q_loss_indices = np.load(q_loss_indices_path)
            self.deim_q_loss_proj_matrix = np.load(q_loss_proj_mat_path)
            print(f"Loaded DEIM Q_loss_Pipe indices ({self.deim_q_loss_indices.shape}) and proj matrix ({self.deim_q_loss_proj_matrix.shape})")
        except FileNotFoundError:
            print(f"ERROR: DEIM Q_loss_Pipe files not found. Searched: {q_loss_indices_path}, {q_loss_proj_mat_path}. Exiting.")
            raise

        # Load DEIM components for viscosity_Nodes
<<<<<<< HEAD
        self.m_viscosity = 5
=======
        self.m_viscosity = 2
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b
        viscosity_indices_path = os.path.join(fom_snapshot_dir, f'deim_viscosity_indices_m{self.m_viscosity}.npy')
        viscosity_proj_mat_path = os.path.join(fom_snapshot_dir, f'deim_viscosity_proj_matrix_m{self.m_viscosity}.npy')
        try:
            self.deim_viscosity_indices = np.load(viscosity_indices_path)
            self.deim_viscosity_proj_matrix = np.load(viscosity_proj_mat_path)
            print(f"Loaded DEIM viscosity indices ({self.deim_viscosity_indices.shape}) and proj matrix ({self.deim_viscosity_proj_matrix.shape})")
        except FileNotFoundError:
            print(f"ERROR: DEIM viscosity files not found. Searched: {viscosity_indices_path}, {viscosity_proj_mat_path}. Exiting.")
            raise
    

        # Pre‑factorize the static linear system used in calcm0 (tree flows).
        # Construct A_0_LE sparsely
        num_left_edges = len(hydraulicMesh.leftEdges)
        num_nodes_A_forest = hydraulicMesh.A_Forest.shape[1] # Should be StaticData.nbNodes
        
        rows_A0LE, cols_A0LE, data_A0LE = [], [], []
        for n in range(num_left_edges):
            rows_A0LE.append(n)
            cols_A0LE.append(hydraulicMesh.leftEdges[n])
            data_A0LE.append(1)
        A_0_LE_sparse = csc_matrix((data_A0LE, (rows_A0LE, cols_A0LE)), shape=(num_left_edges, num_nodes_A_forest))

        # Construct A_0_Pipes sparsely
        num_closed_pipes = len(StaticData.closedPipes)
        rows_A0Pipes, cols_A0Pipes, data_A0Pipes = [], [], []
        for n in range(num_closed_pipes):
            rows_A0Pipes.append(n)
            cols_A0Pipes.append(StaticData.closedPipes[n])
            data_A0Pipes.append(1)
        A_0_Pipes_sparse = csc_matrix((data_A0Pipes, (rows_A0Pipes, cols_A0Pipes)), shape=(num_closed_pipes, num_nodes_A_forest))

        # Ensure A_Forest is sparse (it should be from hydraulicMesh)
        if not isinstance(hydraulicMesh.A_Forest, (csc_matrix, csr_matrix)):
            A_Forest_sparse = csc_matrix(hydraulicMesh.A_Forest)
            print("Warning: hydraulicMesh.A_Forest was not sparse, converted to csc_matrix.")
        else:
            A_Forest_sparse = csc_matrix(hydraulicMesh.A_Forest) # Ensure CSC for vstack consistency if needed

        A_sparse_full = csr_matrix(np.vstack([
            A_Forest_sparse.toarray() if isinstance(A_Forest_sparse, csc_matrix) else A_Forest_sparse.toarray(), # Convert to dense before vstack if different sparse formats
            A_0_LE_sparse.toarray() if isinstance(A_0_LE_sparse, csc_matrix) else A_0_LE_sparse.toarray(),
            A_0_Pipes_sparse.toarray() if isinstance(A_0_Pipes_sparse, csc_matrix) else A_0_Pipes_sparse.toarray()
        ]))
        A_sparse_full = csc_matrix(A_sparse_full) # Convert final to CSC

        # Row deletions for A
        # The indices to delete are static
        # Important: original A had shape (rows_A_forest + num_left_edges + num_closed_pipes, num_nodes_A_forest)
        # The idx for deletion are node indices, but they are applied to rows of A.
        # This means specific equations corresponding to these nodes are removed.
        # A_Forest is (nbNodes, nbNodes), so the first nbNodes rows of A_sparse_full correspond to A_Forest rows.
        
        rows_to_delete_A = np.array([StaticData.vorlauf_SupplieNodes[0], StaticData.ruecklauf_SupplieNodes[0]], dtype=int)
        
        # Create a mask for rows to keep
        num_total_rows_A_full = A_sparse_full.shape[0]
        rows_to_keep_A_mask = np.ones(num_total_rows_A_full, dtype=bool)
        rows_to_keep_A_mask[rows_to_delete_A] = False
        
        A_sparse_deleted_nodes = A_sparse_full[rows_to_keep_A_mask, :]

        # Removing all-zero rows. 
        A_dense_for_nonzero_check = A_sparse_deleted_nodes.toarray()
        non_zero_row_mask_A = ~np.all(A_dense_for_nonzero_check == 0, axis=1)
        
        A_final_sparse = A_sparse_deleted_nodes[non_zero_row_mask_A, :]
        
        if A_final_sparse.shape[0] != A_final_sparse.shape[1]:
            print(f"ERROR in ROMTimeSeries.__init__: Final matrix A for calcm0 is not square after operations! Shape: {A_final_sparse.shape}")
            print("This will cause splu to fail. Check matrix construction and deletion logic.")


        print(f"Shape of final A_sparse for calcm0 factorization: {A_final_sparse.shape}")
        try:
            self._calcm0_A_factorization = splu(A_final_sparse)
            print("Successfully pre-factorized matrix A for calcm0.")
        except Exception as e:
            print(f"ERROR during splu factorization in ROMTimeSeries.__init__: {e}")
            print(f"Final A_sparse for calcm0 shape was: {A_final_sparse.shape}. Is it square and non-singular?")

            self._calcm0_A_factorization = None # Indicate failure
            raise # Re-raise the exception to halt if critical

        # Store the masks used for L, so L can be filtered consistently in calcm0
        self._calcm0_L_rows_to_delete_indices = rows_to_delete_A # These are the original indices applied to L
        self._calcm0_L_non_zero_row_mask = non_zero_row_mask_A # This mask was derived from A's state *after* node deletion.
                                                              # It needs to be applied to L *after* L has also had nodes deleted.

        # Pre-computation for calcp 
        # A_T is static and can be pre-computed
        A_T_full = hydraulicMesh.A_Forest.transpose() * (-1)
        A_T_full = A_T_full[:StaticData.nbPipes,:]
        self._calcp_A_T_sparse = csc_matrix(A_T_full[~np.all(A_T_full == 0, axis=1)])

        # FP and M_S parts of P matrix are also static
        self._calcp_FP_sparse = csc_matrix(np.zeros((StaticData.pFixed.shape[0], StaticData.nbNodes)))
        for n in range(StaticData.pFixed.shape[0]):
            self._calcp_FP_sparse[n, StaticData.pFixed[n, 0]] = 1
        
        self._calcp_M_S_sparse = csc_matrix(StaticData.M_S)
        #  End Pre-computation for calcp

        #  End Pre-factorization 

        # Initialize a sensible full temperature field and project to the ROM.
        T_fom_initial = np.ones(StaticData.nbNodes)*80
        for n in range(StaticData.nbNodes):
            if StaticData.Nodes[n,2]==1:
                T_fom_initial[n]=110
            if StaticData.Nodes[n,2]==-1:
                T_fom_initial[n]=60
        
        # Initialize reduced temperature state T_r by projecting T_fom_initial
        self.T_r = self.V_r.T @ T_fom_initial
        print(f"Initialized reduced state T_r with shape {self.T_r.shape}")

        # Full temperature state T is now reconstructed from T_r
        self.T = self.V_r @ self.T_r
        self.T_old = self.T.copy() # T_old should also be full state for convergence checks
        
        self.p=np.ones(StaticData.nbNodes)*10
        # ROM: Add T_r_old for convergence checks
        self.T_r_old = self.T_r.copy()
        self.viscosity=np.zeros((StaticData.nbPipes))
        self.rho=np.zeros((StaticData.nbPipes))
        self.cp=np.zeros((StaticData.nbPipes))
        self.viscosity_Pumps=np.zeros((StaticData.nbPumps))
        self.rho_Pumps=np.zeros((StaticData.nbPumps))
        self.rho_Supply=np.zeros((StaticData.nbSupplies))
        self.viscosity_Nodes=np.ones(StaticData.nbNodes)*0.0004
        self.rho_Nodes=np.ones(StaticData.nbNodes)*1000
        self.cp_Nodes=np.ones(StaticData.nbNodes)*4190
        self.cp_Pumps=np.zeros((StaticData.nbPumps))
        self.cp_consumer=np.zeros((StaticData.nbOutlets))
        self.cp_Supply=np.zeros((StaticData.nbSupplies))
        self.m_s=np.zeros(len(hydraulicMesh.leftEdges))
        self.m_s_old=np.zeros(len(hydraulicMesh.leftEdges))+1
        self.Pump_inactive=[]
        self.badPoint_violated=[]
        self.Pipes_lambda=np.zeros((StaticData.nbPipes))+0.03
        self.T_out=7
        self.T_mean_yesterday=7
        self.T_mean=7
        self.Ur=np.ones(StaticData.nbPipes) #heat transfer coeeficient
        self.datetimeStationary=[]
        self.datetimeDynamic=[]
        self.time_change_dynamic=datetime.timedelta(seconds=ProgramModes.dynamicTimestep) # timestep length dynamic thermal
        self.time_change_stationary=datetime.timedelta(minutes=ProgramModes.stationaryTimeStep) #timestep length staionary
        # if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
        #     self.time_change_stationary=datetime.timedelta(minutes=15)
        self.StandardConsumptionProfile=StandardConsumptionProfile()
        self.Q=np.zeros(StaticData.nbOutlets)
        self.denom=np.zeros(StaticData.nbOutlets)
        if ProgramModes.TempratureCalculationMode=='dynamic':
            print("Warning: ROMTimeSeries is initializing, but ProgramModes.TempratureCalculationMode is 'dynamic'. Dynamic ROM not fully implemented.")
            self.DynamicTemperature=DynamicTemperature(ProgramModes)
        self.j=0
        self.i=0
        
        self.thermal_conductivities = {
        "steel": 50,
        "insulation": 0.029,
        "casing": 0.40,
        "soil": 1
        }
        #Z_e = 0.2 #Abstand beider Leitungen 
        self.h_e = 1.2 # Überdeckungshöhe
    
        # Debug Logging for Sa/Se stability
        self._previous_buildSaSe_start_debug = None
        self._previous_buildSaSe_target_debug = None
        self._outer_thermic_iter_j_debug = 0 # To track j from within buildSaSe
    
    # Method to store the actual stat
    def save_state(self):
        return copy.deepcopy(self.__dict__)
    
    # Method to conda install y a saved state
    def restore_state(self, saved_state):
        self.__dict__ = copy.deepcopy(saved_state)
        
    def calcThermalResistanceVertical(self,StaticData):
        """Compute vertical thermal resistances for all pipes.

        This collects the individual layer contributions (steel, insulation,
        casing and soil) into the effective resistance used by the heat‑loss
        models. Values are stored in `self.Ur` and reused each timestep.
        """
        R_soil_ver=np.zeros(StaticData.nbPipes)
        R_service_pipe=np.zeros(StaticData.nbPipes)
        R_insulation=np.zeros(StaticData.nbPipes)
        R_casing=np.zeros(StaticData.nbPipes)
        for i in range(StaticData.nbPipes):     
            R_soil_ver[i]=(1 / (2 * np.pi * self.thermal_conductivities["soil"])) * np.arccosh((2 * self.h_e) / StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_cas_out"])
            R_service_pipe[i] = (1 / (2 * np.pi * self.thermal_conductivities["steel"])) * np.log(StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_st_out"] / StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_st_in"])
            R_insulation[i] = (1 / (2 * np.pi * self.thermal_conductivities["insulation"])) * np.log(StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_ins_out"] / StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_st_out"])
            R_casing[i] = (1 / (2 * np.pi * self.thermal_conductivities["casing"])) * np.log(StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_cas_out"] / StaticData.pipe_dimensions[StaticData.Pipes_DN[i]]["d_ins_out"])
        self.Ur=(R_service_pipe + R_insulation + R_casing + R_soil_ver)
    
    def datestart(self,ProgramModes):
        """set Start date of the calculation"""
        self.datetimeStationary.append(ProgramModes.startTime)
        if ProgramModes.TempratureCalculationMode=='dynamic':
            self.datetimeDynamic.append(ProgramModes.startTime)
           
    
    def timestepStationary(self,StaticData,hydraulicMesh,ProgramModes,LoadMeasValues, Tsupply,finalstep):
        """Run one stationary timestep of the hybrid ROM.

        The loop alternates between hydraulic (full‑order) and thermal (ROM)
        updates until convergence of the reduced temperature state `T_r`.
        """
        if self.datetimeStationary[-1]!=ProgramModes.startTime:
            self.datetimeStationary.append(self.datetimeStationary[-1]+self.time_change_stationary)

        self.j=0 # outer thermic loop counter
        self._outer_thermic_iter_j_debug = 0 # Reset for new timestepStationary call

        self.Pump_inactive=[]
        self.badPoint_violated=[]
        # self.T_old=np.ones(StaticData.nbNodes) # FOM version
        # ROM: T_old (full) will be updated inside the loop if needed for other parts
        # Initialize T_old consistent with current T (which is V_r @ T_r)
        self.T_old = self.T.copy() 

        if not(ProgramModes.PowerInput=='None'):
            self.Qstrom=LoadMeasValues.SupplyVolumeflow(self.datetimeStationary[-1])
            self.T_in_meas=LoadMeasValues.SupplyTemperature(self.datetimeStationary[-1])
            self.updateQ(StaticData,ProgramModes)
        else:
            self.updateQ(StaticData,ProgramModes)
        
        # ROM: Convergence check will be on T_r
        # old_T_norm_check = np.linalg.norm(self.T-self.T_old,2)
        # print(f"Initial T diff norm for check: {old_T_norm_check}, Target: {StaticData.nbNodes*ProgramModes.accuracy_factor_thermic}")
        
        # Use a small initial difference to ensure the loop runs at least once if T and T_old are identical initially due to setup
        # Or better, store T_r before the loop and compare T_r with its previous state.
        initial_Tr_for_loop_start = self.T_r.copy() 

        # while np.linalg.norm(self.T-self.T_old,2)>StaticData.nbNodes*ProgramModes.accuracy_factor_thermic and self.j<30: # FOM Condition
        # ROM Condition: Use T_r for convergence and r_dim for scaling the tolerance
        # Ensure T_old is updated for the first iteration's comparison or T_r_old is set before loop
        self.T_r_old = initial_Tr_for_loop_start - (ProgramModes.accuracy_factor_thermic * self.r_dim * 100) # Ensure it's different enough to start

        while np.linalg.norm(self.T_r - self.T_r_old, 2) > self.r_dim * ProgramModes.accuracy_factor_thermic and self.j < 30:
            self._outer_thermic_iter_j_debug = self.j # Update j for buildSaSe logging
            self.T_old = self.T.copy() # Store full T for parts of code that might still use T_old (e.g. if some internal checks were missed)
            self.T_r_old = self.T_r.copy() # Store previous T_r for convergence check
            
            self.waterProperties(StaticData) # Uses DEIM for rho, cp
            self.calcDeltaT(StaticData)
            
            self.calcmfromT(StaticData,ProgramModes)
            self.setL(StaticData)
            self.calcm0(hydraulicMesh,StaticData)
            self.calcm_i(hydraulicMesh,StaticData)
            self.m_s_old=np.zeros(len(hydraulicMesh.leftEdges))+1
            self.i=0
            while np.linalg.norm(self.m_s-self.m_s_old,2)>ProgramModes.accuracy_factor_hydraulic*len(hydraulicMesh.leftEdges) and self.i<30:
                
                self.hydraulicPrep(hydraulicMesh,StaticData,ProgramModes)
                #self.solveLinalgStep(hydraulicMesh,StaticData)
                self.NewtonRaphsonStep()
                self.i=self.i+1
            self.hydraulicPrep(hydraulicMesh,StaticData,ProgramModes)
            self.pressureCond(StaticData)
            self.findPipespg(StaticData)
            self.calcp(hydraulicMesh,StaticData)
            self.thermicPrep(StaticData,ProgramModes,Tsupply)
            
            self.calcT(StaticData)
            self.j=self.j+1
        self.calcPower(StaticData)
        if finalstep==False and self.datetimeStationary[-1]!=ProgramModes.startTime:
            del self.datetimeStationary[-1]
        elif finalstep==True and self.datetimeStationary[-1]==ProgramModes.startTime:
            self.datetimeStationary[-1] += datetime.timedelta(microseconds=1)

            

        
    def timestepDynamicHydraulic(self,StaticData,hydraulicMesh,ProgramModes,LoadMeasValues):
        """One quasi‑dynamic timestep for the hydraulic part (full‑order)."""
        
        if len(self.datetimeStationary)>1 or self.j>0:
            self.datetimeStationary.append(self.datetimeStationary[-1]+self.time_change_stationary)
            
        self.Pump_inactive=[]
        self.badPoint_violated=[]
        
        self.waterProperties(StaticData)

        if not(ProgramModes.PowerInput=='None'):
            self.Qstrom=LoadMeasValues.SupplyVolumeflow(self.datetimeStationary[-1])
            self.T_in_meas=LoadMeasValues.SupplyTemperature(self.datetimeStationary[-1])
            self.updateQ(StaticData,ProgramModes)
        else:
            self.updateQ(StaticData,ProgramModes)
        self.calcDeltaTdyn(StaticData)
        self.calcmfromT(StaticData,ProgramModes)
        self.setL(StaticData)
        self.calcm0(hydraulicMesh,StaticData)
        self.calcm_i(hydraulicMesh,StaticData)
        self.m_s_old=np.zeros(len(hydraulicMesh.leftEdges))+1
        self.i=0
        while np.linalg.norm(self.m_s-self.m_s_old,2)>ProgramModes.accuracy_factor_hydraulic*len(hydraulicMesh.leftEdges) and self.i<30:
            self.hydraulicPrep(hydraulicMesh,StaticData,ProgramModes)
            #self.solveLinalgStep(hydraulicMesh,StaticData)
            self.NewtonRaphsonStep()
            self.i=self.i+1
        self.hydraulicPrep(hydraulicMesh,StaticData,ProgramModes)
        self.pressureCond(StaticData)
        self.findPipespg(StaticData)
        self.calcp(hydraulicMesh,StaticData)
        self.calcPower(StaticData)
        
            
    def timestepDynamic(self,StaticData,ProgramModes,Tsupply):
        """Advance the dynamic thermal model by one step (experimental)."""
        if self.DynamicTemperature.step>0:
            self.T=self.DynamicTemperature.Timestep(StaticData,self,ProgramModes,Tsupply)
            self.datetimeDynamic.append(self.datetimeDynamic[-1]+self.time_change_dynamic)
        else:
            self.DynamicTemperature.Timestep(StaticData,self,ProgramModes,Tsupply)
            
    
    def updateQ(self,StaticData,ProgramModes):
        """Update outlet heat demands using the standard consumption profile."""
        self.Q=np.zeros(StaticData.nbOutlets)
        for n in range(StaticData.nbOutlets):
            if self.denom[n]==0:
                self.denom[n]=self.StandardConsumptionProfile.prepareConsumption(StaticData.Outlets_ConsumerType[n])
            T_weigth=(self.T_mean+0.5*self.T_mean_yesterday)/1.5
            self.Q[n]=self.StandardConsumptionProfile.calcConsumption(StaticData.Outlets_ConsumerType[n], self.denom[n],StaticData.Outlets_Qa[n],T_weigth,self.datetimeStationary[-1].weekday(),self.datetimeStationary[-1].hour)
        self.empty_consumer=np.where(self.Q==0)[0]
        #self.Q[self.empty_consumer]=StaticData.Outlets_Pmax[self.empty_consumer]*1000/200
        
        if ProgramModes.randomConsumerMode in ['Gauss']:
            for n in range(StaticData.nbOutlets):
                self.Q[n]=self.Q[n]*random.gauss(1,0.4)
        self.Q=np.clip(self.Q,StaticData.Outlets_Pmax/200*1000*0,StaticData.Outlets_Pmax*1000)
        #self.Q=StaticData.Outlets_Pmax*1000
        #self.Q[0]=5421891.729210106
        if ProgramModes.PowerInput=='Q':
            self.m_supply_meas=self.Qstrom*self.rho_Nodes[StaticData.vorlauf_SupplieNodes[0]]/3600
            self.Q_meas=(self.T_in_meas-60)*self.m_supply_meas*4190
            self.Qsum=sum(self.Q)
            self.Q_scale=self.Q_meas/self.Qsum
            self.Q=self.Q*self.Q_scale


        
    
    def calcmfromT(self,StaticData,ProgramModes):
        """Compute outlet mass flows from heat demand and temperature drop."""
        mask=self.Delta_T_consumer>0
        self.m_consumer=np.zeros_like(self.Q)
        self.m_consumer[mask]=self.Q[mask]/self.Delta_T_consumer[mask]/self.cp_consumer[mask]
        
        #self.m_consumer[self.empty_consumer]=0.01
        self.m_consumer=np.clip(self.m_consumer,0,500)
        self.m_supply=sum(self.m_consumer)
        
        if ProgramModes.PowerInput=='m':
            self.m_supply_meas=self.Qstrom*self.rho_Nodes[StaticData.vorlauf_SupplieNodes[0]]/3600
        
            self.m_scale=self.m_supply_meas/self.m_supply
            self.m_consumer=self.m_consumer*self.m_scale
        # self.m_consumer=np.clip(self.m_consumer,0,10)
        # self.m_supply=sum(self.m_consumer)
        
    
    def calcTfromm(self):
        self.Delta_T_consumer=self.Q/self.cp_consumer/self.m_consumer
        
    def setL(self,StaticData):
        """Assemble external mass‑flow vector from consumer demands.""" 
        self.L=np.zeros((StaticData.nbNodes))
        for n in range(StaticData.nbOutlets):
            self.L[StaticData.Outlets[n,1]]=self.L[StaticData.Outlets[n,1]]+self.m_consumer[n]
            self.L[StaticData.Outlets[n,2]]=self.L[StaticData.Outlets[n,2]]-self.m_consumer[n]
            
            
            # self.L[StaticData.vorlauf_SupplieNodes[1:]]=-70
            # self.L[StaticData.ruecklauf_SupplieNodes[1:]]=70
        
    def calcm0(self,hydraulicMesh,StaticData):
        """Solve tree (cycle‑free) mass flows using a pre‑factorized system.

        We reapply the same row deletions and masking as in the full model so
        that the ROM produces hydraulics consistent with the FOM structure.
        """
        
        # L is constructed based on current state (self.L)
        L_full = np.concatenate((self.L, np.zeros(len(hydraulicMesh.leftEdges) + len(StaticData.closedPipes))))
        
        # Apply the same row deletions to L that were applied to A during factorization
        # 1. Delete rows corresponding to specific supply nodes
        L_deleted_nodes = np.delete(L_full, self._calcm0_L_rows_to_delete_indices, axis=0)
        
        # 2. Delete rows that became all-zero in A (and thus L should also match this)
        # The mask self._calcm0_L_non_zero_row_mask was created based on A AFTER node deletion.
        # So, it should be applied to L_deleted_nodes.
        L_final = L_deleted_nodes[self._calcm0_L_non_zero_row_mask]

        if self._calcm0_A_factorization is not None:
            try:
                self.m0 = self._calcm0_A_factorization.solve(L_final)
            except Exception as e:
                print(f"ERROR in calcm0 using pre-factorized solve: {e}")
                print(f"Shape of L_final: {L_final.shape}, Factorized A expected compatible shape.")
                print("Falling back to original np.linalg.solve for this step (will be slow).")
                # Fallback to original method if factorized solve fails
                # This requires re-building A as in the original calcm0
                A_0_LE=np.zeros((len(hydraulicMesh.leftEdges),hydraulicMesh.A_Forest.shape[1]))
                for n in range(len(hydraulicMesh.leftEdges)): A_0_LE[n,hydraulicMesh.leftEdges[n]]=1
                A_0_Pipes=np.zeros((len(StaticData.closedPipes),hydraulicMesh.A_Forest.shape[1]))
                for n in range(len(StaticData.closedPipes)): A_0_Pipes[n,StaticData.closedPipes[n]]=1
                A_orig=np.concatenate((hydraulicMesh.A_Forest.toarray(),A_0_LE,A_0_Pipes)) # Assuming A_Forest is sparse
                L_orig_for_fallback=np.concatenate((self.L,np.zeros(len(hydraulicMesh.leftEdges)+len(StaticData.closedPipes))))
                idx_orig=np.array([StaticData.vorlauf_SupplieNodes[0],StaticData.ruecklauf_SupplieNodes[0]])
                L_orig_deleted_nodes = np.delete(L_orig_for_fallback,idx_orig,axis=0)
                A_orig_deleted_nodes = np.delete(A_orig,idx_orig,axis=0)
                L_orig_final = L_orig_deleted_nodes[~np.all(A_orig_deleted_nodes == 0, axis=1)]
                A_orig_final = A_orig_deleted_nodes[~np.all(A_orig_deleted_nodes == 0, axis=1)]
                self.m0=np.linalg.solve(A_orig_final, L_orig_final)
        else:
            # Factorization failed in __init__, use original method
            print("Pre-factorization of A for calcm0 failed in __init__. Using original np.linalg.solve (will be slow).")
            A_0_LE=np.zeros((len(hydraulicMesh.leftEdges),hydraulicMesh.A_Forest.shape[1]))
            for n in range(len(hydraulicMesh.leftEdges)): A_0_LE[n,hydraulicMesh.leftEdges[n]]=1
            A_0_Pipes=np.zeros((len(StaticData.closedPipes),hydraulicMesh.A_Forest.shape[1]))
            for n in range(len(StaticData.closedPipes)): A_0_Pipes[n,StaticData.closedPipes[n]]=1
            A_orig=np.concatenate((hydraulicMesh.A_Forest.toarray(),A_0_LE,A_0_Pipes)) # Assuming A_Forest is sparse
            L_orig_for_fallback=np.concatenate((self.L,np.zeros(len(hydraulicMesh.leftEdges)+len(StaticData.closedPipes))))
            idx_orig=np.array([StaticData.vorlauf_SupplieNodes[0],StaticData.ruecklauf_SupplieNodes[0]])
            L_orig_deleted_nodes = np.delete(L_orig_for_fallback,idx_orig,axis=0)
            A_orig_deleted_nodes = np.delete(A_orig,idx_orig,axis=0)
            L_orig_final = L_orig_deleted_nodes[~np.all(A_orig_deleted_nodes == 0, axis=1)]
            A_orig_final = A_orig_deleted_nodes[~np.all(A_orig_deleted_nodes == 0, axis=1)]
            self.m0=np.linalg.solve(A_orig_final, L_orig_final)

        # Original calcm0 content:
        # A_0_LE=np.zeros((len(hydraulicMesh.leftEdges),hydraulicMesh.A_Forest.shape[1]))
        # for n in range(len(hydraulicMesh.leftEdges)):
        #     A_0_LE[n,hydraulicMesh.leftEdges[n]]=1
            
        # A_0_Pipes=np.zeros((len(StaticData.closedPipes),hydraulicMesh.A_Forest.shape[1]))
        # for n in range(len(StaticData.closedPipes)):
        #     A_0_Pipes[n,StaticData.closedPipes[n]]=1

        #     # A_0_Test=np.zeros((1,hydraulicMesh.A_Forest.shape[1]))
        #     # A_0_Test[0,1]=1
        #    #A_0_Test[1,1742]=1
        
                
        # A=np.concatenate((hydraulicMesh.A_Forest,A_0_LE,A_0_Pipes))
        # L=np.concatenate((self.L,np.zeros(len(hydraulicMesh.leftEdges)+len(StaticData.closedPipes))))
        # #L=L[~np.all(A == 0, axis=1)]
        # #A = A[~np.all(A == 0, axis=1)]
        # idx=np.array([StaticData.vorlauf_SupplieNodes[0],StaticData.ruecklauf_SupplieNodes[0]])
                                                                                      
        # L = np.delete(L,idx,axis=0)
        # A = np.delete(A,idx,axis=0)
        
        # L = L[~np.all(A == 0, axis=1)]
        # A = A[~np.all(A == 0, axis=1)]
        
            
        # #A = A[:,~np.all(A == 0, axis=0)]
        # self.m0=np.linalg.solve(A,L)
        
        
    def calcDeltaT(self,StaticData):
        """Compute consumer temperature deltas for stationary mode."""
        self.Delta_T_consumer=self.T[StaticData.vorlauf_consumptionnode]-StaticData.Outlets_Tref
        
        # if np.any(self.Delta_T_consumer<=10):
        #     idx=np.where(self.Delta_T_consumer<=10)[0]
        #     print(str(idx))
        #     #self.Delta_T_consumer[idx]=10
        
    def calcDeltaTdyn(self,StaticData):
        """Compute consumer temperature deltas for dynamic mode."""
        self.Delta_T_consumer=self.T[StaticData.vorlauf_consumptionnode]-StaticData.Outlets_Tref
        

        
    def hydraulicPrep(self,hydraulicMesh,StaticData,ProgramModes):
        """Prepare hydraulic matrices for the next Newton step."""
        self.calcm_i(hydraulicMesh,StaticData)
        self.resistanceMatrix(StaticData,ProgramModes)
        self.calcF(hydraulicMesh,StaticData)
        self.calcJ(hydraulicMesh,StaticData)
        
    
    def calcm_i(self,hydraulicMesh,StaticData):
        """ Calculates the inner massflows according through the massflopws through the leftedges and the starting massflow
        splits it up into pipes and pumps"""
        self.m_i=self.m0+hydraulicMesh.B.transpose()@self.m_s
        self.m_p=self.m_i[StaticData.nbPipes:]
        self.m_i=self.m_i[:StaticData.nbPipes]
        
    def calcF(self,hydraulicMesh,StaticData): 
        """Assemble residual F (mesh pressure balance) for Newton iteration."""
        self.F=hydraulicMesh.B[:,:StaticData.nbPipes]@(self.R*self.m_i)
        
    def calcJ(self,hydraulicMesh,StaticData):
        """Build Jacobian of F (turbulent friction linearization)."""
        from scipy.sparse import diags
        B_pipes = hydraulicMesh.B[:, :StaticData.nbPipes]
        R_diag_sparse = diags(2 * self.R)
        self.J = (B_pipes @ R_diag_sparse) @ B_pipes.T
        
    def NewtonRaphsonStep(self):
        """Perform one Newton–Raphson step for mesh flows `m_s`."""
        self.m_s_old=self.m_s
        #self.m_s=self.m_s-np.linalg.inv(self.J)@self.F*0.01
        
        #delta_m=np.linalg.solve(-self.J,self.F)
        
        J_sparse = csc_matrix(self.J)
        F_sparse = csc_matrix(self.F).T
        
        delta_m = spsolve(-J_sparse, F_sparse)
    
        # self.m_s=self.m_s+delta_m
        
        #delta_m=np.linalg.solve(-self.J,self.F)
        self.m_s=self.m_s+delta_m

        
    def solveLinalgStep(self,hydraulicMesh,StaticData):
        """Alternative linear solve for the Newton step (legacy)."""
        self.m_s_old=self.m_s
        A=hydraulicMesh.B[:,:StaticData.nbPipes]@np.diag(self.R)@hydraulicMesh.B[:,:StaticData.nbPipes].transpose()
        B=-hydraulicMesh.B[:,:StaticData.nbPipes]@np.diag(self.R*2)@self.m0[:StaticData.nbPipes]
        self.m_s=np.linalg.solve(A,B)
        
    
        
    def pressureCond(self,StaticData):
        """Encode pump, fixed‑pressure and safety constraints for pressure solve."""
        self.BP=np.zeros((StaticData.nbPumps,StaticData.nbNodes))
        self.BPdelta=np.zeros(StaticData.nbPumps)
        for n in range(StaticData.nbPumps-1,-1,-1):
            if n in self.Pump_inactive:
                
                self.BP[n,StaticData.Pumps_startnode[n]]=1
                self.BP[n,StaticData.Pumps_targetnode[n]]=-1
                self.BPdelta[n]=0
                
                
            else:    
                filtered_array = StaticData.pBadPoints[StaticData.pBadPoints[:, 4] == n]
                differences= self.p[(filtered_array[:, 0]).astype(int)] - self.p[(filtered_array[:, 1]).astype(int)]-filtered_array[:, 2]
                try:
                    min_diff_index = np.argmin(differences)
                    row_min=filtered_array[min_diff_index]
                    self.BP[n,row_min[0]]=1
                    self.BP[n,row_min[1]]=-1
                    self.BPdelta[n]=row_min[2]
                except:
                    self.BP = np.delete(self.BP, n, axis=0)
                    self.BPdelta = np.delete(self.BPdelta, n, axis=0)
                    #print("1 Pump has no pressure to control")
            
        
        self.FP=np.zeros((StaticData.pFixed.shape[0],StaticData.nbNodes))
        self.FPValue=np.zeros(StaticData.pFixed.shape[0])
        for n in range(StaticData.pFixed.shape[0]):
            self.FP[n,StaticData.pFixed[n,0]]=1
            self.FPValue[n]=StaticData.pFixed[n,1]
            
    def findPipespg(self,StaticData):
        self.Pipes_pg=9.81*self.rho*StaticData.Pipes_h*1e-5
        
    def calcp(self,hydraulicMesh,StaticData):
        """Solve for nodal pressures subject to pump/fixed constraints.

        Builds a sparse linear system from network topology and constraints and
        solves it with a sparse solver. Pumps that violate flow direction and
        pressure safety constraints are iteratively deactivated until feasible.
        """
        from scipy.sparse import vstack
        while True:

            BP_sparse = csc_matrix(self.BP)
            
            # Use vstack for sparse matrices
            P_sparse = vstack([self._calcp_A_T_sparse, BP_sparse, self._calcp_FP_sparse, self._calcp_M_S_sparse])
            
            delta_p=np.delete(self.R*self.m_i-self.Pipes_pg,np.concatenate((hydraulicMesh.leftEdges,StaticData.closedPipes)).astype(int))
            
            P_X=np.concatenate((delta_p,self.BPdelta,self.FPValue,np.zeros(StaticData.nbSupplies))).astype(float)#.reshape(P.shape[0],1)
            #self.p=np.linalg.solve(P,P_X)
            # Convert the dense matrices to sparse CSC format
            
            P_X_sparse = csc_matrix(P_X).transpose()
            
            # Solve the sparse linear system
            self.p = spsolve(P_sparse, P_X_sparse)
            self.Pump_inactive=[]
            self.badPoint_violated=[]
            for n in range(StaticData.nbPumps):
                if self.p[StaticData.Pumps_startnode[n]]>self.p[StaticData.Pumps_targetnode[n]]:
                    self.Pump_inactive.append(n)
                    
            for n in range(StaticData.pBadPoints.shape[0]):
                if self.p[StaticData.pBadPoints[n,0]]<self.p[StaticData.pBadPoints[n,1]]+StaticData.pBadPoints[n,2]-0.001:
                    self.badPoint_violated.append(n)
            
            if not self.Pump_inactive and not self.badPoint_violated:
                self.p_bad=self.p[StaticData.pBadPoints[n,0]]<self.p[StaticData.pBadPoints[n,1]]
                break
            
            self.pressureCond(StaticData)
                
            
                      
        
                
                
                
                
        
        idx=np.where(self.p<=1)[0]
        if np.any(self.p<=1):
            print('pressure smaller than 1 bar at nodes ' + str(idx))
            self.p[idx]=1
            
        idx=np.where(self.p>=22)[0]
        if np.any(self.p<=1):
            print('pressure larger than 22 bar at nodes ' + str(idx))
            self.p[idx]=22

                
        
        
    def thermicPrep(self,StaticData,ProgramModes,Tsupply):
        """Prepare matrices/vectors for the stationary thermal solve."""
        self.buildSaSe(StaticData)
        self.buildTemperatureMatrix(StaticData,ProgramModes,Tsupply)
        self.buildsetT(StaticData,Tsupply)
        
            
        
    
    def resistanceMatrix(self,StaticData,ProgramModes):
        """Update hydraulic resistances based on current state.

        Uses Reynolds number and a Haaland‑type friction approximation to keep
        the computational cost low while remaining sufficiently accurate.
        """
        self.waterProperties(StaticData) #update Waterproteries maybe change?
        self.v= self.m_i/self.rho/StaticData.Pipes_A #flow velocity v=m/rho/A
        self.Re= self.rho*abs(self.v)*StaticData.Pipes_d/self.viscosity #Reynoldszahl rho/*|v|*d/viscosity
        self.calcLambda(StaticData,ProgramModes.PipeFrictionMethod) #calculate lambda colebrook
        self.R=(self.Pipes_lambda*StaticData.Pipes_cr*np.clip(abs(self.m_i),1,None)*1e-5)/self.rho #calculate  R vector R=lambda*cr/rho*|mi|*1e-5 convert to bar
        
    def waterProperties(self,StaticData):
        """Estimate water properties efficiently with DEIM.

        Instead of evaluating properties at every node, we query the steam
        table only at a few interpolation points and project the results back
        to all nodes through the DEIM projection matrices. This preserves the
        nonlinear dependence on pressure/temperature at a fraction of the cost.
        """
        # self.T is the full temperature state, reconstructed from self.T_r = self.V_r @ self.T_r
        
        # Viscosity (viscosity_Nodes): Approximated using DEIM
        viscosity_at_deim_indices = np.zeros(self.m_viscosity)
        for i, node_idx in enumerate(self.deim_viscosity_indices):
            viscosity_at_deim_indices[i] = self.steamTable.my_pt(self.p[node_idx], self.T[node_idx])
        self.viscosity_Nodes = self.deim_viscosity_proj_matrix @ viscosity_at_deim_indices


        # Density (rho_Nodes): Approximated using DEIM
        rho_at_deim_indices = np.zeros(self.m_rho)
        for i, node_idx in enumerate(self.deim_rho_indices):
            rho_at_deim_indices[i] = self.steamTable.rho_pt(self.p[node_idx], self.T[node_idx])
        self.rho_Nodes = self.deim_rho_proj_matrix @ rho_at_deim_indices
        if self.rho_Nodes.ndim == 1 and StaticData.nbNodes > 1 and self.deim_rho_proj_matrix.shape[0] == StaticData.nbNodes:
             pass 
        elif self.rho_Nodes.shape[0] != StaticData.nbNodes:
            if self.rho_Nodes.size == StaticData.nbNodes:
                self.rho_Nodes = self.rho_Nodes.reshape(StaticData.nbNodes)

        # Specific Heat Capacity (cp_Nodes) - ROM's DEIM version for general use (e.g., in W advection term for calcT)
        cp_at_deim_indices = np.zeros(self.m_cp)
        for i, node_idx in enumerate(self.deim_cp_indices):
            cp_at_deim_indices[i] = self.steamTable.Cp_pt(self.p[node_idx], self.T[node_idx]) * 1000 
        self.cp_Nodes = self.deim_cp_proj_matrix @ cp_at_deim_indices
        if self.cp_Nodes.ndim == 1 and StaticData.nbNodes > 1 and self.deim_cp_proj_matrix.shape[0] == StaticData.nbNodes:
            pass 
        elif self.cp_Nodes.shape[0] != StaticData.nbNodes:
            if self.cp_Nodes.size == StaticData.nbNodes:
                self.cp_Nodes = self.cp_Nodes.reshape(StaticData.nbNodes)

        self.viscosity = _average_properties_numba(self.viscosity_Nodes, StaticData.startnode, StaticData.targetnode)
        self.rho = _average_properties_numba(self.rho_Nodes, StaticData.startnode, StaticData.targetnode)
        self.cp = _average_properties_numba(self.cp_Nodes, StaticData.startnode, StaticData.targetnode)
            
        self.viscosity_Pumps = _average_properties_numba(self.viscosity_Nodes, StaticData.Pumps_startnode, StaticData.Pumps_startnode)
        self.rho_Pumps = _average_properties_numba(self.rho_Nodes, StaticData.Pumps_startnode, StaticData.Pumps_startnode)
        self.cp_Pumps = _average_properties_numba(self.cp_Nodes, StaticData.Pumps_startnode, StaticData.Pumps_startnode)
            
        self.cp_consumer = _average_properties_numba(self.cp_Nodes, StaticData.vorlauf_consumptionnode, StaticData.ruecklauf_consumptionnode)

        self.cp_Supply = _average_properties_numba(self.cp_Nodes, StaticData.vorlauf_SupplieNodes, StaticData.ruecklauf_SupplieNodes)
        self.rho_Supply = _average_properties_numba(self.rho_Nodes, StaticData.vorlauf_SupplieNodes, StaticData.ruecklauf_SupplieNodes)
            
    def calcLambda(self,StaticData,PipeFrictionMethod):
        """calculate the Pipe friction coefficient according to the Reynold-Number etc.
        Different Formulas for the turbulent calculation implemented"""
        self.Pipes_lambda = _calc_lambda_numba(self.Re, StaticData.Pipes_k, StaticData.Pipes_d, self.Pipes_lambda)
        
    def buildSaSe(self,StaticData):
        """Assemble selector matrices for upwind advection.

        `Se` picks the temperature at an element's inflow node and `Sa` at the
        outflow node. Rows tied to nearly zero mass flow are zeroed to avoid
        numerical noise and keep the operators well conditioned.
        """
        from scipy.sparse import coo_matrix, diags
        start=np.concatenate((StaticData.startnode,StaticData.Pumps_startnode,StaticData.vorlauf_consumptionnode, StaticData.ruecklauf_SupplieNodes)) #startnodes Pipes and Pumps ,Outlets and Supplie
        target=np.concatenate((StaticData.targetnode,StaticData.Pumps_targetnode, StaticData.ruecklauf_consumptionnode,StaticData.vorlauf_SupplieNodes)) #targetnodes Pipes and Pumps,,Outlets and Supplie
        idx=np.where(self.m_i<0) #idx with negative massflow --> temperature change in the other direction
        if len(idx[0])>0:
            start[idx],target[idx]=target[idx],start[idx] # switch Massflow direction
        
        # --- Debug Logging for Sa/Se stability ---
        if self._previous_buildSaSe_start_debug is not None and self._previous_buildSaSe_target_debug is not None:
            if (np.array_equal(self._previous_buildSaSe_start_debug, start) and 
                np.array_equal(self._previous_buildSaSe_target_debug, target)):
                print(f"    [j={self._outer_thermic_iter_j_debug}] buildSaSe: start/target arrays UNCHANGED from previous call.")
            else:
                print(f"    [j={self._outer_thermic_iter_j_debug}] buildSaSe: start/target arrays CHANGED.")
        else:
            print(f"    [j={self._outer_thermic_iter_j_debug}] buildSaSe: Initial population or no previous start/target data.")
        
        self._previous_buildSaSe_start_debug = start.copy()
        self._previous_buildSaSe_target_debug = target.copy()
        # --- End Debug Logging ---

        num_edges = start.shape[0]
        edge_indices = np.arange(num_edges)
        data = np.ones(num_edges)

        # Create sparse matrices
        self.Se = coo_matrix((data, (edge_indices, start)), shape=(num_edges, StaticData.nbNodes)).tocsr()
        self.Sa = coo_matrix((data, (edge_indices, target)), shape=(num_edges, StaticData.nbNodes)).tocsr()

        idxdel=np.where(abs(self.m_i)<=1e-6)[0]
        
        # Efficiently zero out rows in CSR matrix
        if len(idxdel) > 0:
            # Create a multiplier vector that is 0 for rows to be deleted, 1 otherwise
            row_multiplier = np.ones(num_edges)
            row_multiplier[idxdel] = 0
            # Create a sparse diagonal matrix from the multiplier
            M = diags(row_multiplier)
            # Multiply the matrices to zero out the desired rows, preserving shape
            self.Sa = M @ self.Sa
            self.Se = M @ self.Se
            
    def buildTemperatureMatrix(self,StaticData,ProgramModes,Tsupply):
        """Build advection/production terms for the thermal ROM.

        Creates the vectors/matrices entering the reduced system, including the
        advection weights `W`, boundary/source vector `D`, and (mode‑dependent)
        heat‑loss terms. Pipe heat losses are approximated via DEIM.
        """
        # self.W is calculated using the ROM's DEIM-based self.cp, self.cp_Pumps, self.cp_consumer, self.cp_Supply for advection consistency with POD basis
        self.W=abs(np.concatenate((self.m_i,self.m_p,self.m_consumer,np.atleast_1d(self.m_supply)),axis=0))*np.concatenate((self.cp, self.cp_Pumps,self.cp_consumer,np.atleast_1d(self.cp_Supply)),axis=0)
        
        # T_pipe_mean is calculated based on the full self.T (reconstructed from T_r)
        self.T_pipe_mean=np.zeros(StaticData.nbPipes)
        for n in range(StaticData.nbPipes):
            self.T_pipe_mean[n]=1/2*(self.T[StaticData.startnode[n]]+self.T[StaticData.targetnode[n]])

        # --- FULL ROM: Calculating Q_loss using DEIM approximation ---
        # print(f"    [ROM DBG, j={self._outer_thermic_iter_j_debug}] Calculating Q_loss_Pipe/Q_loss (W/m) using DEIM. Mode: {ProgramModes.HeatLossMode}")

        if ProgramModes.HeatLossMode in [ 'default','simple']:
            # Calculate Q_loss_Pipe at DEIM indices
            q_loss_pipe_at_deim_indices = np.zeros(self.m_q_loss)
            for i, pipe_idx in enumerate(self.deim_q_loss_indices):
                q_loss_pipe_at_deim_indices[i] = StaticData.Pipes_l[pipe_idx] / self.Ur[pipe_idx] * (self.T_pipe_mean[pipe_idx] - self.T_out)
            
<<<<<<< HEAD
            # Approximate full Q_loss_Pipe using DEIM projection and apply correction
            self.Q_loss_Pipe = self.deim_q_loss_proj_matrix @ q_loss_pipe_at_deim_indices
            self.Q_loss_Pipe = self.Q_loss_Pipe * 5.5787665
=======
            # Approximate full Q_loss_Pipe using DEIM projection
            self.Q_loss_Pipe = self.deim_q_loss_proj_matrix @ q_loss_pipe_at_deim_indices
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b
            
            tempexp=np.exp(-StaticData.Pipes_l/
                           (StaticData.Pipes_A*self.cp*self.rho*self.Ur*np.clip(abs(self.v),0.001,np.inf))).astype(float) # e^(-Ur*l/d/cp/rho/v)
            
            C_Outlets=np.ones(StaticData.nbOutlets)
            # Predefine D_Outlets as zeros
            D_Outlets = np.zeros_like(self.Q)
            
            # Perform the calculation only where m_consumer is non-zero
            non_zero_mask = self.m_consumer != 0
            D_Outlets[non_zero_mask] = -self.Q[non_zero_mask] / (self.m_consumer[non_zero_mask] * self.cp_consumer[non_zero_mask])
            
            C_Supply=np.zeros(StaticData.nbSupplies)
            D_Supply=Tsupply
                
            self.C_diag=np.concatenate((tempexp,np.ones(StaticData.nbPumps),C_Outlets,C_Supply),axis=0)
            
            self.D=np.concatenate((self.T_out*(1-tempexp),np.zeros(StaticData.nbPumps),D_Outlets,D_Supply),axis=0) #D Vektor T
        
        if ProgramModes.HeatLossMode not in [ 'default','simple']:
            T_supply_avg = 0 
            T_return_avg = 0 
            if ProgramModes.HeatLossNextPipe in ['No']:
                if StaticData.nbNodes > 0 and len(StaticData.supplyNodes) > 0 : T_supply_avg=np.average(self.T[StaticData.supplyNodes])
                if StaticData.nbNodes > 0 and len(StaticData.returnNodes) > 0 : T_return_avg=np.average(self.T[StaticData.returnNodes])
            
            q_loss_at_deim_indices = np.zeros(self.m_q_loss)

            for i, pipe_idx in enumerate(self.deim_q_loss_indices):
                current_T_pipe_mean = self.T_pipe_mean[pipe_idx]
                T_next_pipe = 0
                
                if ProgramModes.HeatLossMode  in ['EMCM']:
                    if ProgramModes.HeatLossNextPipe in ['No']:
                        if pipe_idx in StaticData.supply_pipes:
                            q_loss_at_deim_indices[i] = self.heatLossEMCM(current_T_pipe_mean, T_return_avg, self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[pipe_idx]])
                        elif pipe_idx in StaticData.return_pipes:
                            q_loss_at_deim_indices[i] = self.heatLossEMCM(current_T_pipe_mean, T_supply_avg, self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[pipe_idx]])
                        else:
                            q_loss_at_deim_indices[i] = 0 
                    
                elif ProgramModes.HeatLossMode  in ['Norm']:
                    if ProgramModes.HeatLossNextPipe in ['No']:
                        if pipe_idx in StaticData.supply_pipes:
                            q_loss_at_deim_indices[i] = self.heatLossNorm(current_T_pipe_mean, T_return_avg, self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[pipe_idx]]["d_ins_out"], StaticData.pipe_dimensions[StaticData.Pipes_DN[pipe_idx]]["d_st_out"])
                        elif pipe_idx in StaticData.return_pipes:
                            q_loss_at_deim_indices[i] = self.heatLossNorm(current_T_pipe_mean, T_supply_avg, self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[pipe_idx]]["d_ins_out"], StaticData.pipe_dimensions[StaticData.Pipes_DN[pipe_idx]]["d_st_out"])
                        else:
                            q_loss_at_deim_indices[i] = 0
            
            self.Q_loss = self.deim_q_loss_proj_matrix @ q_loss_at_deim_indices
            
            self.Q_loss_Pipe = self.Q_loss * StaticData.Pipes_l
            
<<<<<<< HEAD
            # Apply ROM-specific correction factor
            self.Q_loss_Pipe = self.Q_loss_Pipe * 5.5787665

=======
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b
            C_Outlets=np.ones(StaticData.nbOutlets)

            D_Outlets = np.zeros_like(self.Q)
            
            non_zero_mask = self.m_consumer != 0
            D_Outlets[non_zero_mask] = -self.Q[non_zero_mask] / (self.m_consumer[non_zero_mask] * self.cp_consumer[non_zero_mask])
            
            C_Supply=np.zeros(StaticData.nbSupplies)
            D_Supply=Tsupply
    
            self.C_diag=np.concatenate((np.ones(StaticData.nbPipes+StaticData.nbPumps),C_Outlets,C_Supply),axis=0)
<<<<<<< HEAD
            # Use the corrected Q_loss_Pipe to calculate DPipes for the thermal solve
            DPipes=-self.Q_loss_Pipe/(StaticData.Pipes_A*self.cp*self.rho*np.clip(abs(self.v),0.001,np.inf))
=======
            DPipes=-self.Q_loss*StaticData.Pipes_l/(StaticData.Pipes_A*self.cp*self.rho*np.clip(abs(self.v),0.001,np.inf))
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b
            self.D=np.concatenate((DPipes,np.zeros(StaticData.nbPumps),D_Outlets,D_Supply),axis=0)
                
    def buildsetT(self,StaticData,Tsupply):
        """set the T values at the input of the networks (supply in the supply network and consumer at the return network) """
        setTsupply=np.zeros((StaticData.nbSupplies,StaticData.nbNodes))
        
    def calcT(self,StaticData):
        """Solve the reduced stationary temperature system.

        Projects the full advection operators to the reduced space to form the
        small system `T_A_r @ T_r = T_B_r`, solves it, and reconstructs the
        full temperature field as `T = V_r @ T_r`.
        """
        # Pass CSR matrix components to Numba
        T_A_r, T_B_r = _build_and_project_thermic_system_numba(
            self.Sa.data, self.Sa.indices, self.Sa.indptr,
            self.Se.data, self.Se.indices, self.Se.indptr,
            self.V_r, self.W, self.C_diag, self.D
        )
        
        try:
            self.T_r = lsqr(T_A_r, T_B_r, x0=self.T_r.copy(), atol=1e-12, btol=1e-12)[0]
        except Exception as e:
            print(f"Error during lsqr solve for T_r in ROMTimeSeries.calcT: {e}")
            print(f"Shape of T_A_r: {T_A_r.shape}, Shape of T_B_r: {T_B_r.shape}")
            raise

        self.T = self.V_r @ self.T_r
        
        self.T=np.clip(self.T,30,160)
        
    def heatLossNorm(self,T_f,T_r,T_outside,D_i,d_o):
        #Summer
        #T_f = 80  # Supply temperature (°C)
        #T_r = 60  # Return temperature (°C)
        #T_s = 16  # Undisturbed Groundtemperature in depth Z (°C)
        #lambda_g = 1  # Heat conductivity of the ground (W/(mK))
        
        #Winter
        #T_f = 110  # Supply temperature (°C)
        #T_r = 60  # Return temperature (°C)
        T_s = 4  # Undisturbed Groundtemperature in depth Z (°C)
        lambda_g = self.thermal_conductivities["soil"] # Heat conductivity of the ground (W/(mK))
        
        # Model parameters
        lambda_i = self.thermal_conductivities["insulation"]  # Heat conductivity of the insulation (W/(mK))
        
        #DN50
        D_i = 1000*D_i  
        d_o = 1000*d_o  #from mm in m
        
        #DN150
        #D_i = 242.8  # Outer insulation diameter (mm)
        #d_o = 168.3  # Inner insulation diameter (mm)
        
        #DN400
        #D_i = 548.6  # Outer insulation diameter (mm)
        #d_o = 406.4  # Inner insulation diameter (mm)
        
        
        C = d_o + 200  # Distance between the pipe axes (mm)
        Z = d_o/2 + 1200  # Buried depth (mm)
        R_o = 0.0685
        Z_c = Z + R_o * lambda_g # corrected depth (mm)
        
        # Constants
        pi = math.pi
        
        # Calculations
        # Symmetrical and anti-symmetrical heat losses (Equations 3 and 4)
        
        T_S = (T_f + T_r) / 2 # Symmertrical Temperature
        
        T_a = (T_f - T_r) / 2 # Antisymmertrical Temperature
        
        beta = (lambda_g / lambda_i) * math.log(D_i / d_o) 
        
        
        h_sym = 1 /  (math.log(4 * Z_c / D_i) + beta + math.log(math.sqrt(1 + (2 * Z_c / C)**2))) # Symmetrical heatloss factor
        
        h_a = 1 / (math.log(4 * Z_c / D_i) + beta - math.log(math.sqrt(1 + (2 * Z_c / C)**2))) # Antisymmetrical heatloss factor
        
        #R_sym = (1/ 2* pi* lambda_g) *(math.log(4 * Z_c / D_i) + beta + math.log(math.sqrt(1 + (2 * Z_c / C)**2)))
        #R_a = (1/ 2* pi* lambda_g) *(math.log(4 * Z_c / D_i) + beta - math.log(math.sqrt(1 + (2 * Z_c / C)**2)))
        
        
        q_sym = (T_S - T_s) * 2 * pi * lambda_g * h_sym # symmetrical component of the heat loss
        
        q_a = T_a * 2 * pi * lambda_g * h_a # antisymmetrical component of the heat loss
        
        
        q_f = q_sym + q_a # heat loss of the flow pipe
        
        return q_f
    
    def heatLossEMCM(self,T_flow,T_return,T_out,flow_pipe_dims):
        # Define winter temperatures
        #T_out = 0      # Outside Air Temperature
        #T_flow = 110    # Flow Temperature
        #T_return = 80   # Return Temperature
        
        # # Define summer temperatures
        # T_out = 0      # Outside Air Temperature
        # T_flow = 110    # Flow Temperature
        # T_return = 50   # Return Temperature
        
        # Define constants
        Z_e = 0.2
        h_e = 1.2
        
        # # Define pipe geometries DN50
        # flow_pipe_dims = {
        #    "d_st_in": 0.0545,
        #     "d_st_out": 0.0603,
        #     "d_ins_out": 0.119,
        #     "d_cas_out": 0.125
        # }
        
        
        
        

        
        # Calculate thermal resistance of soil vertically and horizontally
        R_soil_hor = Z_e / self.thermal_conductivities["soil"]
        R_soil_ver = (1 / (2 * np.pi * self.thermal_conductivities["soil"])) * np.arccosh((2 * h_e) / flow_pipe_dims["d_cas_out"])
        
        # Calculate thermal resistance of flow pipe components
        R_service_pipe = (1 / (2 * np.pi * self.thermal_conductivities["steel"])) * np.log(flow_pipe_dims["d_st_out"] / flow_pipe_dims["d_st_in"])
        R_insulation = (1 / (2 * np.pi * self.thermal_conductivities["insulation"])) * np.log(flow_pipe_dims["d_ins_out"] / flow_pipe_dims["d_st_out"])
        R_casing = (1 / (2 * np.pi * self.thermal_conductivities["casing"])) * np.log(flow_pipe_dims["d_cas_out"] / flow_pipe_dims["d_ins_out"])
        
        # Calculate thermal resistance of return pipe components
        # R_service_pipe_return = (1 / (2 * np.pi * thermal_conductivities["steel"])) * np.log(return_pipe_dims["d_st_out"] / return_pipe_dims["d_st_in"])
        # R_insulation_return = (1 / (2 * np.pi * thermal_conductivities["insulation"])) * np.log(return_pipe_dims["d_ins_out"] / return_pipe_dims["d_st_out"])
        # R_casing_return = (1 / (2 * np.pi * thermal_conductivities["casing"])) * np.log(return_pipe_dims["d_cas_out"] / return_pipe_dims["d_ins_out"])
        
        # Equivalent Mesh Current Method
        R1 = R_service_pipe + R_insulation + R_casing + R_soil_ver
        #R2 = R1
        R3 = 2*(R_service_pipe + R_insulation + R_casing) + R_soil_hor 
        
        # Save the delta temperatures as variables
        U1 = T_flow - T_out   # Difference in temperature between flow and ambient
        U2 = T_return - T_out   # Difference in temperature between return and ambient
        #U3 = 0
        
        # Aaron Build the matrices row by row
        # r_matrix = np.array([
        #     [R1, 0, -R1],
        #     [0, R3, -R3],
        #     [-R1, -R3, R1 + R3 + R2]
        # ])
        
        # Aaron Temperatures matrix
        #v_matrix = np.array([U1, U2 - U1, U3])
        
        # Solve for currents (I) using Ohm's Law: U = R * I
        #i_matrix = np.linalg.solve(r_matrix, v_matrix)
        
        #Iloss_flow = i_matrix[0] - i_matrix[1]
        I=(1/R1+(1-U2/U1)/R3)*U1
        #Iloss_return = i_matrix[1] 
        #Total = Iloss_flow + Iloss_return
        
        return I
    
    def calcPumpPower(self,m,p,m_0,p_0,n_0,eta_0,exp,rho):
        n=n_0*(m/m_0*(p/p_0)**0.5)**0.5
        eta=1-(1-eta_0)*(n/n_0)**exp
        P_hyd=m*p/rho*1e5
        P_mech=P_hyd/eta
        P_elek=P_mech/0.95
        return P_elek
        
    
    def calcPower(self,StaticData):
        #self.PumpPower=self.m_p*(self.p[StaticData.Pumps_targetnode]-self.p[StaticData.Pumps_startnode])*1e5/self.rho_Pumps
        self.hydraulicPower=self.calcPumpPower(self.m_p,\
                                               self.p[StaticData.Pumps_targetnode]-self.p[StaticData.Pumps_startnode],\
                                               StaticData.Pumps_m0, StaticData.Pumps_p0,\
                                               StaticData.Pumps_n0,  StaticData.Pumps_eta0,\
                                               0.15, self.rho_Pumps)
        #if np.any((self.p <= 1) | (self.p >= 22)):
        #    self.hydraulicPower=self.hydraulicPower*100
        #cp_Supply=self.cp[StaticData]
        self.HeatPowerSupply=self.m_supply*(self.T[StaticData.vorlauf_SupplieNodes]-self.T[StaticData.ruecklauf_SupplieNodes])*self.cp_Supply
<<<<<<< HEAD
        # Q_loss_Pipe is now corrected inside buildTemperatureMatrix
        self.HeatLoses = np.sum(self.Q_loss_Pipe)
=======
        self.HeatLoses=sum(self.Q_loss_Pipe)*5.5787665
>>>>>>> 8a943aff3b95ac598887419c967bc5fe0537570b
        #self.HeatLosesTest=sum(self.Q_loss*StaticData.Pipes_l)
        #self.HeatLosesTest2=sum(abs(self.m_i*(self.T[StaticData.startnode]*self.cp_Nodes[StaticData.startnode]-self.T[StaticData.targetnode]*self.cp_Nodes[StaticData.targetnode])))
        #self.HeatLosesTest2=sum(self.m_i*(self.T[StaticData.startnode]-self.T[StaticData.targetnode])*self.cp)

