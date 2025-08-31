       # -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:30:41 2023

@author: sander
"""
import pandas as pd
import numpy as np
import math


class StaticData:
    def __init__(self):
        pass

    def csv2numpy(self, x):
        """Read a CSV file and convert it to a NumPy array."""
        try:
            # Read CSV file without a header and using ';' as delimiter
            temp = pd.read_csv(r'FreimannTestVollständig/LEACaGE/staticdata/' + x + '.csv', header=None, delimiter=';')
        except pd.errors.EmptyDataError:
            # If the CSV file is empty, return an empty array with shape (0, 20)
            temp = np.empty((0, 20))
        try:
            # Convert the DataFrame to a NumPy array
            temp = temp.to_numpy()
        except Exception:
            pass
        return temp

    def readData(self, ProgramModes):
        """Read all static data based on the system mode."""
        if ProgramModes.SystemMode in ['default', 'StateEstimation', 'TimeSeries']:
            # Read various CSV files into NumPy arrays
            self.mMeasurements = self.csv2numpy('mMeasurements')
            self.Nodes = self.csv2numpy('Nodes')
            self.OutletsE = self.csv2numpy('OutletsE')
            self.OutletsM = self.csv2numpy('OutletsM')
            self.Outlets = self.csv2numpy('Outlets')
            self.Pipes = self.csv2numpy('Pipes')
            self.pMeasurements = self.csv2numpy('pMeasurements')
            self.Pumps = self.csv2numpy('Pumps')
            self.Sliders = self.csv2numpy('Sliders')
            self.Supplies = self.csv2numpy('Supplies')
            self.TMeasurements = self.csv2numpy('TMeasurements')
            self.Valves = self.csv2numpy('Valves')
            self.pBadPoints = self.csv2numpy('pBadPoints')
            self.pFixed = self.csv2numpy('pFixed')
            self.ValveStates = self.csv2numpy('ValveStates')
            self.SliderStates = self.csv2numpy('SliderStates')

        if ProgramModes.SystemMode == 'NetworkSimulation':
            # Additional data reading for NetworkSimulation mode
            self.mMeasurements = self.csv2numpy('mMeasurements')
            self.Nodes = self.csv2numpy('Nodes')
            self.OutletsE = self.csv2numpy('OutletsE')
            self.OutletsM = self.csv2numpy('OutletsM')
            self.Pipes = self.csv2numpy('Pipes')
            self.pMeasurements = self.csv2numpy('pMeasurements')
            self.Pumps = self.csv2numpy('Pumps')
            self.Sliders = self.csv2numpy('Sliders')
            self.Supplies = self.csv2numpy('Supplies')
            self.TMeasurements = self.csv2numpy('TMeasurements')
            self.Valves = self.csv2numpy('Valves')

    def convertAllData(self):
        """Convert all data."""
        self.findnb()
        self.data_to_system_matrix()
        self.data_to_pump_matrix()
        self.data_to_consumption_matrix()
        self.data_to_consumption()
        self.data_to_valve_matrix()
        self.data_to_supply_matrix()
        self.data_to_null_node_matrix()
        self.data_to_meas_network_matrices()
        self.pipes_data()
        self.pump_data()
        self.supply_or_return()
        self.outlets_data()
        self.valve_state()
        self.DNdims()
        self.findDNpipe()

    def findnb(self):
        """Find the number of components and store in variables."""
        self.nbNodes = self.Nodes.shape[0]
        self.nbPipes = self.Pipes.shape[0]
        self.nbPumps = self.Pumps.shape[0]
        self.nbOutletsE = self.OutletsE.shape[0]
        self.nbOutletsM = self.OutletsM.shape[0]
        self.nbOutlets = self.nbOutletsE + self.nbOutletsM
        self.nbSupplies = self.Supplies.shape[0]
        self.nbValves = self.Valves.shape[0] + self.Sliders.shape[0]
        self.nbComp = self.nbPipes + self.nbValves + self.nbPumps + self.nbOutlets + self.nbSupplies
        self.nbEdges = self.nbPipes + self.nbPumps
        self.nbpMeas = self.pMeasurements.shape[0]
        self.nbmMeas = self.mMeasurements.shape[0]
        self.nbTmeas = self.TMeasurements.shape[0]

    def supply_or_return(self):
        """Find supply or return nodes."""
        self.supplyNodes = np.where(self.Nodes[:, 2] == 1)[0]
        self.returnNodes = np.where(self.Nodes[:, 2] == -1)[0]
        
        self.supply_pipes = []
        self.return_pipes = []
        for n in range(self.nbPipes):
            if self.startnode[n] in self.supplyNodes:
                self.supply_pipes.append(n)
            if self.startnode[n] in self.returnNodes:
                self.return_pipes.append(n)
        self.supply_pipes = np.array(self.supply_pipes)
        self.return_pipes = np.array(self.return_pipes)
                
                
         

    def data_to_system_matrix(self):
        """Build system matrix M from data."""
        self.startnode = self.Pipes[:, 0].astype(int)
        self.targetnode = self.Pipes[:, 1].astype(int)
        self.M_p = np.zeros((self.nbNodes, self.nbPipes))
        self.M_n = np.zeros((self.nbNodes, self.nbPipes))
        for n in range(self.nbPipes):
            self.M_n[int(self.startnode[n]), n] = 1
            self.M_p[int(self.targetnode[n]), n] = 1
        self.M = self.M_p - self.M_n

    def data_to_pump_matrix(self):
        """Build pump matrix P from data."""
        self.Pumps_startnode = self.Pumps[:, 1].astype(int)
        self.Pumps_targetnode = self.Pumps[:, 2].astype(int)
        self.M_Pumps_p = np.zeros((self.nbPumps, self.nbNodes))
        self.M_Pumps_n = np.zeros((self.nbPumps, self.nbNodes))
        self.M_P = np.zeros((2 * self.nbPumps, self.nbNodes))
        for n in range(self.nbPumps):
            self.M_Pumps_n[n, self.Pumps_startnode[n]] = 1
            self.M_Pumps_p[n, self.Pumps_targetnode[n]] = 1
            self.M_P[2 * n, self.Pumps_startnode[n]] = 1
            self.M_P[2 * n + 1, self.Pumps_targetnode[n]] = -1
        self.M_Pumps = self.M_Pumps_p + self.M_Pumps_n
        self.M_Pumps2 = self.M_Pumps_p - self.M_Pumps_n

    def data_to_consumption_matrix(self):
        """Build consumption/outlets matrices M_w and A from data."""
        self.vorlauf_consumptionnode = np.concatenate((self.OutletsM[:, 1], self.OutletsE[:, 1])).astype(int)
        self.ruecklauf_consumptionnode = np.concatenate((self.OutletsM[:, 5], self.OutletsE[:, 2])).astype(int)
        self.A = np.zeros((2 * self.nbOutlets, self.nbNodes))
        for n in range(self.nbOutlets):
            self.A[2 * n, self.vorlauf_consumptionnode[n]] = 1
            self.A[2 * n + 1, self.ruecklauf_consumptionnode[n]] = -1
        self.vorlauf_consumptionnode_unique, idx = np.unique(self.vorlauf_consumptionnode, return_index=True)
        self.ruecklauf_consumptionnode_unique = self.ruecklauf_consumptionnode[idx]
        self.nbOutletsUnique = idx.shape[0]
        self.M_w = np.zeros((self.nbOutletsUnique, self.nbNodes))
        for n in range(self.nbOutletsUnique):
            self.M_w[n, self.vorlauf_consumptionnode_unique[n]] = 1
            self.M_w[n, self.ruecklauf_consumptionnode_unique[n]] = 1
            
        self.Outlets_Tref=self.Outlets[:,9].astype(float)
        self.Outlets_unique_Tref=self.Outlets_Tref[idx]


    def data_to_consumption(self):
        """Build consumption data."""
        self.vorlauf_consumptionnode = self.Outlets[:, 1].astype(int)
        self.ruecklauf_consumptionnode = self.Outlets[:, 2].astype(int)

    def data_to_valve_matrix(self):
        """Build valve matrix V from data."""
        self.ValvesPipes = np.concatenate((self.Valves[:, 1], self.Sliders[:, 1]))

    def data_to_supply_matrix(self):
        """Build supply matrix E from data."""
        self.vorlauf_SupplieNodes = self.Supplies[:, 1].astype(int)
        self.ruecklauf_SupplieNodes = self.Supplies[:, 5].astype(int)
        self.E = np.zeros((2 * self.nbSupplies, self.nbNodes))
        self.M_S = np.zeros((self.nbSupplies, self.nbNodes))
        for n in range(self.nbSupplies):
            self.E[2 * n, self.vorlauf_SupplieNodes[n]] = 1
            self.E[2 * n + 1, self.ruecklauf_SupplieNodes[n]] = -1
            self.M_S[n,self.vorlauf_SupplieNodes[n]] = 1
            self.M_S[n,self.ruecklauf_SupplieNodes[n]] = -1
            
        

    def data_to_null_node_matrix(self):
        """Build null node matrix from data."""
        temp = np.arange(self.nbNodes)
        temp = np.setdiff1d(temp, self.vorlauf_consumptionnode)
        temp = np.setdiff1d(temp, self.ruecklauf_consumptionnode)
        temp = np.setdiff1d(temp, self.vorlauf_SupplieNodes)
        temp = np.setdiff1d(temp, self.ruecklauf_SupplieNodes)
        temp = np.setdiff1d(temp, self.Pumps_startnode)
        temp = np.setdiff1d(temp, self.Pumps_targetnode)
        self.NullNodes = temp

    def data_to_meas_network_matrices(self):
        """Build measurement network matrices from data."""
        self.P=np.zeros((self.nbpMeas,self.nbNodes))
        self.F=np.zeros((self.nbmMeas,self.nbPipes))
        for n in range(self.nbpMeas):
            self.P[n,self.pMeasurements[n,1]]=1
        for n in range(self.nbmMeas):
            self.F[n,self.mMeasurements[n,1]]=1

    def pipes_data(self):
        """Process pipes data."""
        self.Pipes_d=self.Pipes[:,3] #inner diameter
        self.Pipes_l=self.Pipes[:,2] #length
        self.Pipes_k=self.Pipes[:,4] #roughness
        self.Pipes_A=math.pi/4*np.power(self.Pipes_d,2) # inner Area
        self.Pipes_cr=8/np.power(math.pi,2)*self.Pipes_l/np.power(self.Pipes_d,5) #geometrischer Beiwert
        self.Pipes_h=np.matmul(-self.M.transpose(),self.Nodes[:,1]).astype(float) #heigth difference

        #self.Pipes_next=self.Pipes[:,5].astype(int)

    def outlets_data(self):
        """Process outlets data."""
        self.Outlets_ConsumerType=self.Outlets[:,4]
        self.Outlets_ConsumerType[self.Outlets_ConsumerType=='MFH_O_RW']='MultiFamilyHouse.RoomAndWaterOld'
        self.Outlets_ConsumerType[self.Outlets_ConsumerType=='OFH_N_RW']='OneFamilyHouse.RoomAndWaterNew'
        self.Outlets_ConsumerType[self.Outlets_ConsumerType=='OFH_O_RW']='OneFamilyHouse.RoomAndWaterOld'
        self.Outlets_ConsumerType[self.Outlets_ConsumerType=='WHS']='Wholesale'
        self.Outlets_ConsumerType[self.Outlets_ConsumerType=='CRI']='CreditInstitutes'
        self.Outlets_ConsumerType[self.Outlets_ConsumerType=='LDG']='Lodging'
        
        self.Outlets_Qa=self.Outlets[:,5].astype(float)
        self.Outlets_Pmax=self.Outlets[:,6].astype(float)
        
    def pump_data(self):
        """Process Supply data."""
        self.Pumps_m0=self.Pumps[:,8].astype(float)
        self.Pumps_p0=self.Pumps[:,9].astype(float)
        self.Pumps_n0=self.Pumps[:,10].astype(float)
        self.Pumps_eta0=self.Pumps[:,11].astype(float)


    def valve_state(self):
        """Process valve state data."""
        self.V=np.zeros((self.nbValves,self.nbPipes))
        delete=[]
        self.closedPipes=[]
        self.closedValves=self.nbValves
        closed=np.concatenate((self.ValveStates[:,2], self.SliderStates[:,1]))
        for n in range(self.nbValves):
            if closed[n]==0  :
                delete.append(n)
                self.closedValves=self.closedValves-1
            elif closed[n]==1:
                self.V[n,self.ValvesPipes[n]]=1
                self.closedPipes.append(self.ValvesPipes[n])
            else:
                pass
        self.V=np.delete(self.V, delete,axis=0)
    
    def DNdims(self):

        
        self.pipe_dimensions = {
            "DN0020": {"d_st_in": 0.0217, "d_st_out": 0.0269, "d_ins_out": 0.084, "d_cas_out": 0.09},
            "DN0025": {"d_st_in": 0.0285, "d_st_out": 0.0337, "d_ins_out": 0.084, "d_cas_out": 0.09},
            "DN0032": {"d_st_in": 0.0372, "d_st_out": 0.0424, "d_ins_out": 0.104, "d_cas_out": 0.11},
            "DN0040": {"d_st_in": 0.0431, "d_st_out": 0.0483, "d_ins_out": 0.104, "d_cas_out": 0.11},
            "DN0050": {"d_st_in": 0.0545, "d_st_out": 0.0603, "d_ins_out": 0.119, "d_cas_out": 0.125},
            "DN0065": {"d_st_in": 0.0703, "d_st_out": 0.0761, "d_ins_out": 0.134, "d_cas_out": 0.14},
            "DN0080": {"d_st_in": 0.0825, "d_st_out": 0.0889, "d_ins_out": 0.154, "d_cas_out": 0.16},
            "DN0100": {"d_st_in": 0.1071, "d_st_out": 0.1143, "d_ins_out": 0.1936, "d_cas_out": 0.2},
            "DN0125": {"d_st_in": 0.1325, "d_st_out": 0.1397, "d_ins_out": 0.2182, "d_cas_out": 0.225},
            "DN0150": {"d_st_in": 0.1603, "d_st_out": 0.1683, "d_ins_out": 0.2428, "d_cas_out": 0.25},
            "DN0200": {"d_st_in": 0.2101, "d_st_out": 0.2191, "d_ins_out": 0.3068, "d_cas_out": 0.315},
            "DN0250": {"d_st_in": 0.263, "d_st_out": 0.273, "d_ins_out": 0.3904, "d_cas_out": 0.4},
            "DN0300": {"d_st_in": 0.3127, "d_st_out": 0.3239, "d_ins_out": 0.4396, "d_cas_out": 0.45},
            "DN0350": {"d_st_in": 0.3444, "d_st_out": 0.3556, "d_ins_out": 0.4888, "d_cas_out": 0.5},
            "DN0400": {"d_st_in": 0.3938, "d_st_out": 0.4064, "d_ins_out": 0.5486, "d_cas_out": 0.56},
            "DN0450": {"d_st_in": 0.4444, "d_st_out": 0.457, "d_ins_out": 0.618, "d_cas_out": 0.63},
            "DN0500": {"d_st_in": 0.4954, "d_st_out": 0.508, "d_ins_out": 0.6968, "d_cas_out": 0.71},
            "DN0600": {"d_st_in": 0.5958, "d_st_out": 0.61, "d_ins_out": 0.7844, "d_cas_out": 0.8},
            "DN0700": {"d_st_in": 0.695, "d_st_out": 0.711, "d_ins_out": 0.8826, "d_cas_out": 0.9},
            "DN0800": {"d_st_in": 0.7954, "d_st_out": 0.813, "d_ins_out": 0.9812, "d_cas_out": 1.0},
            "DN0900": {"d_st_in": 0.894, "d_st_out": 0.914, "d_ins_out": 1.0796, "d_cas_out": 1.1},
            "DN1000": {"d_st_in": 0.994, "d_st_out": 1.016, "d_ins_out": 1.178, "d_cas_out": 1.2},
            "DN1100": {"d_st_in": 1.096, "d_st_out": 1.118, "d_ins_out": 1.2764, "d_cas_out": 1.3},
            "DN1200": {"d_st_in": 1.194, "d_st_out": 1.219, "d_ins_out": 1.375, "d_cas_out": 1.4}
        }



        

    
    def findDNpipe(self):
        self.Pipes_DN = np.empty(self.nbPipes, dtype=object)              
        for i in range(len(self.Pipes_d)):
            if self.Pipes_d[i] == self.pipe_dimensions["DN0020"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0020"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0025"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0025"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0032"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0032"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0040"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0040"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0050"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0050"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0065"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0065"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0080"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0080"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0100"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0100"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0125"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0125"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0150"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0150"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0200"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0200"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0250"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0250"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0300"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0300"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0350"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0350"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0400"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0400"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0450"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0450"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0500"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0500"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0600"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0600"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0700"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0700"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0800"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0800"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN0900"]["d_st_in"]:
                self.Pipes_DN[i] = "DN0900"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN1000"]["d_st_in"]:
                self.Pipes_DN[i] = "DN1000"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN1100"]["d_st_in"]:
                self.Pipes_DN[i] = "DN1100"
            elif self.Pipes_d[i] == self.pipe_dimensions["DN1200"]["d_st_in"]:
                self.Pipes_DN[i] = "DN1200"

            else:
                #print("Pipe"+ str(i)+" has no DN fitting DN Series")
                closest_dn = None
                min_diff = float('inf')
                for dn, dim in self.pipe_dimensions.items():
                    diff = abs(dim["d_st_in"] - self.Pipes_d[i])
                    if diff < min_diff:
                        min_diff = diff
                        closest_dn = dn
                #print(min_diff)
                self.Pipes_DN[i] = closest_dn
                        
                        
                
            
        
