# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 15:40:29 2024

@author: sander
"""
import numpy as np
import math
from pyXSteam.XSteam import XSteam
import datetime
from StandardConsumptionProfile import StandardConsumptionProfile
from DynamicTemperature import DynamicTemperature
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve,inv, lsqr
from scipy.integrate import quad
import random
import copy



class TimeSeries:
    """Class for handling a TimeSeries Calculation into the future"""
    
    def __init__(self,StaticData,hydraulicMesh,ProgramModes):
        """Initalizes all necassary Variables for the calculation"""
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.T=np.ones(StaticData.nbNodes)*80
        self.T_old=np.ones(StaticData.nbNodes)
        for n in range(StaticData.nbNodes):
            if StaticData.Nodes[n,2]==1:
                self.T[n]=110
            if StaticData.Nodes[n,2]==-1:
                self.T[n]=60
        self.p=np.ones(StaticData.nbNodes)*10
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
    
    # Method to store the actual stat
    def save_state(self):
        return copy.deepcopy(self.__dict__)
    
    # Method to conda install y a saved state
    def restore_state(self, saved_state):
        self.__dict__ = copy.deepcopy(saved_state)


        
        
    def calcThermalResistanceVertical(self,StaticData):
        
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
        """One timestep for a completly stationary Calculation """
        if self.datetimeStationary[-1]!=ProgramModes.startTime:
            self.datetimeStationary.append(self.datetimeStationary[-1]+self.time_change_stationary)

        self.j=0
        self.Pump_inactive=[]
        self.badPoint_violated=[]
        self.T_old=np.ones(StaticData.nbNodes)
        if not(ProgramModes.PowerInput=='None'):
            self.Qstrom=LoadMeasValues.SupplyVolumeflow(self.datetimeStationary[-1])
            self.T_in_meas=LoadMeasValues.SupplyTemperature(self.datetimeStationary[-1])
            self.updateQ(StaticData,ProgramModes)
        else:
            self.updateQ(StaticData,ProgramModes)
        
        while np.linalg.norm(self.T-self.T_old,2)>StaticData.nbNodes*ProgramModes.accuracy_factor_thermic and self.j<30:
            self.T_old=self.T
            
            self.waterProperties(StaticData)
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
        
        """ One Timestep for a quasi dynamic calculation """
        
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
        """ One dynamic thermal timestep"""
        if self.DynamicTemperature.step>0:
            self.T=self.DynamicTemperature.Timestep(StaticData,self,ProgramModes,Tsupply)
            self.datetimeDynamic.append(self.datetimeDynamic[-1]+self.time_change_dynamic)
        else:
            self.DynamicTemperature.Timestep(StaticData,self,ProgramModes,Tsupply)
            
    
    def updateQ(self,StaticData,ProgramModes):
        """ Updates the needed heating power for outlets regarding Standardlastprofil"""
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
        """Calculates the massflow at the Outlets depending on the power and the tempreture difference"""
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
        """" sets the outer massflows vector according to the consumer massflows and the free (n-1) Suppliers""" 
        self.L=np.zeros((StaticData.nbNodes))
        for n in range(StaticData.nbOutlets):
            self.L[StaticData.Outlets[n,1]]=self.L[StaticData.Outlets[n,1]]+self.m_consumer[n]
            self.L[StaticData.Outlets[n,2]]=self.L[StaticData.Outlets[n,2]]-self.m_consumer[n]
            
            
            # self.L[StaticData.vorlauf_SupplieNodes[1:]]=-70
            # self.L[StaticData.ruecklauf_SupplieNodes[1:]]=70
        
    def calcm0(self,hydraulicMesh,StaticData):
        """Calculates the massflows zero through the Forest of the network (without meshes)"""
        A_0_LE=np.zeros((len(hydraulicMesh.leftEdges),hydraulicMesh.A_Forest.shape[1]))
        for n in range(len(hydraulicMesh.leftEdges)):
            A_0_LE[n,hydraulicMesh.leftEdges[n]]=1
            
        A_0_Pipes=np.zeros((len(StaticData.closedPipes),hydraulicMesh.A_Forest.shape[1]))
        for n in range(len(StaticData.closedPipes)):
            A_0_Pipes[n,StaticData.closedPipes[n]]=1

            # A_0_Test=np.zeros((1,hydraulicMesh.A_Forest.shape[1]))
            # A_0_Test[0,1]=1
           #A_0_Test[1,1742]=1
        
                
        A=np.concatenate((hydraulicMesh.A_Forest,A_0_LE,A_0_Pipes))
        L=np.concatenate((self.L,np.zeros(len(hydraulicMesh.leftEdges)+len(StaticData.closedPipes))))
        #L=L[~np.all(A == 0, axis=1)]
        #A = A[~np.all(A == 0, axis=1)]
        idx=np.array([StaticData.vorlauf_SupplieNodes[0],StaticData.ruecklauf_SupplieNodes[0]])
                                                                                      
        L = np.delete(L,idx,axis=0)
        A = np.delete(A,idx,axis=0)
        
        L = L[~np.all(A == 0, axis=1)]
        A = A[~np.all(A == 0, axis=1)]
        
            
        #A = A[:,~np.all(A == 0, axis=0)]
        self.m0=np.linalg.solve(A,L)
        
        
    def calcDeltaT(self,StaticData):
        """Calculactes the Tempreture Difference at the Consumers accoording to the system states"""
        self.Delta_T_consumer=self.T[StaticData.vorlauf_consumptionnode]-StaticData.Outlets_Tref
        
        # if np.any(self.Delta_T_consumer<=10):
        #     idx=np.where(self.Delta_T_consumer<=10)[0]
        #     print(str(idx))
        #     #self.Delta_T_consumer[idx]=10
        
    def calcDeltaTdyn(self,StaticData):
        """Calculactes the Tempreture Difference at the Consumers accoording to the system states"""
        self.Delta_T_consumer=self.T[StaticData.vorlauf_consumptionnode]-StaticData.Outlets_Tref
        

        
    def hydraulicPrep(self,hydraulicMesh,StaticData,ProgramModes):
        """Does a hydraulic Preparation to perform the next Newton Step"""
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
        """Calculate the function F (sum of pressure differences through all meshes) which should be zero at the end"""
        self.F=hydraulicMesh.B[:,:StaticData.nbPipes]@(self.R*self.m_i)
        
    def calcJ(self,hydraulicMesh,StaticData):
        """approximate the derivative of F, excluding Pipe Friction"""
        self.J=hydraulicMesh.B[:,:StaticData.nbPipes]@np.diag(2*self.R)@hydraulicMesh.B[:,:StaticData.nbPipes].transpose()
        
    def NewtonRaphsonStep(self):
        """ Perform one NewtonRaphson Step"""
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
        """Other Method for NewtonRaphson Step"""
        self.m_s_old=self.m_s
        A=hydraulicMesh.B[:,:StaticData.nbPipes]@np.diag(self.R)@hydraulicMesh.B[:,:StaticData.nbPipes].transpose()
        B=-hydraulicMesh.B[:,:StaticData.nbPipes]@np.diag(self.R*2)@self.m0[:StaticData.nbPipes]
        self.m_s=np.linalg.solve(A,B)
        
    
        
    def pressureCond(self,StaticData):
        """Setting conditions for pressure calculations into matrices"""
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
        """calculating the pressure in the node"""
        while True:
            try:
                self.A_T
            except:
                self.A_T=hydraulicMesh.A_Forest.transpose()*(-1)
                self.A_T=self.A_T[:StaticData.nbPipes,:]
                self.A_T = self.A_T[~np.all(self.A_T == 0, axis=1)]
            P=np.concatenate((self.A_T,self.BP,self.FP,StaticData.M_S))
            delta_p=np.delete(self.R*self.m_i-self.Pipes_pg,np.concatenate((hydraulicMesh.leftEdges,StaticData.closedPipes)).astype(int))
            
            P_X=np.concatenate((delta_p,self.BPdelta,self.FPValue,np.zeros(StaticData.nbSupplies))).astype(float)#.reshape(P.shape[0],1)
            #self.p=np.linalg.solve(P,P_X)
            # Convert the dense matrices to sparse CSC format
            P_sparse = csc_matrix(P)
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
        """preparing stationary thermic caclvulation"""
        self.buildSaSe(StaticData)
        self.buildTemperatureMatrix(StaticData,ProgramModes,Tsupply)
        self.buildsetT(StaticData,Tsupply)
        
            
        
    
    def resistanceMatrix(self,StaticData,ProgramModes):
        """setting new Resistance entries according to states at the current iteration"""
        self.waterProperties(StaticData) #update Waterproteries maybe change?
        self.v= self.m_i/self.rho/StaticData.Pipes_A #flow velocity v=m/rho/A
        self.Re= self.rho*abs(self.v)*StaticData.Pipes_d/self.viscosity #Reynoldszahl rho/*|v|*d/viscosity
        self.calcLambda(StaticData,ProgramModes.PipeFrictionMethod) #calculate lambda colebrook
        self.R=(self.Pipes_lambda*StaticData.Pipes_cr*np.clip(abs(self.m_i),1,None)*1e-5)/self.rho #calculate  R vector R=lambda*cr/rho*|mi|*1e-5 convert to bar
        
    def waterProperties(self,StaticData):
        """calculating the water Properties (viscosity, density & spezific heat capacity) according to pressure and tempreture in the iteration"""
        # self.viscosity_Nodes=np.zeros(StaticData.nbNodes)
        # self.rho_Nodes=np.zeros(StaticData.nbNodes)
        # self.cp_Nodes=np.zeros(StaticData.nbNodes)
        for n in range(StaticData.nbNodes):
            self.viscosity_Nodes[n]=self.steamTable.my_pt(self.p[n],self.T[n])
            self.rho_Nodes[n]=self.steamTable.rho_pt(self.p[n],self.T[n])
            #self.cp_Nodes[n]=self.steamTable.Cp_pt(self.p[n],self.T[n])*1000
            

        
        for n in range(StaticData.nbPipes):
            self.viscosity[n]=1/2*(self.viscosity_Nodes[StaticData.startnode[n]]+self.viscosity_Nodes[StaticData.targetnode[n]])
            self.rho[n]=1/2*(self.rho_Nodes[StaticData.startnode[n]]+self.rho_Nodes[StaticData.targetnode[n]])
            self.cp[n]=1/2*(self.cp_Nodes[StaticData.startnode[n]]+self.cp_Nodes[StaticData.targetnode[n]])

            

            
        for n in range(StaticData.nbPumps):
            
            self.viscosity_Pumps[n]=1/2*(self.viscosity_Nodes[StaticData.Pumps_startnode[n]]+self.viscosity_Nodes[StaticData.Pumps_startnode[n]])
            self.rho_Pumps[n]=1/2*(self.rho_Nodes[StaticData.Pumps_startnode[n]]+self.rho_Nodes[StaticData.Pumps_startnode[n]])
            self.cp_Pumps[n]=1/2*(self.cp_Nodes[StaticData.Pumps_startnode[n]]+self.cp_Nodes[StaticData.Pumps_startnode[n]])
            
        for n in range(StaticData.nbOutlets):
            self.cp_consumer[n]=1/2*(self.cp_Nodes[StaticData.vorlauf_consumptionnode[n]]+self.cp_Nodes[StaticData.ruecklauf_consumptionnode[n]])
        for n in range(StaticData.nbSupplies):
            self.cp_Supply[n]=1/2*(self.cp_Nodes[StaticData.vorlauf_SupplieNodes[n]]+self.cp_Nodes[StaticData.ruecklauf_SupplieNodes[n]])
            self.rho_Supply[n]=1/2*(self.rho_Nodes[StaticData.vorlauf_SupplieNodes[n]]+self.rho_Nodes[StaticData.ruecklauf_SupplieNodes[n]])
   

        
            
    def calcLambda(self,StaticData,PipeFrictionMethod):
        """calculate the Pipe friction coefficient according to the Reynold-Number etc.
        Different Formulas for the turbulent calculation implemented"""
        idx_laminar= np.where(self.Re<=2320) #pipe idx where flow is laminar
        idx_turbulent= np.where(self.Re>2320) #pipe idx where flow is turbulent 
        with np.errstate(divide='ignore', invalid='ignore'):
            self.Pipes_lambda[idx_laminar]=np.clip(64/self.Re[idx_laminar],None,1) #for laminar flow is lambda=64/Re   #maxed at lambda =1 for slow massflows
        
        lambda_new=self.Pipes_lambda[idx_turbulent] #set lamdda new
        Re=self.Re[idx_turbulent] 
        d=StaticData.Pipes_d[idx_turbulent]
        k=StaticData.Pipes_k[idx_turbulent]
        
        if len(idx_turbulent[0])>0:            
            match PipeFrictionMethod:
                
                
                case 'ColebrookWhite':
                    lambda_pre=0.01+lambda_new #set lamda pre different to lamda_new fo start
                    for n in range(len(idx_turbulent[0])):                    
                        lambda_new[n]=math.pow(1.8*math.log10(6.9/Re[n]+math.pow(k[n]/3.7/d[n], 1.11)), -2)#Haaland als Initial Guees
                      
                        #hydraulic even
                        if Re[n]<math.pow(k[n]/d[n],8/7):  
                            while abs(lambda_new[n]-lambda_pre[n])>0.000001:
                                  lambda_pre[n]=lambda_new[n]
                                  lambda_new[n]=math.pow((-2.035*math.log10(2.503/Re[n]/math.sqrt(lambda_pre[n]))),-2)
                                  
                          #normal colebrook white
                        elif Re[n]>math.pow(k[n]/d[n],8/7) and Re[n]<= -40000/k[n]*d[n]*math.log10(k[n]/d[n]/3.729): 
                            while abs(lambda_new[n]-lambda_pre[n])>0.000001:
                                  lambda_pre[n]=lambda_new[n]
                                  lambda_new[n]=math.pow((-2.035*math.log10(2.503/Re[n]/math.sqrt(lambda_pre[n])+k[n]/3.729/d[n])),-2)
                              
                          #hydraulic rough
                        else:
                              lambda_new[n]=math.pow(-2.035*math.log10(k[n]/d[n]/3.729),-2)
                              
                        self.Pipes_lambda[idx_turbulent[0]]=lambda_new


                case 'Alashkar':
                    A= k/d/3.7065
                    B= 2.5226/Re
                    self.Pipes_lambda[idx_turbulent[0]]= 1.325474505 * np.power(np.log(A-0.8686068432*B*np.log(A-0.8784893582*B*np.log(A+np.power(1.665368035*B, 0.8373492157)))), -2)
            
            
                case 'Haaland':
                    self.Pipes_lambda[idx_turbulent[0]]=np.power(1.8*np.log10(6.9/Re+np.power(k/3.7/d, 1.11)), -2)
            
                case 'GoudarSonnad':
                    A=2/np.log(10)
                    B=k/d/3.7
                    D=np.log(10)*Re/5.02
                    S=B*D+np.log(D)
                    Q=np.power(S,S/(S+1))
                    G= B*D+np.log(D/Q)
                    Z=np.log(Q/G)
                    DLA=Z*G/(G+1)
                    DCFA=DLA*(1+Z/2/(np.power(G+1,2)+Z/3*(2*G-1)))
                    self.Pipes_lambda[idx_turbulent[0]]=np.power(A*np.log(D/Q)+DCFA,-2)
                    
            
                case 'Tkachenko':
                    A0=-0.79638*np.log(k/d/8.208+7.3357/Re)
                    A1=Re*k/d+9.3120665*A0
                    self.Pipes_lambda[idx_turbulent[0]]=np.power((8.128943+A1)/(8.128943*A0-0.86859209*A1*np.log(A1/3.7099535/Re)),2)
                
                case 'Avci':
                    A= 5000*np.power(k/d,3)
                    B=10*np.sqrt(k/d)/(1+225*np.power(k/d,2))
                    self.Pipes_lambda[idx_turbulent[0]]=6.4/np.power(np.log(1/Re+0.01*k/d*(1+A+B)),2.4)
                    
                    
                
                    
        else:    
            pass
        
    def buildSaSe(self,StaticData):
        """build massflow direction dependend matrices input and output for each edge for themic calculation"""
        start=np.concatenate((StaticData.startnode,StaticData.Pumps_startnode,StaticData.vorlauf_consumptionnode, StaticData.ruecklauf_SupplieNodes)) #startnodes Pipes and Pumps ,Outlets and Supplie
        target=np.concatenate((StaticData.targetnode,StaticData.Pumps_targetnode, StaticData.ruecklauf_consumptionnode,StaticData.vorlauf_SupplieNodes)) #targetnodes Pipes and Pumps,,Outlets and Supplie
        idx=np.where(self.m_i<0) #idx with negative massflow --> temperature change in the other direction
        if len(idx[0])>0:
            start[idx],target[idx]=target[idx],start[idx] # switch Massflow direction
        self.Se=np.zeros((start.shape[0],StaticData.nbNodes)) #Kanten-Eintritts-Matrix
        self.Sa=np.zeros((target.shape[0],StaticData.nbNodes)) #Kanten-Austritts-Matrix 
        for n in range(start.shape[0]): # iterarte through start/target maybe check if same length else error
            self.Se[n,start[n]]=1 #set Se entries
            self.Sa[n,target[n]]=1 #set Sa entries
        idxdel=np.where(abs(self.m_i)<=1e-6)[0]
        self.Sa[idxdel,:]=0
        self.Se[idxdel,:]=0
            
    def buildTemperatureMatrix(self,StaticData,ProgramModes,Tsupply):
        """build the heatflow matrices etc. for thermic calculation"""
        self.W=abs(np.concatenate((self.m_i,self.m_p,self.m_consumer,np.atleast_1d(self.m_supply)),axis=0))*np.concatenate((self.cp, self.cp_Pumps,self.cp_consumer,np.atleast_1d(self.cp_Supply)),axis=0) #vektor (mi*cp)
        if ProgramModes.HeatLossMode in [ 'default','simple']:
            self.T_pipe_mean=np.zeros(StaticData.nbPipes)
            for n in range(StaticData.nbPipes):
                self.T_pipe_mean[n]=1/2*(self.T[StaticData.startnode[n]]+self.T[StaticData.targetnode[n]])
            self.Q_loss_Pipe=StaticData.Pipes_l/self.Ur*(self.T_pipe_mean-self.T_out)
            tempexp=np.exp(-StaticData.Pipes_l/\
                           (StaticData.Pipes_A*self.cp*self.rho*self.Ur*np.clip(abs(self.v),0.001,np.inf))).astype(float) # e^(-Ur*l/d/cp/rho/v)
            
            C_Outlets=np.ones(StaticData.nbOutlets)
            # Predefine D_Outlets as zeros
            D_Outlets = np.zeros_like(self.Q)
            
            # Perform the calculation only where m_consumer is non-zero
            non_zero_mask = self.m_consumer != 0
            D_Outlets[non_zero_mask] = -self.Q[non_zero_mask] / (self.m_consumer[non_zero_mask] * self.cp_consumer[non_zero_mask])
            
            C_Supply=np.zeros(StaticData.nbSupplies)
            D_Supply=Tsupply
                
            self.C=np.diag(np.concatenate((tempexp,np.ones(StaticData.nbPumps),C_Outlets,C_Supply),axis=0)) #C matrix Diagonal tempexp für Leitungen und 1 für Pumpen
            
            self.D=np.concatenate((self.T_out*(1-tempexp),np.zeros(StaticData.nbPumps),D_Outlets,D_Supply),axis=0) #D Vektor T
        
        if ProgramModes.HeatLossMode not in [ 'default','simple']:
            self.T_pipe_mean=np.zeros(StaticData.nbPipes)
            for n in range(StaticData.nbPipes):
                self.T_pipe_mean[n]=1/2*(self.T[StaticData.startnode[n]]+self.T[StaticData.targetnode[n]])
            self.Q_loss=np.zeros(StaticData.nbPipes)
            
            if ProgramModes.HeatLossNextPipe in ['No']:
                T_supply_avg=np.average(self.T[StaticData.supplyNodes])
                T_return_avg=np.average(self.T[StaticData.returnNodes])
            if ProgramModes.HeatLossNextPipe in ['Yes']:
                T_next=self.T_pipe_mean[StaticData.Pipes_next]
                
                
            if ProgramModes.HeatLossMode  in ['EMCM']:                    
                
                if ProgramModes.HeatLossNextPipe in ['No']:
                #calculate Heat Loss in W/m per Pipe avg Tempreture of other Side next to it with EMCM(Equivalent Mesh current method)
                    for n in range(StaticData.nbPipes):
                        
                        if n in StaticData.supply_pipes:
                            self.Q_loss[n]=self.heatLossEMCM(self.T_pipe_mean[n],T_return_avg,self.T_out,StaticData.pipe_dimensions[StaticData.Pipes_DN[n]])
                        if n in StaticData.return_pipes:
                            self.Q_loss[n]=self.heatLossEMCM(self.T_pipe_mean[n],T_supply_avg,self.T_out,StaticData.pipe_dimensions[StaticData.Pipes_DN[n]])
                        
                if ProgramModes.HeatLossNextPipe in ['Yes']:
                    for n in range(StaticData.nbPipes):
                        self.Q_loss[n]=self.heatLossEMCM(self.T_pipe_mean[n],T_next[n],self.T_out,StaticData.pipe_dimensions[StaticData.Pipes_DN[n]])
                        
                 
            if ProgramModes.HeatLossMode  in ['Norm']:
                if ProgramModes.HeatLossNextPipe in ['No']:
                    #calculate Heat Loss in W/m per Pipe avg Tempreture of other Side next to it with Norm
                    for n in range(StaticData.nbPipes):
                        if n in StaticData.supply_pipes:
                            self.Q_loss[n]=self.heatLossNorm(self.T_pipe_mean[n],T_return_avg,self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[n]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[n]]["d_st_out"]) 
                        if n in StaticData.return_pipes:
                            self.Q_loss[n]=self.heatLossNorm(self.T_pipe_mean[n],T_supply_avg,self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[n]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[n]]["d_st_out"]) 
                
                if ProgramModes.HeatLossNextPipe in ['Yes']:
                    for n in range(StaticData.nbPipes):
                        self.Q_loss[n]=self.heatLossNorm(self.T_pipe_mean[n],T_next[n],self.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[n]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[n]]["d_st_out"])
            self.Q_loss_Pipe=self.Q_loss*StaticData.Pipes_l
            
            C_Outlets=np.ones(StaticData.nbOutlets)
            # Predefine D_Outlets as zeros
            D_Outlets = np.zeros_like(self.Q)
            
            # Perform the calculation only where m_consumer is non-zero
            non_zero_mask = self.m_consumer != 0
            D_Outlets[non_zero_mask] = -self.Q[non_zero_mask] / (self.m_consumer[non_zero_mask] * self.cp_consumer[non_zero_mask])
            
            C_Supply=np.zeros(StaticData.nbSupplies)
            D_Supply=Tsupply
    
            self.C=np.diag(np.concatenate((np.ones(StaticData.nbPipes+StaticData.nbPumps),C_Outlets,C_Supply),axis=0))
            #DPipes=-self.Q_loss*StaticData.Pipes_l/(self.cp*abs(np.clip(self.m_i,0.1,np.inf)))
            DPipes=-self.Q_loss*StaticData.Pipes_l/(StaticData.Pipes_A*self.cp*self.rho*np.clip(abs(self.v),0.001,np.inf))
            self.D=np.concatenate((DPipes,np.zeros(StaticData.nbPumps),D_Outlets,D_Supply),axis=0)
        # Convert dense matrices to sparse
        Sa_sparse = csr_matrix(self.Sa)
        Se_sparse = csr_matrix(self.Se)
        C_sparse = csr_matrix(self.C)
        D_sparse = csr_matrix(self.D.astype(float))
        
        # Efficient diagonal multiplication
        W_diag = csr_matrix(np.diag(self.W.astype(float)))  # Alternatively, handle W without constructing a diagonal matrix
        
        # Compute T_A and T_B using sparse operations
        T_A = Sa_sparse.T @ W_diag @ C_sparse @ Se_sparse - Se_sparse.T @ W_diag @ Se_sparse
        T_B = -Sa_sparse.T @ W_diag @ D_sparse.T
        
        # Ensure T_A and T_B are sparse matrices
        if not isinstance(T_A, csr_matrix):
            T_A = csr_matrix(T_A)
        if not isinstance(T_B, csr_matrix):
            T_B = csr_matrix(T_B)
        
        # Create mask for rows with non-zero entries
        mask = T_A.getnnz(axis=1) > 0
        

        
        # Apply mask
        T_A = T_A[mask]
        T_B = T_B[mask]
        
        # Convert sparse matrices to dense arrays
        self.T_A = T_A.toarray() if isinstance(T_A, csr_matrix) else np.array(T_A)
        self.T_B = T_B.toarray() if isinstance(T_B, csr_matrix) else np.array(T_B)
        
                
    def buildsetT(self,StaticData,Tsupply):
        """set the T values at the input of the networks (supply in the supply network and consumer at the return network) """
        setTsupply=np.zeros((StaticData.nbSupplies,StaticData.nbNodes))
        # for n in range(StaticData.nbSupplies):
        #     setTsupply[n,StaticData.vorlauf_SupplieNodes[n]]=1
            
        # setTconsumer=np.zeros((StaticData.nbOutletsUnique,StaticData.nbNodes))
        # for n in range(StaticData.nbOutletsUnique):
        #     setTconsumer[n,StaticData.ruecklauf_consumptionnode_unique[n]]=1
        # self.Tset=np.concatenate((setTsupply,setTconsumer))
        
        # Tconsumer=StaticData.Outlets_unique_Tref
        
        # self.setTValues=np.concatenate((Tsupply,Tconsumer))
        
    def calcT(self,StaticData):
        """calc stationary T  """
        # T=np.concatenate((self.T_A,self.Tset))
        # T_X=np.concatenate((self.T_B,self.setTValues))
        T=np.array(self.T_A).astype(float)
        T_X=np.array(self.T_B).astype(float).flatten()
        idx=np.where(sum(abs(T))==0)[0]
        T0=np.zeros((len(idx),StaticData.nbNodes))
        for n in range(len(idx)):
            T0[n,idx[n]]=1
        T=np.concatenate((T,T0))
        T_X=np.concatenate((T_X,30*np.ones(len(idx))))
        
        #self.T=np.linalg.inv(T.transpose()@T)@T.transpose()@T_X
        
        norms = np.linalg.norm(T, axis=1, keepdims=True)
        T_normalized = T / norms
        
        # Convert dense matrices to sparse matrices
        T_sparse = csc_matrix(T_normalized)
        T_X_normalized= T_X/ norms.squeeze()
        
        # # Compute (T^T * T)
        # TtT_sparse = T_sparse.transpose() @ T_sparse
        
        # # Solve (T^T * T) * X = (T^T * T_X) for X
        # TtT_T_X_sparse = T_sparse.transpose() @ T_X_normalized.transpose()
        # X_sparse = spsolve(TtT_sparse, TtT_T_X_sparse)
        
        # # Assign to self.T
        # self.T = X_sparse
        self.T = lsqr(T_sparse,T_X_normalized,x0=self.T,atol=1e-12,btol=1e-12)[0]
        #self.T = lsqr(T,T_X,x0=self.T,atol=1e-12,btol=1e-12)[0]
        
        
        
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
        self.HeatLoses=sum(self.Q_loss_Pipe)
        #self.HeatLosesTest=sum(self.Q_loss*StaticData.Pipes_l)
        #self.HeatLosesTest2=sum(abs(self.m_i*(self.T[StaticData.startnode]*self.cp_Nodes[StaticData.startnode]-self.T[StaticData.targetnode]*self.cp_Nodes[StaticData.targetnode])))
        #self.HeatLosesTest2=sum(self.m_i*(self.T[StaticData.startnode]-self.T[StaticData.targetnode])*self.cp)

