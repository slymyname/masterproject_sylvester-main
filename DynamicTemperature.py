# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 13:57:18 2023

@author: sander
"""

from scipy.interpolate import BSpline, make_interp_spline, lagrange
import numpy as np
import matplotlib as mpl

class DynamicTemperature:
    
    def __init__(self,ProgramModes):
        self.ti=ProgramModes.dynamicTimestep
        self.int_faster=5
        self.T_dynamic=[]
        self.step=int(0)
        self.flowtime=0
        self.v=[]
        self.delta_x_max=ProgramModes.delta_x_max
        self.TPipe=[]
        
    def Timestep(self,StaticData,TimeSeries,ProgramModes,Tsupply):
        """ Dynamic TempretureTimeStep for different Methods""" 
        if ProgramModes.TemperatureDynamicMode in [ 'Icking']:
            self.v.append(TimeSeries.v)
            self.calculatetime(StaticData,TimeSeries)
            self.buildSE(StaticData,TimeSeries)
            self.calculateTE(StaticData,TimeSeries)
            self.calculateTcoeff(StaticData,TimeSeries)
            self.calcTk(StaticData,TimeSeries)
            self.setTValues(StaticData,TimeSeries,Tsupply)
            T=self.finishstep()
            return T
        elif ProgramModes.TemperatureDynamicMode in ['default','ImplicitUpwind']:
            self.spanControlVolumes(StaticData)
            self.buildSE(StaticData,TimeSeries)
            if self.step==0:
                self.StartValues(StaticData, TimeSeries)
            else:
                self.setFirstControlElement(StaticData,TimeSeries)
                self.calcImplicitUpwind(StaticData, TimeSeries,ProgramModes)
            self.calcT(StaticData,TimeSeries)
            self.setTValues(StaticData,TimeSeries,Tsupply)
            T=self.finishstep()
            return T
        
        elif ProgramModes.TemperatureDynamicMode in ['ImplicitUpwind2']:
            self.spanControlVolumes(StaticData)
            self.buildSE(StaticData,TimeSeries)
            if self.step==0:
                self.StartValues(StaticData, TimeSeries)
            else:
                self.setFirstControlElement(StaticData,TimeSeries)
                self.calcImplicitUpwind2(StaticData, TimeSeries,ProgramModes)
            self.calcT(StaticData,TimeSeries)
            self.setTValues(StaticData,TimeSeries,Tsupply)
            T=self.finishstep()
            return T
        
            
            
    def spanControlVolumes(self,StaticData):
        """spans Controls Volume Grid over all Pipes at least 2 per Pipe"""
        self.nbControlVolumes=np.ceil(StaticData.Pipes_l/self.delta_x_max).astype(int)+1
        self.ControlLength=StaticData.Pipes_l/(self.nbControlVolumes-1)
        
        
    def StartValues(self,StaticData,TimeSeries):
        """ Define the start values in the ControlVolumes"""
        TE=TimeSeries.T[StaticData.startnode]
        TA=TimeSeries.T[StaticData.targetnode]
        for p in range(len(TE)):
            self.TPipe.append([])
            start_temp = TE[p]
            target_temp = TA[p]
            for i in range(self.nbControlVolumes[p]):
                self.TPipe[p].append([])
                # Interpolation factor (i / (num_volumes - 1)) ensures values are spaced evenly between start and target
                interpolation_factor = i / (self.nbControlVolumes[p] - 1) if self.nbControlVolumes[p] > 1 else 0
            
                interpolated_temp = start_temp + interpolation_factor * (target_temp - start_temp)                
                self.TPipe[p][i]=interpolated_temp
        self.TPipe_old=self.TPipe
        self.TPipe_older=self.TPipe
                
    def calcT(self,StaticData,TimeSeries):
        """ Calculate the tempreture at the Nodes from the last ControlVolumes in the pippes"""
        TA=np.zeros(StaticData.nbPipes)
        mask = TimeSeries.v > 0
        TA[mask]=np.array(([last for *_, last in [item for item, keep in zip(self.TPipe, mask) if keep]])).flatten()
        TA[~mask]=np.array(([first for first,*_ in [item for item, keep in zip(self.TPipe, ~mask) if keep]])).flatten()
        
        TA=np.concatenate((TA,TimeSeries.T[self.start[StaticData.nbPipes:]]))
        # Compute the numerator and denominator
        numerator = self.Se @ (TimeSeries.W[:StaticData.nbPipes + StaticData.nbPumps] * TA).astype(float)
        denominator = self.Se @ TimeSeries.W[:StaticData.nbPipes + StaticData.nbPumps].astype(float)
        
        # Safeguard the denominator to avoid division by zero
        denominator = np.where(denominator == 0, np.inf, denominator).astype(float)
        
        # Perform the division
        self.Tk = numerator / denominator
        
    def setFirstControlElement(self,StaticData,TimeSeries):
        """sets First Element in Pipe from Node before"""
        for p in range(StaticData.nbPipes):
            if TimeSeries.v[p]>=0:
                self.TPipe[p][0]=TimeSeries.T[self.start[p]]
                                            
            else:
                self.TPipe[p][-1]=TimeSeries.T[self.start[p]]
        
            

        
    def calcImplicitUpwind(self,StaticData,TimeSeries,ProgramModes):
        """Calculates one ImplicitUpwind step"""
        if ProgramModes.HeatLossMode  in [ 'default','simple']:
            self.T_pipe_mean=np.zeros(StaticData.nbPipes)
            for n in range(StaticData.nbPipes):
                self.T_pipe_mean[n]=1/2*(TimeSeries.T[StaticData.startnode[n]]+TimeSeries.T[StaticData.targetnode[n]])
            TimeSeries.Q_loss_Pipe=StaticData.Pipes_l/TimeSeries.Ur*(self.T_pipe_mean-TimeSeries.T_out)
            k=1/(StaticData.Pipes_A*TimeSeries.rho*TimeSeries.cp)/TimeSeries.Ur
            for p in range(StaticData.nbPipes):
                if TimeSeries.v[p]>=0:
                    for i in range(1,self.nbControlVolumes[p]):
        
                        self.TPipe[p][i]=\
                        (self.TPipe[p][i] +\
                        TimeSeries.v[p] * self.ti / self.ControlLength[p] * self.TPipe[p][i-1] +\
                        self.ti*k[p]*TimeSeries.T_out)/\
                        (1+
                         TimeSeries.v[p]* self.ti / self.ControlLength[p]+
                         k[p]*self.ti)
                else:
                    for i in reversed(range(0,self.nbControlVolumes[p]-1)):
                        self.TPipe[p][i]=\
                        (self.TPipe[p][i] -\
                        TimeSeries.v[p] * self.ti / self.ControlLength[p] * self.TPipe[p][i+1] +\
                        self.ti*k[p]*TimeSeries.T_out)/\
                        (1-
                         TimeSeries.v[p]* self.ti / self.ControlLength[p]+
                         k[p]*self.ti)
                            
        
        if ProgramModes.HeatLossMode not in [ 'default','simple']:
            
            self.T_pipe_mean=np.zeros(StaticData.nbPipes)
            for n in range(StaticData.nbPipes):
                self.T_pipe_mean[n]=1/2*(self.TPipe[n][0]+self.TPipe[n][-1])
            
            if ProgramModes.HeatLossNextPipe in ['No']:
                T_supply_avg=np.average(self.Tk[StaticData.supplyNodes])
                T_return_avg=np.average(self.Tk[StaticData.returnNodes])
            if ProgramModes.HeatLossNextPipe in ['Yes']:
                T_next=self.T_pipe_mean[StaticData.Pipes_next]
            
            #Q_loss=np.zeros(StaticData.nbPipes)
            
            for p in range(StaticData.nbPipes):
                
                if ProgramModes.HeatLossMode  in ['Norm']:
                    if ProgramModes.HeatLossNextPipe in ['No']:
                        if p in StaticData.supply_pipes:
                            Q_loss=TimeSeries.heatLossNorm(np.average(self.TPipe[p][:]),T_return_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_st_out"]) 
                        if p in StaticData.return_pipes:
                            Q_loss=TimeSeries.heatLossNorm(np.average(self.TPipe[p][:]),T_supply_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_st_out"]) 
                     
                    if ProgramModes.HeatLossNextPipe in ['Yes']:
                        Q_loss=TimeSeries.heatLossNorm(self.T_pipe_mean[p],T_next[p],TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_st_out"])
                        
                    
                    
                
                if ProgramModes.HeatLossMode  in ['EMCM']:
                    if ProgramModes.HeatLossNextPipe in ['No']:
                        if p in StaticData.supply_pipes:
                            Q_loss=TimeSeries.heatLossEMCM(np.average(self.TPipe[p][:]),T_return_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]) 
                        if p in StaticData.return_pipes:
                            Q_loss=TimeSeries.heatLossEMCM(np.average(self.TPipe[p][:]),T_supply_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]) 
                        
                        
                    if ProgramModes.HeatLossNextPipe in ['Yes']:
                        Q_loss=TimeSeries.heatLossEMCM(self.T_pipe_mean[p],T_next[p],TimeSeries.T_out,StaticData.pipe_dimensions[StaticData.Pipes_DN[p]])
                            
                S=-Q_loss/(TimeSeries.cp[p]*TimeSeries.rho[p]*StaticData.Pipes_A[p])
                TimeSeries.Q_loss_Pipe[p]=Q_loss*StaticData.Pipes_l[p]
                    
                    
                    
                if TimeSeries.v[p]>=0:
                    for i in range(1,self.nbControlVolumes[p]):
                        self.TPipe[p][i]=\
                        (self.TPipe[p][i] *self.ControlLength[p] +\
                        S * self.ti *self.ControlLength[p] +\
                        self.TPipe[p][i-1]*TimeSeries.v[p]*self.ti)/\
                        (self.ti*TimeSeries.v[p]+self.ControlLength[p])
                        
                else:
                    for i in reversed(range(0,self.nbControlVolumes[p]-1)):
                        self.TPipe[p][i]=\
                        (self.TPipe[p][i]*self.ControlLength[p] +\
                        S * self.ti *self.ControlLength[p] -\
                        self.TPipe[p][i+1]*TimeSeries.v[p]*self.ti)/\
                        (-self.ti*TimeSeries.v[p]+self.ControlLength[p])
                        
    def calcImplicitUpwind2(self,StaticData,TimeSeries,ProgramModes):
        """Calculates one ImplicitUpwind step"""
        if ProgramModes.HeatLossMode  in [ 'default','simple']:
            self.T_pipe_mean=np.zeros(StaticData.nbPipes)
            for n in range(StaticData.nbPipes):
                self.T_pipe_mean[n]=1/2*(TimeSeries.T[StaticData.startnode[n]]+TimeSeries.T[StaticData.targetnode[n]])
                TimeSeries.Q_loss_Pipe=StaticData.Pipes_l/TimeSeries.Ur*(self.T_pipe_mean-TimeSeries.T_out)
            k=1/(StaticData.Pipes_A*TimeSeries.rho*TimeSeries.cp)/TimeSeries.Ur
            for p in range(StaticData.nbPipes):
                if TimeSeries.v[p] >= 0:  # Positive velocity, use backward upwind in space
                    # First-order upwind scheme for the second control volume (i = 1)
                    self.TPipe[p][1] = \
                        (
                            4 * self.TPipe_old[p][1] - self.TPipe_older[p][1] + \
                            2 * TimeSeries.v[p] * self.ti / self.ControlLength[p] * \
                            self.TPipe[p][0]  + \
                            2 * self.ti * k[p] * TimeSeries.T_out
                        ) / \
                        (3 + 2 * k[p] * self.ti + 2* TimeSeries.v[p] * self.ti / self.ControlLength[p])
                    
                    # Second-order upwind scheme for control volumes (i >= 2)
                    for i in range(2, self.nbControlVolumes[p]):
                        self.TPipe[p][i] = \
                            (
                                4 * self.TPipe_old[p][i] - self.TPipe_older[p][i] + \
                                TimeSeries.v[p] * self.ti / self.ControlLength[p] * \
                                (4 * self.TPipe[p][i-1] - self.TPipe[p][i-2]) + \
                                2 * self.ti * k[p] * TimeSeries.T_out
                            ) / \
                            (3 + 2 * k[p] * self.ti + 3 * TimeSeries.v[p] * self.ti / self.ControlLength[p])
            
                else:  # Negative velocity, use forward upwind in space
                    # First-order upwind scheme for the second control volume (i = 1)
                    self.TPipe[p][self.nbControlVolumes[p] - 2] = \
                        (
                            4 * self.TPipe_old[p][self.nbControlVolumes[p] - 2] - \
                            self.TPipe_older[p][self.nbControlVolumes[p] - 2] + \
                            2 * abs(TimeSeries.v[p]) * self.ti / self.ControlLength[p] * \
                            self.TPipe[p][self.nbControlVolumes[p] - 1] + \
                            2 * self.ti * k[p] * TimeSeries.T_out
                        ) / \
                        (3 + 2 * k[p] * self.ti + 2* abs(TimeSeries.v[p]) * self.ti / self.ControlLength[p])
            
                    # Second-order upwind scheme for control volumes (i <= self.nbControlVolumes[p] - 3)
                    for i in reversed(range(0, self.nbControlVolumes[p] - 2)):
                        self.TPipe[p][i] = \
                            (
                                4 * self.TPipe_old[p][i] - self.TPipe_older[p][i] + \
                                abs(TimeSeries.v[p]) * self.ti / self.ControlLength[p] * \
                                (4 * self.TPipe[p][i+1] - self.TPipe[p][i+2]) + \
                                2 * self.ti * k[p] * TimeSeries.T_out
                            ) / \
                            (3 + 2 * k[p] * self.ti + 3 * abs(TimeSeries.v[p]) * self.ti / self.ControlLength[p])



                

                            
        
        if ProgramModes.HeatLossMode not in [ 'default','simple']:
            
            self.T_pipe_mean=np.zeros(StaticData.nbPipes)
            for n in range(StaticData.nbPipes):
                self.T_pipe_mean[n]=1/2*(self.TPipe[n][0]+self.TPipe[n][-1])
            
            if ProgramModes.HeatLossNextPipe in ['No']:
                T_supply_avg=np.average(self.Tk[StaticData.supplyNodes])
                T_return_avg=np.average(self.Tk[StaticData.returnNodes])
            if ProgramModes.HeatLossNextPipe in ['Yes']:
                T_next=self.T_pipe_mean[StaticData.Pipes_next]
            
            #Q_loss=np.zeros(StaticData.nbPipes)
            
            for p in range(StaticData.nbPipes):
                
                if ProgramModes.HeatLossMode  in ['Norm']:
                    if ProgramModes.HeatLossNextPipe in ['No']:
                        if p in StaticData.supply_pipes:
                            Q_loss=TimeSeries.heatLossNorm(np.average(self.TPipe[p][:]),T_return_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_st_out"]) 
                        if p in StaticData.return_pipes:
                            Q_loss=TimeSeries.heatLossNorm(np.average(self.TPipe[p][:]),T_supply_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_st_out"]) 
                     
                    if ProgramModes.HeatLossNextPipe in ['Yes']:
                        Q_loss=TimeSeries.heatLossNorm(self.T_pipe_mean[p],T_next[p],TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_ins_out"],StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]["d_st_out"])
                        
                    
                    
                
                if ProgramModes.HeatLossMode  in ['EMCM']:
                    if ProgramModes.HeatLossNextPipe in ['No']:
                        if p in StaticData.supply_pipes:
                            Q_loss=TimeSeries.heatLossEMCM(np.average(self.TPipe[p][:]),T_return_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]) 
                        if p in StaticData.return_pipes:
                            Q_loss=TimeSeries.heatLossEMCM(np.average(self.TPipe[p][:]),T_supply_avg,TimeSeries.T_out, StaticData.pipe_dimensions[StaticData.Pipes_DN[p]]) 
                        
                        
                    if ProgramModes.HeatLossNextPipe in ['Yes']:
                        Q_loss=TimeSeries.heatLossEMCM(self.T_pipe_mean[p],T_next[p],TimeSeries.T_out,StaticData.pipe_dimensions[StaticData.Pipes_DN[p]])
                            
                S=-Q_loss/(TimeSeries.cp[p]*TimeSeries.rho[p]*StaticData.Pipes_A[p])
                TimeSeries.Q_loss_Pipe[p]=Q_loss*StaticData.Pipes_l[p]  
                    
                    
                    
                if TimeSeries.v[p] >= 0:  # Positive velocity (backward upwind)
                
                    # First-order upwind for i = 1 (second control volume)
                    self.TPipe[p][1] = (
                        (4 * self.TPipe_old[p][1] - self.TPipe_older[p][1]) * self.ControlLength[p] +  # BDF2 time step
                        2 * S * self.ti * self.ControlLength[p] +  # Source term
                        2 * self.TPipe[p][0] * TimeSeries.v[p] * self.ti  # 1st order upwind in space (i = 1)
                    ) / (
                        (2 * self.ti * TimeSeries.v[p] + 3 * self.ControlLength[p])  # Denominator
                    )
            
                    # Second-order upwind for control volumes i >= 2
                    for i in range(2, self.nbControlVolumes[p]):
                        self.TPipe[p][i] = (
                            (4 * self.TPipe_old[p][i] - self.TPipe_older[p][i]) * self.ControlLength[p] +  # BDF2 time step
                            2 * S * self.ti * self.ControlLength[p] +  # Source term
                            (4 * self.TPipe[p][i-1] - self.TPipe[p][i-2]) * TimeSeries.v[p] * self.ti  # 2nd order upwind in space
                        ) / (
                            (3 * TimeSeries.v[p] * self.ti + 3 * self.ControlLength[p])  # Denominator
                        )
                
                else:  # Negative velocity (forward upwind)
                
                    # First-order upwind for i = nbControlVolumes[p] - 2 (second-to-last control volume)
                    self.TPipe[p][self.nbControlVolumes[p] - 2] = (
                        (4 * self.TPipe_old[p][self.nbControlVolumes[p] - 2] - self.TPipe_older[p][self.nbControlVolumes[p] - 2]) * self.ControlLength[p] +  # BDF2 time step
                        2 * S * self.ti * self.ControlLength[p] +  # Source term
                        2* self.TPipe[p][self.nbControlVolumes[p] - 1] * abs(TimeSeries.v[p]) * self.ti  # 1st order upwind in space
                    ) / (
                        (2 * abs(TimeSeries.v[p]) * self.ti + 3 * self.ControlLength[p])  # Denominator
                    )
            
                    # Second-order upwind for control volumes i <= nbControlVolumes[p] - 3
                    for i in reversed(range(0, self.nbControlVolumes[p] - 2)):
                        self.TPipe[p][i] = (
                            (4 * self.TPipe_old[p][i] - self.TPipe_older[p][i]) * self.ControlLength[p] +  # BDF2 time step
                            2 * S * self.ti * self.ControlLength[p] +  # Source term
                            (4 * self.TPipe[p][i+1] - self.TPipe[p][i+2]) * abs(TimeSeries.v[p]) * self.ti  # 2nd order upwind in space
                        ) / (
                            (3 * abs(TimeSeries.v[p]) * self.ti + 3 * self.ControlLength[p])  # Denominator
                        )
    
        self.TPipe_older = [pipe.copy() for pipe in self.TPipe_old]
        self.TPipe_old = [pipe.copy() for pipe in self.TPipe]
                            
    def calculatetime(self,StaticData,TimeSeries):
        """ Calculates the time where the Medium at every was last at a node"""
        # if self.step>0:
        #     idx=np.where(self.flowtime<=self.ti)[0]
        #     not_idx=np.where(self.flowtime>self.ti)[0]
        #     self.flowtime[idx]=StaticData.Pipes_l[idx]/abs(TimeSeries.v[idx])
        #     self.flowtime[not_idx]=(1-self.ti/self.flowtime[not_idx])*self.flowtime[not_idx]+self.ti*StaticData.Pipes_l[not_idx]/self.flowtime[not_idx]/abs(TimeSeries.v[not_idx])
            
        # else:
        #     self.flowtime=StaticData.Pipes_l/abs(TimeSeries.v)
            
        # self.te=self.step*self.ti-self.flowtime
        
        # self.s= self.ti*np.floor(self.te/self.ti)
        
        self.s=np.ones(StaticData.nbPipes)*self.step*self.ti
        self.te=np.zeros(StaticData.nbPipes)
        self.change_direction=[]
        for i in range(StaticData.nbPipes):
            if self.v[-1][i]==0:
                self.te[i]=self.step*self.ti
                continue
            Z=-StaticData.Pipes_l[i]
            while Z<0:
                if int(self.s[i]/self.ti)>=0:
                    if self.v[int(self.s[i]/self.ti)][i]==0:
                        self.te[i]=self.step*self.ti
                        break
                    Z=Z+self.ti*self.int_faster*self.v[int(self.s[i]/self.ti)][i]*np.sign(self.v[-1][i])
                else:
                    if self.v[0][i]==0:
                        self.te[i]=self.step*self.ti
                        break
                    Z=Z+self.ti*self.int_faster*abs(self.v[0][i])
                

                    
                self.s[i]=self.s[i]-self.ti*self.int_faster
                if Z<-StaticData.Pipes_l[i]:
                    print('change direction of pipe', i, 'at timestep', self.step)
                    self.change_direction.append(i)
                    break


            if Z>0:
                self.te[i]=self.s[i]+Z/abs(self.v[-1][i])
            else:
                self.te[i]=self.s[i]+abs(StaticData.Pipes_l[i]+Z)/abs(self.v[-1][i])
            
    def interpolate_temperature(self,te, t_values, T_values):
        """Interpolates the tempreture for time te from TimeValues where Tempretures exists"""
        poly_coeff=np.polyfit(t_values, T_values,1)
        return np.poly1d(poly_coeff)(te)
        
        #return lagrange(t_values, T_values)(te)
            
    def calculateTE(self,StaticData,TimeSeries):
        """Calculate the input tempreture of every Pipe"""
        self.TE=np.zeros(StaticData.nbPipes+StaticData.nbPumps)
        if self.step==0:
            self.T_dynamic.append(TimeSeries.T[:StaticData.nbNodes])
        
        for i in range(StaticData.nbPipes): 
            t_m1=self.s[i]
            t_m2=self.s[i]+self.ti*self.int_faster
            t_strich=(t_m1+t_m2)/2
            
            
            
            if round((t_m1 - 1 * self.ti*self.int_faster) / self.ti) > 0 and round((t_m1 + self.ti*self.int_faster) / self.ti) <= self.step - 1:
                if self.te[i] > t_strich:
                    
                    # t_values = [round((t_m2 - 4 * self.ti) / self.ti), round((t_m2 - 3 * self.ti) / self.ti), round((t_m2 - 2 * self.ti) / self.ti), round((t_m2 - self.ti) / self.ti), round(t_m2 / self.ti), round((t_m2 + self.ti) / self.ti)]
                    # T_values = [self.T_dynamic[round((t_m2 - 4 * self.ti) / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 - 3 * self.ti) / self.ti)][self.start[i]],
                    #             self.T_dynamic[round((t_m2 - 2 * self.ti) / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 - self.ti) / self.ti)][self.start[i]],
                    #             self.T_dynamic[round(t_m2 / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 + self.ti) / self.ti)][self.start[i]]
                    #             ]
                    # t_values = [round((t_m2 - 2 * self.ti) / self.ti), round((t_m2 - self.ti) / self.ti), round(t_m2 / self.ti), round((t_m2 + self.ti) / self.ti)]
                    # T_values = [self.T_dynamic[round((t_m2 - 2 * self.ti) / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 - self.ti) / self.ti)][self.start[i]],
                    #             self.T_dynamic[round(t_m2 / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 + self.ti) / self.ti)][self.start[i]]
                    #             ]

                    t_values = [round((t_m2 - 1 * self.ti*self.int_faster) / self.ti), round((t_m2 - 0 * self.ti*self.int_faster) / self.ti), round((t_m2 + 1 * self.ti*self.int_faster) / self.ti)]
                    T_values = [self.T_dynamic[round((t_m2 - 1 * self.ti*self.int_faster) / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 - 0 * self.ti*self.int_faster) / self.ti)][self.start[i]],
                                self.T_dynamic[round((t_m2 + 1 * self.ti*self.int_faster) / self.ti)][self.start[i]]
                                ]
                    
                    # t_values = [round((t_m2 - 1 * self.ti*self.int_faster) / self.ti), round((t_m2 - 0 * self.ti*self.int_faster) / self.ti)]
                    # T_values = [self.T_dynamic[round((t_m2 - 1 * self.ti*self.int_faster) / self.ti)][self.start[i]], self.T_dynamic[round((t_m2 - 0 * self.ti*self.int_faster) / self.ti)][self.start[i]]   
                    #             ]
            
                if self.te[i] <= t_strich:
                    # t_values = [round((t_m1 - 4 * self.ti) / self.ti), round((t_m1 - 3 * self.ti) / self.ti), round((t_m1 - 2 * self.ti) / self.ti), round((t_m1 - self.ti) / self.ti), round(t_m1 / self.ti), round((t_m1 + self.ti) / self.ti)]
                    # T_values = [self.T_dynamic[round((t_m1 - 4 * self.ti) / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 - 3 * self.ti) / self.ti)][self.start[i]],
                    #             self.T_dynamic[round((t_m1 - 2 * self.ti) / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 - self.ti) / self.ti)][self.start[i]],
                    #             self.T_dynamic[round(t_m1 / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 + self.ti) / self.ti)][self.start[i]]
                    #           ]
                    # t_values = [round((t_m1 - 2 * self.ti) / self.ti), round((t_m1 - self.ti) / self.ti), round(t_m1 / self.ti), round((t_m1 + self.ti) / self.ti)]
                    # T_values = [self.T_dynamic[round((t_m1 - 2 * self.ti) / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 - self.ti) / self.ti)][self.start[i]],
                    #             self.T_dynamic[round(t_m1 / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 + self.ti) / self.ti)][self.start[i]]
                    #             ]

                    t_values = [round((t_m1 - 1 * self.ti*self.int_faster) / self.ti), round((t_m1 - 0 * self.ti*self.int_faster) / self.ti), round((t_m1 + 1 * self.ti*self.int_faster) / self.ti)]
                    T_values = [self.T_dynamic[round((t_m1 - 1 * self.ti*self.int_faster) / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 - 0 * self.ti*self.int_faster) / self.ti)][self.start[i]],
                                self.T_dynamic[round((t_m1 + 1 * self.ti*self.int_faster) / self.ti)][self.start[i]]
                                ]
                    # t_values = [round((t_m1 - 1 * self.ti*self.int_faster) / self.ti), round((t_m1 - 0 * self.ti*self.int_faster) / self.ti)]
                    # T_values = [self.T_dynamic[round((t_m1 - 1 * self.ti*self.int_faster) / self.ti)][self.start[i]], self.T_dynamic[round((t_m1 - 0 * self.ti*self.int_faster) / self.ti)][self.start[i]]   
                    #             ]
                
                try:
                    self.TE[i] = self.interpolate_temperature(self.te[i]/self.ti, t_values, T_values)
                except:
                    self.TE[i] = self.T_dynamic[0][self.target[i]]
                    
            elif round((t_m1 + self.ti*self.int_faster) / self.ti) > self.step - 1:
                self.TE[i]=self.T_dynamic[-1][self.start[i]]
                
                
            else:
                self.TE[i]=self.T_dynamic[0][self.target[i]]
            
        self.TE[StaticData.nbPipes:]=self.T_dynamic[-1][self.start[StaticData.nbPipes:]]
        
        self.TE=np.clip(self.TE, 0, 150)
            
    def calculateTcoeff(self,StaticData,TimeSeries):
        """ Calculate the T at the end of each * heat flow """ 
        self.ta=self.step*self.ti
        exp=np.exp(-1*(self.ta-self.te)/\
        (StaticData.Pipes_A*TimeSeries.cp*TimeSeries.rho*TimeSeries.Ur)).astype(float) # e^(-Ur*l/d/cp/rho/v)
        exp=np.concatenate((exp, np.ones(StaticData.nbPumps)),axis=0)
        self.Tcoeff=(TimeSeries.T_out+(self.TE-TimeSeries.T_out)*exp)*abs(TimeSeries.W[:StaticData.nbPipes+StaticData.nbPumps])
        
    def buildSE(self,StaticData,TimeSeries):
        """ get the diretion of flow in each Pipe and get the Node entry matrix """
        self.start=np.concatenate((StaticData.startnode,StaticData.Pumps_startnode)) #startnodes Pipes and Pumps
        self.target=np.concatenate((StaticData.targetnode,StaticData.Pumps_targetnode)) #self.targetnodes Pipes and Pumps
        idx=np.where(TimeSeries.m_i<0) #idx[0] with negative massflow --> temperature change in the other direction
        if len(idx[0])>0:
            self.start[idx[0]],self.target[idx[0]]=self.target[idx[0]],self.start[idx[0]]
        try:
            self.start[self.change_direction]=self.target[self.change_direction]
        except:
            pass

        
        self.Se=np.zeros((StaticData.nbNodes,self.target.shape[0]))
        
        for n in range(self.target.shape[0]): # iterarte through start/self.target maybe check if same length else error
            self.Se[self.target[n],n]=1
            
    def calcTk(self,StaticData,TimeSeries):
        """ calulate the Tempreture at each Node"""
        # Compute numerator and denominator
        numerator = self.Se @ self.Tcoeff
        denominator = self.Se @ TimeSeries.W[:StaticData.nbPipes + StaticData.nbPumps]
        
        # Safeguard the denominator to avoid division by zero
        denominator = np.where(denominator == 0, np.inf, denominator)
        
        # Perform the division
        self.Tk = (numerator / denominator).astype(float)
        
    def setTValues(self,StaticData,TimeSeries,Tsupply):
        """ Set The Tempreture Values at Network flow in Nodes  Suplliers Supply Network and consumers in the rerutnm network"""
        
        Tset=np.concatenate((StaticData.vorlauf_SupplieNodes,StaticData.ruecklauf_consumptionnode_unique))
        # Tsupply=110*np.ones(StaticData.nbSupplies)
        # if self.step>1000:
            
        #     Tsupply[0]=Tsupply[0]#+20#*np.sin((self.step)/80)
        Tconsumer=StaticData.Outlets_Tref
        W=TimeSeries.cp_consumer*TimeSeries.m_consumer
        # Create a mask for non-zero elements in W
        non_zero_mask = W != 0
        
        # Initialize safe_division with zeros or a fallback value
        safe_division = np.zeros_like(TimeSeries.Q)
        
        # Perform division only for non-zero W
        safe_division[non_zero_mask] = TimeSeries.Q[non_zero_mask] / W[non_zero_mask]
        
        # Compute the mask
        mask = self.Tk[StaticData.vorlauf_consumptionnode] - Tconsumer < safe_division -5
        
        # Explicitly set mask to True where W is zero
        mask[~non_zero_mask] = False
        

        Tconsumer[mask]=self.Tk[StaticData.vorlauf_consumptionnode][mask]-TimeSeries.Q[mask]/W[mask]
        
        # Check if self.Se_consumer already exists
        if not hasattr(self, 'Se_consumer'):
            # Initialize Se_consumer as a zero matrix
            self.Se_consumer = np.zeros((StaticData.nbNodes, StaticData.nbOutlets))
        
            # Iterate through StaticData.nbOutlets
            for n in range(StaticData.nbOutlets):
                # Assign 1 to the appropriate element if indices match
                self.Se_consumer[StaticData.ruecklauf_consumptionnode[n], n] = 1
        setTValues_idx=np.concatenate((StaticData.vorlauf_SupplieNodes,StaticData.ruecklauf_consumptionnode_unique))
        # Compute numerator and denominator
        numerator = self.Se_consumer @ (W*Tconsumer)
        denominator = self.Se_consumer @ W
        
        # Safeguard the denominator to avoid division by zero
        denominator = np.where(denominator == 0, np.inf, denominator)
        
        # Perform the division
        setTValues = numerator / denominator
        setTValues[StaticData.vorlauf_SupplieNodes]=Tsupply
        self.Tk[setTValues_idx]=setTValues[setTValues_idx]
        self.Tk[np.isnan(self.Tk)]=30
        

            
        
        
    def finishstep(self):
        """ end the Timestep reurn Node tempreture
        inkrement step """
        self.step+=1
        self.T_dynamic.append(self.Tk)
        return self.Tk
        
            