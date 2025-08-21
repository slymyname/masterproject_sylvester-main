# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:32:40 2024

@author: sander
"""

from scipy.optimize import minimize

class OptimizeStationary:
    def __init__(self):
        self.weight_Heat=1
        self.weight_Pump=1/5
    
    def f(self,PumpPower,HeatLoss):
        
        f=self.weight_Pump * PumpPower + self.weight_Heat * HeatLoss
        return f/100000
    

    
    def f_eval(self,TSupply, *args):
        # Unpack the additional arguments
        StaticData, TimeSeries, hydraulicMesh, ProgramModes, LoadMeasValues = args
        
        # Call updateResults
        self.updateResults(StaticData, TimeSeries, TSupply, hydraulicMesh, ProgramModes, LoadMeasValues)
        
        # Return the result of the objective function
        return self.f(sum(TimeSeries.hydraulicPower), TimeSeries.HeatLoses)

    
    def updateResults(self,StaticData,TimeSeries,TSupply,hydraulicMesh,ProgramModes,LoadMeasValues):
        if abs(TimeSeries.T[StaticData.vorlauf_SupplieNodes]- TSupply)>0.01:    
            TimeSeries.timestepStationary(StaticData,hydraulicMesh,ProgramModes,LoadMeasValues,TSupply,False)

        
    def optimize(self,TSupply,StaticData, TimeSeries, hydraulicMesh, ProgramModes,LoadMeasValues):
        self.Result = minimize(self.f_eval, TSupply, args=(StaticData, TimeSeries, hydraulicMesh, ProgramModes, LoadMeasValues),
        method='SLSQP', bounds=[(90, 140)],
        jac='2-point',
        #constraints = nonlinearconstraints,
        options = {'ftol': 1e-16, 'finite_diff_rel_step':1e-3, 'disp': True, 'maxiter': 120})
        return self.Result.x
