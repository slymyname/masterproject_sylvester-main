# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:57:46 2024

@author: sander
"""
import time
import datetime

import numpy as np
import matplotlib.pyplot as plt
from StaticData import StaticData

from ProgramModes import ProgramModes
from hydraulicMesh import hydraulicMesh
from TimeSeries import TimeSeries
from LoadMeasValues import LoadMeasValues
from optimize import OptimizeStationary

from PlotResults import PlotResult

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
ProgramModes.TInput='None' # for Measurement Value at Supply 'Yes' or 'None' 
ProgramModes.Optimize='Yes' # 'Yes' or 'No'

ProgramModes.startTime=datetime.datetime(2022, 2, 1, 1, 0, 0, 0) #change startdate
ProgramModes.endTime=datetime.datetime(2022, 2, 1, 1 , 30 ,0, 0) #change enddate


StaticData=StaticData()
StaticData.readData(ProgramModes)
StaticData.convertAllData()

hydraulicMesh=hydraulicMesh()
hydraulicMesh.MeshPrep(StaticData)

TimeSeries=TimeSeries(StaticData,hydraulicMesh,ProgramModes)
TimeSeries.calcThermalResistanceVertical(StaticData)
TimeSeries.datestart(ProgramModes)

LoadMeasValues=LoadMeasValues()
if not(ProgramModes.PowerInput=='None' and ProgramModes.TInput=='None'):
    LoadMeasValues.loadmat(TimeSeries.datetimeStationary)
if ProgramModes.Optimize=='Yes':
    Optimize=OptimizeStationary()

#PlotResult.plot_meas_sim( TimeSeries, LoadMeasValues, StaticData, 'm', 0)

PlotResult=PlotResult(ProgramModes)

LoadMeasValues.TemperatureOutside(r'prm/produkt_tu_stunde_19970701_20231231_03379.txt')



Tsupply=np.array([120])  #here the Supply Tempreture can be set

if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
     
        
    # finalstep=False
    while TimeSeries.datetimeStationary[-1]<ProgramModes.endTime:
        
        LoadMeasValues.updateTout(TimeSeries, ProgramModes)
        if not(ProgramModes.PowerInput=='None' and ProgramModes.TInput=='None'):
            LoadMeasValues.loadmat(TimeSeries.datetimeStationary)
                
        if ProgramModes.TInput=='Yes':
            Tsupply=np.array([LoadMeasValues.SupplyTemperature(TimeSeries.datetimeStationary[-1])])
            
        if ProgramModes.Optimize=='Yes':
            Tsupply=Optimize.optimize(Tsupply, StaticData, TimeSeries, hydraulicMesh, ProgramModes,LoadMeasValues)
        
        TimeSeries.timestepStationary(StaticData,hydraulicMesh,ProgramModes,LoadMeasValues, Tsupply,True)
        PlotResult.T.append(TimeSeries.T)
        PlotResult.m_i.append(TimeSeries.m_i)
        PlotResult.p.append(TimeSeries.p)
        PlotResult.m_consumer.append(TimeSeries.m_consumer)
        PlotResult.HL1.append(TimeSeries.HeatLoses)
        PlotResult.PumpPower.append(TimeSeries.hydraulicPower)

        
        
        

        

elif ProgramModes.TempratureCalculationMode=='dynamic':
    LoadMeasValues.updateTout(TimeSeries, ProgramModes)
    if ProgramModes.TInput=='Yes':
        Tsupply=np.array([LoadMeasValues.SupplyTemperature(TimeSeries.datetimeStationary[-1])]) 
    TimeSeries.timestepStationary(StaticData,hydraulicMesh,ProgramModes,LoadMeasValues,Tsupply,False)
    PlotResult.m_i.append(TimeSeries.m_i.copy())
    PlotResult.p.append(TimeSeries.p.copy())
    PlotResult.T_s.append(TimeSeries.T.copy())
    PlotResult.HL1.append(TimeSeries.HeatLoses)
    PlotResult.PumpPower.append(TimeSeries.hydraulicPower)
    while TimeSeries.datetimeStationary[-1]<ProgramModes.endTime:


        while TimeSeries.datetimeDynamic[-1]<TimeSeries.datetimeStationary[-1]+TimeSeries.time_change_stationary:
            if ProgramModes.TInput=='Yes':
                Tsupply=np.array([LoadMeasValues.SupplyTemperature(TimeSeries.datetimeDynamic[-1])]) 
            TimeSeries.timestepDynamic(StaticData,ProgramModes,Tsupply)
            PlotResult.T.append(TimeSeries.T.copy())
            
        LoadMeasValues.updateTout(TimeSeries, ProgramModes)
        if not(ProgramModes.PowerInput=='None' and ProgramModes.TInput=='None'):
            LoadMeasValues.loadmat(TimeSeries.datetimeStationary)           
        if len(TimeSeries.datetimeStationary)>0:
            if ProgramModes.TInput=='Yes':
                Tsupply=np.array([LoadMeasValues.SupplyTemperature(TimeSeries.datetimeDynamic[-1])]) 
            TimeSeries.timestepDynamicHydraulic(StaticData, hydraulicMesh,ProgramModes,LoadMeasValues,Tsupply)
            PlotResult.m_i.append(TimeSeries.m_i.copy())
            PlotResult.p.append(TimeSeries.p.copy())
            PlotResult.m_consumer.append(TimeSeries.m_consumer)
            PlotResult.T_s.append(TimeSeries.T.copy())
            PlotResult.HL1.append(TimeSeries.HeatLoses)
            PlotResult.PumpPower.append(TimeSeries.hydraulicPower)
            
        
# record end time
end = time.time()
 
# print the difference between start 
# and end time in milli. secs
print("The time of execution of above program is :",
      (end-start), "s")


PlotResult.plot_combined(TimeSeries,ProgramModes)
PlotResult.plot_mpc_states(TimeSeries, StaticData)
print(sum(PlotResult.HL1)*ProgramModes.stationaryTimeStep/60*1e-6,'MWh')
