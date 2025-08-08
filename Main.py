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
ProgramModes.Optimize='No' # 'Yes' or 'No'

ProgramModes.startTime=datetime.datetime(2022, 2, 1, 0, 0, 0, 0) #change startdate
ProgramModes.endTime=datetime.datetime(2022, 2, 3, 0, 0, 0, 0)


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



Tsupply_values = [110.50538848, 110.30381051, 111.06077353, 110.59445023, 110.67591816,
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
Tsupply_gen=iter(Tsupply_values)

if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
     
        
    # finalstep=False
    while TimeSeries.datetimeStationary[-1]<ProgramModes.endTime:
        
        LoadMeasValues.updateTout(TimeSeries, ProgramModes)
        if not(ProgramModes.PowerInput=='None' and ProgramModes.TInput=='None'):
            LoadMeasValues.loadmat(TimeSeries.datetimeStationary)
                
        if ProgramModes.TInput=='Yes':
            Tsupply=np.array([LoadMeasValues.SupplyTemperature(TimeSeries.datetimeStationary[-1])])
        else:
            try:
                Tsupply=np.array([next(Tsupply_gen)])
            except:
                Tsupply=np.array([120])
            
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

print(sum(PlotResult.HL1)*ProgramModes.stationaryTimeStep/60*1e-6,'MWh')
PlotResult.plot_combined(TimeSeries,ProgramModes)
print(sum(PlotResult.HL1)*ProgramModes.stationaryTimeStep/60*1e-6,'MWh')