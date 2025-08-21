# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 09:57:27 2023

@author: sander
"""
import datetime
class ProgramModes:
    
    def __init__(self):
        self.SystemMode='default'
        self.MeasurementMode='default'
        self.NetworkMode='default'
        self.TempratureCalculationMode='default'
        self.TemperatureDynamicMode='default'
        self.HeatLossMode='default'
        self.HeatLossNextPipe='No'
        self.randomConsumerMode='default'
        self.stationaryTimeStep=int(15) #in minutes
        self.dynamicTimestep=int(20) # in seconds
        self.delta_x_max=30 #in m
        self.startTime=datetime.datetime(2024, 1, 1, 0, 0, 0, 0)
        self.endTime=datetime.datetime(2024, 1, 2, 0, 0, 0, 0)
        self.PipeFrictionMethod='ColebrookWhite'
        #self.Input='default'
        self.accuracy_factor_hydraulic=0.005
        self.accuracy_factor_thermic=0.005
        self.SimulationEngine='FOM'
        self.EnableSnapshotSaving = True
        
    