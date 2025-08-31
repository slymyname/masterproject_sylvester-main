# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 09:49:29 2024

@author: sander
"""

import matplotlib.pyplot as plt
import numpy as np

class PlotResult:
    
    def __init__(self,ProgramModes):
        self.T=[]
        self.m_i=[]
        self.m_consumer=[]
        self.p=[]
        self.T_s=[]
        self.HL1=[]
        self.HL2=[]
        self.PumpPower=[]
    def plotT(self, ax, TimeSeries, ProgramModes):
        arrayT = np.vstack(self.T)
        if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
            ax.plot(TimeSeries.datetimeStationary[:], arrayT[:, :],drawstyle='steps-post')
        elif ProgramModes.TempratureCalculationMode=='dynamic':
            ax.plot(TimeSeries.datetimeDynamic[:], arrayT[:, :],drawstyle='steps-post')
        ax.set_title('Temperature Plot')
        ax.set_ylabel(r'Temperature in $°C$')
        ax.grid(True)
        ax.set_ylim([np.min(arrayT)*0.9,np.max(arrayT)*1.1])
        ax.set_xlim(TimeSeries.datetimeStationary[0],TimeSeries.datetimeStationary[-1])
        ax.axhline(y=80, color='g', linestyle='--', label='Min Supply Temperature (80 °C)')
        ax.axhline(y=150, color='m', linestyle='--', label='Max Supply Temperature (150 °C)')
        ax.legend()

    def plotp(self, ax, TimeSeries, ProgramModes):
        arrayp = np.vstack(self.p)
        ax.plot(TimeSeries.datetimeStationary[:], arrayp[:, :],drawstyle='steps-post')
        ax.set_title('Pressure Plot')
        ax.set_ylabel(r'Pressure in $bar$')
        ax.grid(True)
        ax.set_ylim([np.min(arrayp)*0.9,np.max(arrayp)*1.1])
        ax.axhline(y=17.5, color='r', linestyle='--', label='Max Outlet Pressure (17.5 bar)')
        ax.axhline(y=4, color='b', linestyle='--', label='Min Inlet Pressure (4 bar)')
        ax.legend()
        
    def plotm_i(self, ax, TimeSeries, ProgramModes):
        arraym_i = np.vstack(self.m_i)
        ax.plot(TimeSeries.datetimeStationary[:], arraym_i[:, :],drawstyle='steps-post')
        ax.set_title('Massflow Plot')
        ax.set_ylabel(r'Massflow in $\frac{kg}{s}$')
        ax.grid(True)
        ax.set_ylim([np.min(arraym_i)*0.9,np.max(arraym_i)*1.1])
        
    def plotHeatLosses(self, ax, TimeSeries, ProgramModes):
        arrayQloss = np.vstack(self.HL1)
        arrayPower = np.vstack(self.PumpPower)
        ax.plot(TimeSeries.datetimeStationary[:], arrayQloss[:, :],drawstyle='steps-post')
        ax.plot(TimeSeries.datetimeStationary[:], arrayPower[:, :],drawstyle='steps-post')
        ax.set_title('HeatLosses Plot')
        ax.set_ylabel(r'HeatLoss in $W$')
        ax.grid(True)
        ax.set_ylim([-1e5,1e6])
        
        
    def plot_combined(self, TimeSeries, ProgramModes):
        fig, axes = plt.subplots(2, 2, sharex=True, figsize=(8, 24))  # Create 2x2 subplots
        
        # Flatten the axes array to make it a 1D array
        axes = axes.flatten()
    
        # Call the separate plotting functions
        self.plotT(axes[0], TimeSeries, ProgramModes)
        self.plotp(axes[1], TimeSeries, ProgramModes)
        self.plotm_i(axes[2], TimeSeries, ProgramModes)
        self.plotHeatLosses(axes[3], TimeSeries, ProgramModes)
        
        
        # Customize shared x-ticks
        # ticks = range(1, len(TimeSeries.datetimeDynamic), int(60 * 60 / ProgramModes.dynamicTimestep))
        # labels = [TimeSeries.datetimeDynamic[i].strftime("%H") for i in ticks]
        # ax2.set_xticks(ticks)
        # ax2.set_xticklabels(labels)
        
        axes[0].set_xlabel('Time')  # Set x-label for the combined plot
    
        plt.tight_layout()  # Adjust subplots to fit into figure area.
        plt.show()
     
        
    def plot_single(self, plot_type, TimeSeries, ProgramModes):
        fig, ax = plt.subplots(figsize=(16, 9))  # Create a single subplot
        
        if plot_type == 'T':
            self.plotT(ax, TimeSeries, ProgramModes)
        elif plot_type == 'p':
            self.plotp(ax, TimeSeries, ProgramModes)
        elif plot_type == 'm':
            self.plotm_i(ax, TimeSeries, ProgramModes)
        elif plot_type == 'heatloss':
            self.plotHeatLosses(ax, TimeSeries, ProgramModes)
        else:
            raise ValueError("Invalid plot_type. Choose from 'T', 'p', 'm', or 'heatloss'.")

        ax.set_xlabel('Time')  # Set x-label for single plot
        plt.tight_layout()     # Adjust layout
        plt.show()
            
    def plot_single_choose(self, plot_type, TimeSeries, ProgramModes, line_indices=None):
        fig, ax = plt.subplots(figsize=(16, 9))  # Create a single subplot
        
        # Select and plot the desired data based on plot_type
        if plot_type == 'T':
            arrayT = np.vstack(self.T)
            selected_data = arrayT[:, line_indices] if line_indices else arrayT
            if ProgramModes.TempratureCalculationMode == 'default' or ProgramModes.TempratureCalculationMode == 'stationary':
                ax.plot(TimeSeries.datetimeStationary[:], selected_data, drawstyle='steps-post')
            elif ProgramModes.TempratureCalculationMode == 'dynamic':
                ax.plot(TimeSeries.datetimeDynamic[:], selected_data, drawstyle='steps-post')
            ax.set_title('Temperature Plot')
            ax.set_ylabel(r'Temperature in $°C$')
            
        elif plot_type == 'p':
            arrayp = np.vstack(self.p)
            selected_data = arrayp[:, line_indices] if line_indices else arrayp
            ax.plot(TimeSeries.datetimeStationary[:], selected_data, drawstyle='steps-post')
            ax.set_title('Pressure Plot')
            ax.set_ylabel(r'Pressure in $bar$')
            
        elif plot_type == 'm':
            arraym_i = np.vstack(self.m_i)
            selected_data = arraym_i[:, line_indices] if line_indices else arraym_i
            ax.plot(TimeSeries.datetimeStationary[:], selected_data, drawstyle='steps-post')
            ax.set_title('Massflow Plot')
            ax.set_ylabel(r'Massflow in $\frac{kg}{s}$')
            
        elif plot_type == 'heatloss':
            arrayQloss = np.vstack(self.HL1)
            arrayPower = np.vstack(self.PumpPower)
            
            # Separate handling for two arrays if line indices are provided
            selected_Qloss = arrayQloss[:, line_indices] if line_indices else arrayQloss
            selected_Power = arrayPower[:, line_indices] if line_indices else arrayPower
            
            ax.plot(TimeSeries.datetimeStationary[:], selected_Qloss, drawstyle='steps-post', label="Heat Loss")
            ax.plot(TimeSeries.datetimeStationary[:], selected_Power, drawstyle='steps-post', label="Pump Power")
            ax.set_title('HeatLosses Plot')
            ax.set_ylabel(r'HeatLoss in $W$')
            ax.legend()
            
        else:
            raise ValueError("Invalid plot_type. Choose from 'T', 'p', 'm', or 'heatloss'.")
        
        ax.set_xlabel('Time')
        ax.grid(True)
        plt.tight_layout()
        plt.show()
    
    def plot_meas_sim(self,TimeSeries,LoadMeasValues,StaticData,ProgramModes,measType,no):
        plt.figure()
        for Data in LoadMeasValues.Data:
            if measType=='m':
                AKZ=StaticData.mMeasurements[no,0]+'XQ20'
                for i, entry in enumerate(Data):
                    if entry.get("AKZ") == AKZ:
                        idx = i
                        break
                    
                ValuesCalc=np.array([m[StaticData.mMeasurements[no,1]] for m in self.m_i])
                TimeCalc=TimeSeries.datetimeStationary
                
            if measType=='T':
                AKZ=StaticData.TMeasurements[no,0]+'XQ01'
                for i, entry in enumerate(Data):
                    if entry.get("AKZ") == AKZ:
                        idx = i
                        break
                    
                ValuesCalc=np.array([T[StaticData.TMeasurements[no,1]] for T in self.T])
                if ProgramModes.TempratureCalculationMode=='default' or ProgramModes.TempratureCalculationMode=='stationary':
                    TimeCalc=TimeSeries.datetimeStationary
                elif ProgramModes.TempratureCalculationMode=='dynamic':
                    TimeCalc=TimeSeries.datetimeDynamic
                
                
            if measType=='p':
                AKZ=StaticData.pMeasurements[no,0]+'XQ01'
                for i, entry in enumerate(Data):
                    if entry.get("AKZ") == AKZ:
                        idx = i
                        break
                    
                ValuesCalc=np.array([p[StaticData.pMeasurements[no,1]] for p in self.p])
                TimeCalc=TimeSeries.datetimeStationary
                   
            ValuesMeas=Data[idx]['value']
            TimeMeas=Data[idx]['time']
            if measType=='m':
                ValuesMeas=ValuesMeas/3.6
            
    
            
            plt.plot(TimeMeas,ValuesMeas)
        plt.plot(TimeCalc,ValuesCalc)
        plt.xlabel('Datetime')
        plt.ylabel('Values')
        plt.grid(True)
        plt.xticks(rotation=45)
        plt.legend()
        
        # Show plot
        plt.show()
        
        
        
        
    
