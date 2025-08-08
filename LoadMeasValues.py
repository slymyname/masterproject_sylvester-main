# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 09:06:14 2024

@author: sander
"""

from scipy.io import loadmat
from scipy.io.matlab import mat_struct
from datetime import datetime
import csv
import numpy as np

class LoadMeasValues:
    
    def __init__(self):
        self.idx_Nord_supply_time = 0
        self.last_idx_Nord_supply_time = 0
        self.idx_Nord_return_time = 0
        self.last_idx_Nord_return_time = 0
        
        self.idx_West_supply_time = 0
        self.last_idx_West_supply_time = 0
        self.idx_West_return_time = 0
        self.last_idx_West_return_time = 0
        
        self.idx_Ost_supply_time = 0
        self.last_idx_Ost_supply_time = 0
        self.idx_Ost_return_time = 0
        self.last_idx_Ost_return_time = 0
        
        self.idx_T_supply_time=1
        self.last_idx_T_supply_time=1
        
        self.Data=[]
    
    def loadmat(self,DateTimes):  
        if len(DateTimes)>1:
            if DateTimes[-1].day!=DateTimes[-2].day:
                current_day=DateTimes[-1]
                
                self.idx_Nord_supply_time = 0
                self.last_idx_Nord_supply_time = 0
                self.idx_Nord_return_time = 0
                self.last_idx_Nord_return_time = 0
                
                self.idx_West_supply_time = 0
                self.last_idx_West_supply_time = 0
                self.idx_West_return_time = 0
                self.last_idx_West_return_time = 0
                
                self.idx_Ost_supply_time = 0
                self.last_idx_Ost_supply_time = 0
                self.idx_Ost_return_time = 0
                self.last_idx_Ost_return_time = 0
                
                self.idx_T_supply_time=1
                self.last_idx_T_supply_time=1
                del self.SupplyVolumeflow_Nord_supply_index,self.SupplyVolumeflow_Nord_return_index 
                del self.SupplyVolumeflow_West_supply_index, self.SupplyVolumeflow_West_return_index
                del self.SupplyVolumeflow_Ost_supply_index,self.SupplyVolumeflow_Ost_return_index 
                del self.Temperature_Supply_index
                
                try:
                    datestring=current_day.strftime("%Y_%m_%d")
                    #Use the absolute path to your .mat file
                    mat_data = loadmat('D:\EnEffDaten\Data' + datestring + '.mat', struct_as_record=False, squeeze_me=True)
                    matobj=mat_data['Data']
                        # Convert the entire cell array of structs
                    self.Data.append( [self.matstruct_to_dict(s) for s in matobj] )
                except ValueError:
                    print("The current day is not found in the Data.")
        elif len(self.Data)==0:
            try:
                current_day=DateTimes[-1]
                datestring=current_day.strftime("%Y_%m_%d")
                #Use the absolute path to your .mat file
                mat_data = loadmat('D:\EnEffDaten\Data' + datestring + '.mat', struct_as_record=False, squeeze_me=True)
                matobj=mat_data['Data']
                    # Convert the entire cell array of structs
                self.Data.append( [self.matstruct_to_dict(s) for s in matobj] )
            except ValueError:
                print("The current day is not found in the Data.")
                
       
        
    def matstruct_to_dict(self,matobj):
        """
        Recursively converts a scipy.io.matlab.mio5_params.mat_struct object to a nested dictionary,
        with special handling for MATLAB datetime types.
        """
        if isinstance(matobj, mat_struct):
            result = {}
            for fieldname in matobj._fieldnames:
                element = getattr(matobj, fieldname)
                if isinstance(element, mat_struct):
                    result[fieldname] = self.matstruct_to_dict(element)
                elif isinstance(element, np.ndarray):
                    # Check for datetime conversion
                    if fieldname=='time':
                        result[fieldname] = [datetime.strptime(e, '%Y-%m-%d %H:%M:%S.%f') for e in element]
                        
                    else:
                        result[fieldname] = element
                else:
                    result[fieldname] = element
            return result
        elif isinstance(matobj, np.ndarray):
            if matobj.dtype == 'O':
                return [self.matstruct_to_dict(item) for item in matobj]
            else:
                return matobj
        else:
            return matobj
        
        
    def find_datetime(self,datetime_list, given_datetime, last_found_index):
        # Start from the last found index and move forward
        for i in range(last_found_index, len(datetime_list)):
            if datetime_list[i] >= given_datetime:
                return  0 if i - 1 < 0 else i - 1  # Return the index of the element just before the given datetime
        
        return len(datetime_list) - 1  # If no greater element is found, return the last index
        
    def SupplyVolumeflow(self,datetime):
        # Nord
        AKZ_Nord_supply = 'UFM61UM10F001XQ20'
        try:
            self.SupplyVolumeflow_Nord_supply_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_Nord_supply:
                    self.SupplyVolumeflow_Nord_supply_index = i
                    break
        
        AKZ_Nord_return = 'UFM61UM40F001XQ20'
        try:
            self.SupplyVolumeflow_Nord_return_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_Nord_return:
                    self.SupplyVolumeflow_Nord_return_index = i
                    break
        
        # West
        AKZ_West_supply = 'UFM63UM10F001XQ20'
        try:
            self.SupplyVolumeflow_West_supply_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_West_supply:
                    self.SupplyVolumeflow_West_supply_index = i
                    break
        
        AKZ_West_return = 'UFM63UM40F001XQ20'
        try:
            self.SupplyVolumeflow_West_return_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_West_return:
                    self.SupplyVolumeflow_West_return_index = i
                    break
        
        # Ost
        AKZ_Ost_supply = 'UFM62UM10F001XQ20'
        try:
            self.SupplyVolumeflow_Ost_supply_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_Ost_supply:
                    self.SupplyVolumeflow_Ost_supply_index = i
                    break
        
        AKZ_Ost_return = 'UFM62UM40F001XQ20'
        try:
            self.SupplyVolumeflow_Ost_return_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_Ost_return:
                    self.SupplyVolumeflow_Ost_return_index = i
                    break
                
        # Nord - Supply
        self.idx_Nord_supply_time = self.find_datetime(self.Data[-1][self.SupplyVolumeflow_Nord_supply_index ]['time'],datetime, self.last_idx_Nord_supply_time)
        self.last_idx_Nord_supply_time = self.idx_Nord_supply_time
        self.Qstrom_Nord_supply = self.get_non_zero_value(self.Data[-1][self.SupplyVolumeflow_Nord_supply_index]['value'], self.idx_Nord_supply_time)

        
        # Nord - Return
        self.idx_Nord_return_time = self.find_datetime(self.Data[-1][self.SupplyVolumeflow_Nord_return_index ]['time'],datetime, self.last_idx_Nord_return_time)
        self.last_idx_Nord_return_time = self.idx_Nord_return_time
        self.Qstrom_Nord_return = self.get_non_zero_value(self.Data[-1][self.SupplyVolumeflow_Nord_return_index]['value'], self.idx_Nord_return_time)

        
        # West - Supply
        self.idx_West_supply_time = self.find_datetime(self.Data[-1][self.SupplyVolumeflow_West_supply_index ]['time'],datetime, self.last_idx_West_supply_time)
        self.last_idx_West_supply_time = self.idx_West_supply_time
        self.Qstrom_West_supply = self.get_non_zero_value(self.Data[-1][self.SupplyVolumeflow_West_supply_index]['value'], self.idx_West_supply_time)

        
        # West - Return
        self.idx_West_return_time = self.find_datetime(self.Data[-1][self.SupplyVolumeflow_West_return_index ]['time'], datetime, self.last_idx_West_return_time)
        self.last_idx_West_return_time = self.idx_West_return_time
        self.Qstrom_West_return = self.get_non_zero_value(self.Data[-1][self.SupplyVolumeflow_West_return_index]['value'], self.idx_West_return_time)

        
        # Ost - Supply
        self.idx_Ost_supply_time = self.find_datetime(self.Data[-1][self.SupplyVolumeflow_Ost_supply_index ]['time'],datetime, self.last_idx_Ost_supply_time)
        self.last_idx_Ost_supply_time = self.idx_Ost_supply_time
        self.Qstrom_Ost_supply=self.get_non_zero_value(self.Data[-1][self.SupplyVolumeflow_Ost_supply_index]['value'],self.idx_Ost_supply_time)
        
        # Ost - Return
        self.idx_Ost_return_time = self.find_datetime(self.Data[-1][self.SupplyVolumeflow_Ost_return_index ]['time'],datetime, self.last_idx_Ost_return_time)
        self.last_idx_Ost_return_time = self.idx_Ost_return_time
        self.Qstrom_Ost_return=self.get_non_zero_value(self.Data[-1][self.SupplyVolumeflow_Ost_return_index]['value'],self.idx_Ost_return_time)
        
        self.Q_strom_Nord=0.5*(self.Qstrom_Nord_supply+self.Qstrom_Nord_return)
        self.Q_strom_West=0.5*(self.Qstrom_West_supply+self.Qstrom_West_return)
        self.Q_strom_Ost=0.5*(self.Qstrom_Ost_supply+self.Qstrom_Ost_return)
        
        self.Q_strom_gesamt=self.Q_strom_Nord+self.Q_strom_West+self.Q_strom_Ost
        
        return self.Q_strom_gesamt
        

  
        
    

        
        
    def SupplyTemperature(self,datetime):
        AKZ_Supply_T = 'MFM60UM20T001XQ01'
        try:
            self.Temperature_Supply_index
        except:
            for i, entry in enumerate(self.Data[-1]):
                if entry.get("AKZ") == AKZ_Supply_T:
                    self.Temperature_Supply_index = i
                    break
                
        self.idx_T_supply_time = self.find_datetime(self.Data[-1][self.Temperature_Supply_index]['time'],datetime, self.last_idx_T_supply_time)
        self.last_idx_T_supply_time = self.idx_T_supply_time
        self.T_supply = self.get_non_zero_value(self.Data[-1][self.Temperature_Supply_index]['value'],self.idx_T_supply_time)
        return self.T_supply
    
    def get_non_zero_value(self,data,index):
        """Helper function to get the most recent non-zero value."""
        while index >= 0:
            value = data[index]
            if value != 0:
                return value
            index -= 1
        return 0  # Fallback if no non-zero value is found
    
    def TemperatureOutside(self,file_path):
        self.T_out_datetimes = []
        self.T_out = []
        
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter=';')
            
            # Skip the header
            next(reader)
            
            for row in reader:
                # MESS_DATUM is the second column (index 1)
                mess_datum = row[1].strip()
                
                # TT_TU is the fourth column (index 3)
                tt_tu = float(row[3].strip())
                
                # Convert MESS_DATUM to a datetime object
                # MESS_DATUM is formatted as YYYYMMDDHH
                date_time = datetime.strptime(mess_datum, '%Y%m%d%H')
                
                # Append the datetime object to the list
                self.T_out_datetimes.append(date_time)
                
                # Append the temperature to the list
                self.T_out.append(tt_tu)
                
    def updateTout(self,TimeSeries,ProgramModes):
        #Tsupply=Tsupply+0.1 #here the Supply Tempreture can be changed
        if len(TimeSeries.datetimeStationary)>1:
            if TimeSeries.datetimeStationary[-1].hour!=TimeSeries.datetimeStationary[-2].hour:
                current_hour=TimeSeries.datetimeStationary[-1].replace(minute=0, second=0, microsecond=0)
                try:
                    idx = self.T_out_datetimes.index(current_hour)
                    TimeSeries.T_out=self.T_out[idx]
                except ValueError:
                    print("The current hour is not found in the list.")
        else:
            current_hour=TimeSeries.datetimeStationary[-1].replace(minute=0, second=0, microsecond=0)
            try:
                idx = self.T_out_datetimes.index(current_hour)
                TimeSeries.T_out=self.T_out[idx]
            except ValueError:
                print("The current hour is not found in the list.")
        
                
