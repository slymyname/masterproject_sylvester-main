# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:39:21 2024

@author: sander
"""

import pandas as pd
import numpy as np


class StandardConsumptionProfile:
    def readConsumptioncsv(self):
        """Read the standard consumption profiles from a CSV file."""
        temp = pd.read_csv(r'prm/StandardConsumptionProfiles.csv', header=None, delimiter=';', names=range(25))
        temp = temp.to_numpy()

        SPL = []
        idx_start = []
        for i in range(temp.shape[0]):
            try:
                temp[i, :].astype(float)
                SPL.append(temp[i, :].astype(float))
                SPL[-1] = SPL[-1][~np.isnan(SPL[-1])]
            except ValueError:
                SPL.append(temp[i, :])
                idx_start.append(i)
            
        idx_start.append(temp.shape[0] + 1)
        self.dictSPL = {}
        for i in range(len(idx_start) - 1):
            self.dictSPL[SPL[idx_start[i]][0]] = SPL[idx_start[i] + 1:idx_start[i + 1]][0]

    def readWeathercsv(self):
        """Read weather history from a CSV file."""
        temp = pd.read_csv(r'prm/WeatherHistory.csv', header=None, delimiter=',')
        temp = temp.to_numpy()
        temp = temp[365:365 + 365, 3].astype(float)
        day = np.tile(np.arange(0, 7), 53)
        return temp, day

    def getSigmoidCoeff(self, ConsumerType):
        """Get the sigmoid coefficients for a given consumer type."""
        Coeff = self.dictSPL[ConsumerType + '.SigmoidCoefficients']
        return Coeff[0], Coeff[1], Coeff[2], Coeff[3]

    def calcSigmoid(self, T, A, B, C, D):
        """Calculate the sigmoid function value."""
        h = A / (1 + (B / (T - 40)) ** C) + D
        return h

    def daySigmoid(self, h, ConsumerType, weekday):
        """Adjust the sigmoid value based on the weekday."""
        h_day = h * self.dictSPL[ConsumerType + '.WeekdayFactors'][weekday]
        return h_day

    def calcDenom(self, weather, day, A, B, C, D, ConsumerType):
        """Calculate the denominator for energy calculation."""
        denom = 0
        h = self.getSigmoidCoeff(StandardConsumptionProfile,ConsumerType)
        for i in range(len(weather)):
            h = self.calcSigmoid(StandardConsumptionProfile,weather[i], A, B, C, D)
            h_day = self.daySigmoid(StandardConsumptionProfile,h, ConsumerType, day[i])
            denom += h
        return denom

    def calcEnergy(self, h_day, denom, W):
        """Calculate the energy consumption."""
        E = h_day * W / denom * 1000
        return E

    def calcP(self, E, ConsumerType, hour, weekday):
        """Calculate the power consumption."""
        try:
            P = E * self.dictSPL[ConsumerType + '.DaytimeFactors'][hour] / 100
        except KeyError:
            if weekday == 0:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.Mo'][hour] / 100
            elif weekday == 1:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.Di'][hour] / 100
            elif weekday == 2:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.Mi'][hour] / 100
            elif weekday == 3:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.Do'][hour] / 100
            elif weekday == 4:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.Fr'][hour] / 100
            elif weekday == 5:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.Sa'][hour] / 100
            elif weekday == 6:
                P = E * self.dictSPL[ConsumerType + '.DaytimeFactors.So'][hour] / 100
        return P

    def prepareConsumption(self, ConsumerType):
        """Prepare the consumption data."""
        StandardConsumptionProfile.readConsumptioncsv(StandardConsumptionProfile)
        weather, day = StandardConsumptionProfile.readWeathercsv(StandardConsumptionProfile)
        A, B, C, D = StandardConsumptionProfile.getSigmoidCoeff(StandardConsumptionProfile,ConsumerType)
        denom = StandardConsumptionProfile.calcDenom(StandardConsumptionProfile,weather, day, A, B, C, D, ConsumerType)
        return denom

    def calcConsumption(self, ConsumerType, denom, Q, T, weekday, hour):
        """Calculate the consumption for a given hour and weekday."""
        A, B, C, D = StandardConsumptionProfile.getSigmoidCoeff(StandardConsumptionProfile,ConsumerType)
        h = StandardConsumptionProfile.calcSigmoid(StandardConsumptionProfile,T, A, B, C, D)
        h_day = StandardConsumptionProfile.daySigmoid(StandardConsumptionProfile,h, ConsumerType, weekday)
        E = StandardConsumptionProfile.calcEnergy(StandardConsumptionProfile,h_day, denom, Q)
        P = StandardConsumptionProfile.calcP(StandardConsumptionProfile,E, ConsumerType, hour, weekday)
        return P


# Example usage:
# denom = StandardConsumptionProfile().prepareConsumption('OneFamilyHouse.RoomAndWaterNew')
# P = StandardConsumptionProfile().calcConsumption('OneFamilyHouse.RoomAndWaterNew', denom, 3e8, T, 1, hour)
# 'OneFamilyHouse.RoomAndWaterNew'




    
    

        
