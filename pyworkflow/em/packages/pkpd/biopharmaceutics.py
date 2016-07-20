# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (info@kinestat.com)
# *
# * Kinestat Pharma
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************
"""
Biopharmaceutics: Drug sources and how they dissolve
"""
import math
import numpy as np
from pyworkflow.em.pkpd_units import PKPDUnit

class BiopharmaceuticsModel:
    def __init__(self):
        self.parameters = []

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def getDescription(self):
        pass

    def getParameterNames(self):
        pass

    def calculateParameterUnits(self,sample):
        pass

    def setParameters(self, parameters):
        self.parameters = parameters

    def getAg(self,t):
        # Total amount of drug absorbed from time 0 to time t
        return 0.0

class BiopharmaceuticsModelOrder0(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Constant absorption rate']

    def getParameterNames(self):
        return ['K']

    def calculateParameterUnits(self,sample):
        return [PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN]

    def getAg(self,t):
        if t<=0:
            return 0.0
        K = self.parameters[0]
        return K*t

class BiopharmaceuticsModelOrder1(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Initial amount','Absorption rate']

    def getParameterNames(self):
        return ['Amax','Ka']

    def calculateParameterUnits(self,sample):
        return [PKPDUnit.UNIT_WEIGHT_mg,PKPDUnit.UNIT_INVTIME_MIN]

    def getAg(self,t):
        if t<=0:
            return 0.0
        Amax = self.parameters[0]
        Ka = self.parameters[1]
        return Amax*(1-math.exp(-Ka*t))

class BiopharmaceuticsModelOrderFractional(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Initial amount','Constant absorption rate', 'alpha']

    def getParameterNames(self):
        return ['Amax','K','alpha']

    def calculateParameterUnits(self,sample):
        return [PKPDUnit.UNIT_WEIGHT_mg,PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN,PKPDUnit.UNIT_NONE]

    def getAg(self,t):
        if t<=0:
            return 0.0
        Amax = self.parameters[0]
        K = self.parameters[1]
        alpha = self.parameters[2]
        aux = alpha*K*t
        if aux>Amax:
            return Amax
        return Amax-math.pow(math.pow(Amax,alpha)-aux,1.0/alpha)

class DrugSource:
    IV = 0
    EV = 1

    def __init__(self):
        self.type = None
        self.parsedDoseList = []
        self.tlag = 0
        self.evProfile = None

    def getDoseAt(self,t0,dt=0.5):
        doseAmount = 0.0
        for dose in self.parsedDoseList:
            doseAmount += dose.getDoseAt(t0,dt)
        return doseAmount

    def getCumulatedDose(self,t0,tF):
        return self.getDoseAt(t0,tF-t0)

    def getDoseUnits(self):
        return self.parsedDoseList[0].dunits.unit

    def getAmountReleasedAt(self,t0,dt=0.5):
        if self.type == DrugSource.IV:
            return self.getDoseAt(t0,dt)
        else:
            return self.evProfile.getAg(t0+dt)-self.evProfile.getAg(t0)

    def getAmountReleasedBetween(self,t0,tF):
        if self.type == DrugSource.IV:
            return self.getDoseAt(t0,tF-t0)
        else:
            return self.evProfile.getAg(tF)-self.evProfile.getAg(t0)

