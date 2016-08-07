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
import copy
import math
import numpy as np
from pyworkflow.em.pkpd_units import PKPDUnit
from pyworkflow.em.data import PKPDDose

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

    def getEquation(self):
        return ""

    def getModelEquation(self):
        return ""

    def getDescription(self):
        return ""

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(len(self.parameters)):
            lower = lowerBound[i]
            upper = upperBound[i]
            if lower<0 and upper>0:
                retval.append("False")
            elif lower>0 or upper<0:
                retval.append("True")
            else:
                retval.append("NA")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0

class BiopharmaceuticsModelOrder0(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Constant absorption rate']

    def getParameterNames(self):
        return ['K']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        K = self.parameters[0]
        return K*t

    def getEquation(self):
        K = self.parameters[0]
        return "D(t)=(%f)*t"%K

    def getModelEquation(self):
        return "D(t)=K*t"

    def getDescription(self):
        return "Zero order absorption (%s)"%self.__class__.__name__

    def areParametersSignificant(self, lowerBound, upperBound):
        return True

class BiopharmaceuticsModelOrder1(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Absorption rate']

    def getParameterNames(self):
        return ['Ka']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_INVTIME_MIN]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        Ka = self.parameters[0]
        return self.Amax*math.exp(-Ka*t)

    def getEquation(self):
        Ka = self.parameters[0]
        return "D(t)=(%f)*(1-exp(-(%f)*t)"%(self.Amax,Ka)

    def getModelEquation(self):
        return "D(t)=Amax*(1-exp(-Ka*t))"

    def getDescription(self):
        return "First order absorption (%s)"%self.__class__.__name__

class BiopharmaceuticsModelOrderFractional(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Initial amount','Constant absorption rate', 'alpha']

    def getParameterNames(self):
        return ['Amax','K','alpha']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_WEIGHT_mg,PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN,PKPDUnit.UNIT_NONE]
        return self.paramterUnits

    def getAg(self,t):
        # COSS Hay que pensar si esto es correcto
        if t<=0:
            return 0.0
        Amax = self.parameters[0]
        K = self.parameters[1]
        alpha = self.parameters[2]
        aux = alpha*K*t
        if aux>Amax:
            return Amax
        return Amax-math.pow(math.pow(Amax,alpha)-aux,1.0/alpha)

    def getEquation(self):
        Amax = self.parameters[0]
        K = self.parameters[1]
        alpha = self.parameters[2]
        return "D(t)=(%f)-((%f)^(%f)-(%f)*(%f)*t)^(1/(%f))"%(Amax,Amax,alpha,K,alpha,alpha)

    def getModelEquation(self):
        return "D(t)=Amax-(Amax^alpha-alpha*K*t)^(1/alpha)"

    def getDescription(self):
        return "Fractional order absorption (%s)"%self.__class__.__name__

    def areParametersValid(self, p):
        return np.sum(p<0)==0 and p[2]>0 and p[2]<1

class DrugSource:
    IV = 0
    EV = 1

    def __init__(self):
        self.type = None
        self.parsedDoseList = []
        self.tlag = 0
        self.evProfile = None

    def setDoses(self, parsedDoseList, t0, tF):
        self.parsedDoseList = []
        for dose in parsedDoseList:
            if dose.doseType != PKPDDose.TYPE_REPEATED_BOLUS:
                self.parsedDoseList.append(dose)
            else:
                for t in np.arange(dose.t0,dose.tF,dose.every):
                    if t0<=t and t<=tF:
                        newDose = copy.copy(dose)
                        newDose.doseType = PKPDDose.TYPE_BOLUS
                        newDose.t0 = t
                        self.parsedDoseList.append(newDose)

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
        doseAmount = 0.0
        for dose in self.parsedDoseList:
            if self.type == DrugSource.IV:
                doseAmount += dose.getDoseAt(t0-dose.t0-self.tlag,dt)
            else:
                if dose.doseType!=PKPDDose.TYPE_INFUSION:
                    self.evProfile.Amax = dose.doseAmount
                    doseAmount += self.evProfile.getAg(t0-dose.t0-self.tlag)-self.evProfile.getAg(t0-dose.t0-self.tlag+dt)
                else:
                    raise Exception("getAmountReleasedAt not implemented for infusion")
        if doseAmount<0:
            doseAmount=0
        return doseAmount

    def getAmountReleasedBetween(self,t0,tF):
        return self.getAmountReleasedAt(t0,tF-t0)

    def getEquation(self):
        if self.type == DrugSource.IV:
            retval = "D=Div(t)"
        else:
            retval = self.evProfile.getEquation()
        return "%s (tlag=%f)"%(retval,self.tlag)

    def getModelEquation(self):
        if self.type == DrugSource.IV:
            return "D=Div(t)"
        else:
            return self.evProfile.getModelEquation()

    def getDescription(self):
        if self.type == DrugSource.IV:
            return "Intravenous dose"
        else:
            return self.evProfile.getDescription()

    def getParameterNames(self):
        if self.type == DrugSource.IV:
            return []
        else:
            return self.evProfile.getParameterNames()

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def calculateParameterUnits(self,sample):
        if self.type == DrugSource.IV:
            return []
        else:
            return self.evProfile.calculateParameterUnits(sample)

    def areParametersSignificant(self, lowerBound, upperBound):
        if self.type == DrugSource.IV:
            return True
        else:
            return self.evProfile.areParametersSignificant(lowerBound, upperBound)

    def areParametersValid(self, p):
        if self.type == DrugSource.IV:
            return True
        else:
            return self.evProfile.areParametersValid(p)
