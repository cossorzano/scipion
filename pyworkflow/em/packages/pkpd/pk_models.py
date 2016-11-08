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
PK models
"""
import math

import numpy as np

from pyworkflow.em.data import PKPDModel, PKPDODEModel
from pyworkflow.em.pkpd_units import inverseUnits, divideUnits, unitFromString, PKPDUnit


class PKModel(PKPDModel):
    pass

class PKGenericModel(PKModel):
    pass

class PKPDExponentialModel(PKGenericModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        self.yPredicted = np.zeros(x.shape[0])
        proceed=True
        for k in range(1,self.Nexp):
            ck = parameters[2*k]
            ck_1 = parameters[2*(k-1)]
            if ck_1<ck:
                proceed=False
                self.yPredicted = -1000*np.ones(x.shape[0])
                break
        if proceed:
            for k in range(0,self.Nexp):
                ck = parameters[2*k]
                lk = parameters[2*k+1]
                self.yPredicted += ck*np.exp(-lk*x)
        return self.yPredicted

    def getDescription(self):
        return "Sum of exponentials (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            p = np.polyfit(self.x,self.ylog,1)
            print("First estimate of 1 exponential term: ")
            print("Y=%f*exp(-%f*X)"%(math.exp(p[1]),-p[0]))

            cBound = (math.exp(p[1])*0.01,math.exp(p[1])*100.0)
            lambdaBound = (-p[0]*0.01,-p[0]*100.0)
            self.bounds = []
            for i in range(self.Nexp):
                self.bounds.append(cBound)
                self.bounds.append(lambdaBound)

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Number of exponentials: "+str(self.Nexp))
        print("Bounds: "+str(self.bounds))

    def getModelEquation(self):
        return "Y=sum_i c_i*exp(-lambda_i * X)"

    def getEquation(self):
        toPrint="Y="
        for i in range(self.Nexp):
            toPrint+= "+[%f*exp(-%f*X)]"%(self.parameters[2*i],self.parameters[2*i+1])
        return toPrint

    def getParameterNames(self):
        parameterList = []
        for i in range(self.Nexp):
            parameterList.append('c%d'%(i+1))
            parameterList.append('lambda%d'%(i+1))
        return parameterList

    def calculateParameterUnits(self,sample):
        xunits = self.experiment.getVarUnits(self.xName)
        yunits = self.experiment.getVarUnits(self.yName)
        self.parameterUnits = [yunits,inverseUnits(xunits)]*self.Nexp
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(self.Nexp):
            cLower = lowerBound[2*i]
            cUpper = upperBound[2*i]
            if cLower<0 and cUpper>0:
                retval.append("False")
            elif cLower>0:
                retval.append("True")
            elif cUpper<0:
                retval.append("Suspicious, this term may be negative")
            else:
                retval.append("NA")

            decayLower = lowerBound[2*i+1]
            decayUpper = upperBound[2*i+1]
            if decayLower<0 and decayUpper>0:
                retval.append("Suspicious, looks like a constant")
            elif decayLower<0:
                retval.append("Suspicious, this term may be unstable")
            elif decayLower>0:
                retval.append("True")
            else:
                retval.append("NA")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0

class PKPDSimpleEVModel(PKModel):
    def __init__(self, includeTlag=True):
        PKModel.__init__(self)
        self.includeTlag = includeTlag

    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x

        Ka=parameters[0]
        Vd=parameters[1]
        if self.includeTlag:
            tlag=parameters[2]
        else:
            tlag = 0.0

        self.yPredicted = np.zeros(len(x))
        for i in range(len(x)):
            if x[i]>=tlag:
                td=x[i]-tlag
                self.yPredicted[i] = Ka*self.F*self.D/(Vd*(Ka-self.Ke))*(np.exp(-self.Ke*td)-np.exp(-Ka*td))
        return self.yPredicted

    def getDescription(self):
        return "Simple non-iv model (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            # Keep only the ascending part of the curve
            idx = []
            for i in range(1,self.y.shape[0]):
                idx.append(i-1)
                if self.y[i-1]>self.y[i]:
                    break
            xAscending = self.x[idx]
            yAscending = self.y[idx]

            ylogE = np.polyval(np.asarray([-self.Ke,math.log(self.C0)],np.double),xAscending)
            ylogToFit = np.log(yAscending)-ylogE
            p = np.polyfit(xAscending,-ylogToFit,1)
            Ka = -p[0]
            Vd = self.D / np.max(self.y) # Incorrect: (Ka*self.F*self.D)/(math.exp(p[1])*(Ka-self.Ke))
            self.bounds=[(Ka*0.2,Ka*5),(0.2*Vd,5*Vd)]
            print("First estimate of Ka: %f"%Ka)
            print("First estimate of Vd: %f"%Vd)
            if self.includeTlag:
                tlag = (p[1]-math.log(self.C0))/(Ka-self.Ke)
                print("First estimate of tlag: %f"%tlag)
                self.bounds.append((-5*abs(tlag),5*abs(tlag)))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: "+str(self.bounds))

    def getModelEquation(self):
        if self.includeTlag:
            return "Y=Ka*F*D/(Vd*(Ka-Ke))*(exp(-Ke*(t-tlag))-exp(-Ka*(t-tlag))"
        else:
            return "Y=Ka*F*D/(Vd*(Ka-Ke))*(exp(-Ke*t)-exp(-Ka*t)"

    def getEquation(self):
        Ka=self.parameters[0]
        Vd=self.parameters[1]
        if self.includeTlag:
            tlag=self.parameters[2]
            return "Y=%f*%f*D/(%f*(%f-%f))*(exp(-%f*(t-(%f)))-exp(-%f*(t-(%f))))"%(Ka,self.F,Vd,Ka,self.Ke,self.Ke,tlag,Ka,tlag)
        else:
            return "Y=%f*%f*D/(%f*(%f-%f))*(exp(-%f*t)-exp(-%f*t))"%(Ka,self.F,Vd,Ka,self.Ke,self.Ke,Ka)

    def getParameterNames(self):
        if self.includeTlag:
            return ['Ka','Vd','tlag']
        else:
            return ['Ka','Vd']

    def calculateParameterUnits(self,sample):
        xunits = self.experiment.getVarUnits(self.xName)
        Vunits = divideUnits(self.Dunits,self.C0units)
        self.parameterUnits = [inverseUnits(xunits),Vunits]
        if self.includeTlag:
            self.parameterUnits.append(xunits)
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Ka
        decayLower = lowerBound[0]
        decayUpper = upperBound[0]
        if decayLower<0 and decayUpper>0:
            retval.append("Suspicious, Ka looks like a constant")
        elif decayLower<0:
            retval.append("Suspicious, Ka may be unstable")
        elif decayLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vd
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vd looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vd seems to be negative")
        else:
            retval.append("True")

        # tlag
        if self.includeTlag:
            tLower = lowerBound[2]
            tUpper = upperBound[2]
            if tLower<0 and tUpper>0:
                retval.append("Suspicious, tlag looks like 0")
            elif tUpper<0:
                retval.append("Suspicious, tlag seems to be negative")
            else:
                retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_Monocompartment(PKPDODEModel):
    def F(self, t, y):
        Cl=self.parameters[0]
        V=self.parameters[1]
        return -Cl/V*y

    def G(self, t, dD):
        V=self.parameters[1]
        return dD/V

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 1

    def getDescription(self):
        return "Monocompartmental model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C + 1/V * dD/dt"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        return "dC/dt = -(%f)/(%f) * C + 1/(%f) dD/dt"%(Cl,V,V)

    def getParameterNames(self):
        return ['Cl','V']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Cl looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Cl may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")
        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_Twocompartments(PKPDODEModel):
    def F(self, t, y):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        C=y[0]
        Cp=y[1]

        Q12 = Clp * (C-Cp)
        return np.array([-(Cl*C - Q12)/V, Q12/Vp],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0],np.double)

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 2

    def getDescription(self):
        return "Two-compartments model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C - Clp/V * (C-Cp) + 1/V * dD/dt and dCp/dt = Clp/Vp * (C-Cp)"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        return "dC/dt = -(%f)/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt; dCp/dt = (%f)/(%f) * (C-Cp)"%(Cl,V,Clp,V,V,Clp,Vp)

    def getParameterNames(self):
        return ['Cl','V','Clp','Vp']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Cl looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Cl may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        # Clp
        ClLower = lowerBound[2]
        ClUpper = upperBound[2]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clp looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vp
        VLower = lowerBound[3]
        VUpper = upperBound[3]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vp looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vp seems to be negative")
        else:
            retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_MonocompartmentUrine(PKPDODEModel):
    def F(self, t, y):
        C = y[0]
        Cl=self.parameters[0]
        V=self.parameters[1]
        fe=self.parameters[2]
        return np.array([-Cl/V*C, fe*Cl*C],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0],np.double)

    def getResponseDimension(self):
        return 2

    def getStateDimension(self):
        return 2

    def getDescription(self):
        return "Monocompartmental model urine (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C + 1/V * dD/dt and dAu/dt=fe*Cl*C"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        fe=self.parameters[2]
        return "dC/dt = -(%f)/(%f) * C + 1/(%f) dD/dt and dAu/dt = (%f)*(%f)*C"%(Cl,V,V,fe,Cl)

    def getParameterNames(self):
        return ['Cl','V','fe']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Ka looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Ka may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        feLower = lowerBound[2]
        feUpper = upperBound[2]
        if feLower<0 and feUpper>0:
            retval.append("Suspicious, fe looks like 0")
        elif feUpper<0:
            retval.append("Suspicious, fe seems to be negative")
        elif feLower>1:
            retval.append("Suspicious, fe seems to be larger than 1")
        else:
            retval.append("True")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0 and p[2]<1
