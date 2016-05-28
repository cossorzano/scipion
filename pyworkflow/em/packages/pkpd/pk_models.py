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
PD models
"""
import numpy as np
from pyworkflow.em.data import PKPDModel
import math

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
        print("Model: Y=sum_i c_i*exp(-lambda_i * X)")
        print("Number of exponentials: "+str(self.Nexp))
        print("Bounds: "+str(self.bounds))

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
