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

import math
import numpy as np
from pyworkflow.em.data import PKPDModel
from pyworkflow.em.pkpd_units import PKPDUnit

class PDModel(PKPDModel):
    pass

class PDGenericModel(PDModel):
    pass

class PDLinear(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        self.yPredicted = np.zeros(x.shape[0])
        e0 = parameters[0]
        s = parameters[1]
        self.yPredicted = e0+s*x
        return self.yPredicted

    def getDescription(self):
        return "Linear (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            p = np.polyfit(self.x,self.y,1)
            print("First estimate of linear term: ")
            print("Y=%f+%f*X"%(p[1],p[0]))

            self.bounds = []
            self.bounds.append((0.1*p[1],10*p[1]))
            self.bounds.append((0.1*p[0],10*p[0]))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: "+str(self.bounds))

    def getModelEquation(self):
        return "Y=e0+s*X"

    def getEquation(self):
        toPrint="Y=(%f)+(%f)*X"%(self.parameters[0],self.parameters[1])
        return toPrint

    def getParameterNames(self):
        return ['e0','s']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0+s*X']*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        self.parameterUnits=[PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        retval.append(lowerBound[0]>0 or upperBound[0]<0)
        retval.append(lowerBound[0]>1 or upperBound[0]<1)
        return retval

    def areParametersValid(self, p):
        return True



class PDLogLinear(PDGenericModel):

    # voy a ultilizar la formula E = m*log(C(t) + c0) como  E = m*log(C(t) + 10^(e0/m) )
    # suponemos que m esta dentro de parameters en parameters[0]

    def forwardModel(self, parameters, x=None):
        if x==None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])

        m = parameters[0]
        C0 = parameters[1]
        self.yPredicted = [m*math.log(xi - C0) if np.isfinite(xi) and xi>C0 else float("inf") for xi in x]
        return self.yPredicted

    def getDescription(self):
        return "Log-Linear (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            Cmin=np.min(self.x)
            xprime = self.x - 0.9*Cmin
            p = np.polyfit(np.log(xprime), self.y, 1)
            C0=0.9*Cmin
            m=p[0]
            print("First estimate of log-linear term: ")
            print("Y=(%f)*log(X - (%f))" % (m, C0))

            self.bounds = []
            self.bounds.append((0.1 * m, 10 * m))
            self.bounds.append((-9 * C0, 11 * C0)) # C0 +- 10* C0

    def printSetup(self):
        print ("Model: %s " %self.getModelEquation())
        print ("Bounds:  " + str(self.bounds))

    def getModelEquation(self):
        return "Y = m*log(X - C0)"

    def getEquation(self):
        toPrint = "Y=(%f)*log(X - (%f))" % (self.parameters[0], self.parameters[1])
        return toPrint

    def getParameterNames(self):
        return ['m','C0']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y = m*log(X - C0)'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE] # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        return retval

    def areParametersValid(self, p):
        return True
