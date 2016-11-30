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
from pyworkflow.em.pkpd_units import inverseUnits, divideUnits, unitFromString, PKPDUnit

import math

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
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        sunits = divideUnits(yunits,xunits)
        self.parameterUnits=[yunits, sunits]

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        retval.append(lowerBound[0]>0 or upperBound[0]<0)
        retval.append(lowerBound[1]>1 or upperBound[1]<1)
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
            idx = np.where(xprime>0)[0]
            p = np.polyfit(np.log(xprime[idx]), self.y[idx], 1)
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
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, xunits]

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDSaturated(PDGenericModel):

    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])

        e0 = parameters[0]
        emax = parameters[1]
        eC50 = parameters[2]

        self.yPredicted = e0 + (emax*x / (eC50 + x))

        return self.yPredicted

    def getDescription(self):
        return "Saturated (%s)"%self.__class__.__name__ #no se reconoce

    def prepare(self):
        if self.bounds == None:
            e0 = np.min(self.y)
            emax = np.max(self.y-e0)
            eC50 = 0.5*(np.max(self.x)+np.min(self.y))
            print("First estimate of saturated term: ")
            print("Y=(%f) + ( (%f)*X / ((%f) + X) )" % (e0, emax, eC50))

            self.bounds = []
            self.bounds.append((0.1 * e0, 10 * e0))
            self.bounds.append((0.1 * emax, 10 * emax))
            self.bounds.append((0.1 * eC50, 10 * eC50))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + (emax*X / (eC50 + X))"

    def getEquation(self):
        toPrint = "Y = (%f) + ( (%f)*X / ((%f) + X) )"%(self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'emax', 'eC50']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y = e0 + (emax*X / (eC50 + X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, xunits]

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(True)
        retval.append(True)
        return retval

    def areParametersValid(self, p):
        return True


class PDSigmoid(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        e0 = parameters[0]
        emax = parameters[1]
        eC50 = parameters[2]
        h = parameters[3]
        eC50prime = eC50**h
        xprime = x**h
        self.yPredicted = e0 - ( (emax*(xprime)) / ( (eC50prime) + (xprime)))
        return self.yPredicted

    def getDescription(self):
        return "Sigmoid (%s)"%self.__class__.__name__
        #def prepare(self): #trabajar

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 - ( emax*(X**h) / ( (eC50**h) + (X**h)))"

    def getEquation(self):
        toPrint = "Y = (%f) - ( (%f)*(X**(%f)) / ( ((%f)**(%f)) + (X**(%f))))"%(self.parameters[0], self.parameters[1], self.parameters[3], self.parameters[2], self.parameters[3], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'emax', 'eC50', 'h']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form  Y = e0 - ( emax*(X**h) / ( (eC50**h) + (X**h)) )'] * self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, xunits]

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDGompertz(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]

        d = np.exp(b - (g*x))
        self.yPredicted = a / (np.exp(d))

        return self.yPredicted


    def getDescription(self):
        return "Gompertz (%s)"%self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = a*exp(-exp(b-g*X))"

    def getEquation(self):
        toPrint = "Y = (%f)*exp(-exp((%f)-(%f)*X))"%(self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=a*exp(-exp(b-g*X))']*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits)]

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDLogistic1(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]

        d = np.exp(b - (g * x))

        self.yPredicted = a / (1 + d)

        return self.yPredicted

    def getDescription(self):
        return "Logistic1 (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = a/(1+exp(b-g*X))"

    def getEquation(self):
        toPrint = "Y = (%f)/(1+exp((%f)-(%f)*X))" % (self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=a/(1+exp(b-g*X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDLogistic2(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]

        d = np.exp(b - (g * x))

        self.yPredicted = 1 / (a + d)

        return self.yPredicted

    def getDescription(self):
        return "Logistic2 (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = 1/(a+exp(b-g*X))"

    def getEquation(self):
        toPrint = "Y = 1/((%f)+exp((%f)-(%f)*X))" % (self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=1/(a+exp(b-g*X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDLogistic3(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]

        d = np.exp(-(g * x))

        self.yPredicted = a / (1 + b*d)

        return self.yPredicted

    def getDescription(self):
        return "Logistic3 (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = a/(1+b*exp(-g*X))"

    def getEquation(self):
        toPrint = "Y = (%f)/(1+(%f)*exp(-(%f)*X))" % (self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=a/(1+b*exp(-g*X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDLogistic4(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]

        d = np.exp(-(g * x))

        self.yPredicted = 1 / (a + b*d)

        return self.yPredicted

    def getDescription(self):
        return "Logistic4 (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = 1/(a+b*exp(-g*X))"

    def getEquation(self):
        toPrint = "Y = 1/(%f)+(%f)*exp(-(%f)*X))" % (self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=1/(a+b*exp(-g*X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDRichards(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]
        d = parameters[3]

        p = np.exp(b-(g * x))

        self.yPredicted = a / ((1+p)**(1/d))

        return self.yPredicted

    def getDescription(self):
        return "Richards (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = a/ ((1 + exp(b - g*X))^(1/d))"

    def getEquation(self):
        toPrint = "Y = (%f)/ ((1 + exp((%f) - (%f)*X))^(1/(%f)))" % (self.parameters[0], self.parameters[1], self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y = a/ ((1 + exp(b - g*X))^(1/d))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDMorgan(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        b = parameters[0]
        g = parameters[1]
        a = parameters[2]
        d = parameters[3]

        xprime = x**d

        self.yPredicted = ((b*g) + (a*xprime)) / (g + xprime)

        return self.yPredicted

    def getDescription(self):
        return "Morgan-Mercer-Flodin (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = (b*g + a*(X^d)) / (g + (X^d))"

    def getEquation(self):
        toPrint = "Y = ((%f)*(%f) + (%f)*(X^(%f))) / ((%f) + (X^(%f)))" % (self.parameters[0], self.parameters[1], self.parameters[2], self.parameters[3], self.parameters[1], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['b', 'g', 'a', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=(b*g+a*(X^d))/(g+(X^d))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True



class PDWeibull(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        self.yPredicted = np.zeros(x.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]
        d = parameters[3]

        xprime = x**d

        self.yPredicted = a - (b*(np.exp(- g * xprime)))

        return self.yPredicted

    def getDescription(self):
        return "Wibull (%s)" % self.__class__.__name__

    # def prepare(self):

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = a - b*exp(-g*(X^d))"

    def getEquation(self):
        toPrint = "Y  = (%f) - (%f)*exp(-(%f)*(X^(%f)))" % (self.parameters[0], self.parameters[1], self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=a-b*exp(-g*(X^d))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]  # COSS: Buscar unidades de C

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True
