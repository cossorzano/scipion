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

    def getNumberOfParameters(self):
        return 2

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
        print("Model: Y=e0+s*X")
        print("Bounds: "+str(self.bounds))

    def getEquation(self):
        toPrint="Y=(%f)+(%f)*X"%(self.parameters[0],self.parameters[1])
        return toPrint

    def getParameterNames(self):
        return ['e0','s']

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        retval.append(lowerBound[0]>0 or upperBound[0]<0)
        retval.append(lowerBound[0]>1 or upperBound[0]<1)
        return retval

    def areParametersValid(self, p):
        return True
