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
Signal Analysis models
"""
import numpy as np
from pyworkflow.em.data import PKPDModelBase

class SAModel(PKPDModelBase):
    def calculateParameters(self, show=True):
        pass

class NCAObsIVModel(SAModel):
    def getDescription(self):
        return "Non-compartmental Analysis based on observations (%s)"%self.__class__.__name__

    def getParameterNames(self):
        return ['AUC_0t','AUC_0inf']

    def calculateParameters(self, show=True):
        t = self.x
        C = self.y
        AUC0t = 0
        for i in range(len(C)-1):
            AUC0t += (t[i+1]-t[i])*(C[i]+C[i+1])
        AUC0t*=0.5
        if show:
            print("AUC(0-t) = %f"%AUC0t)

        self.parameters = []
        self.parameters.append(AUC0t)
        self.parameters.append(0.0)