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
import math
from pyworkflow.em.data import PKPDModelBase

class SAModel(PKPDModelBase):
    def calculateParameters(self, show=True):
        pass

class NCAObsIVModel(SAModel):
    def getDescription(self):
        return "Non-compartmental Analysis based on observations (%s)"%self.__class__.__name__

    def getParameterNames(self):
        return ['AUC_0t','AUC_0inf','AUMC_0t','AUMC_0inf','MRT','Vd_0t','Vd_0inf','Vss','CL_0t','CL_0inf', 'thalf']

    def calculateParameters(self, show=True):
        t = self.x
        C = self.y

        # AUC0t, AUMC0t
        AUC0t = 0
        AUMC0t = 0
        for i in range(len(C)-1):
            dt = (t[i+1]-t[i])
            AUC0t  += dt*(C[i]+C[i+1])
            AUMC0t += dt*(C[i]*t[i]+C[i+1]*t[i+1])
        AUC0t*=0.5
        AUMC0t*=0.5

        # AUC0inf, AUMC0inf
        AUC0inf = AUC0t+C[-1]/self.lambdaz
        AUMC0inf = AUMC0t+C[-1]*(t[-1]+1/self.lambdaz)/self.lambdaz

        # MRT
        MRT = AUMC0inf/AUC0inf

        # Volumes
        Vd0t = self.F*self.D/(AUC0t*self.lambdaz)
        Vd0inf = self.F*self.D/(AUC0inf*self.lambdaz)
        Vss = self.D*AUMC0inf/(AUC0inf*AUC0inf)

        # Clearances
        CL0t = self.F*self.D/AUC0t
        CL0inf = self.F*self.D/AUC0inf

        # Thalf
        thalf = math.log(2.0)*Vd0inf/CL0inf

        # Finish
        if show:
            print("AUC(0-t) = %f"%AUC0t)
            print("AUC(0-inf) = %f"%AUC0inf)
            print("AUMC(0-t) = %f"%AUMC0t)
            print("AUMC(0-inf) = %f"%AUMC0inf)
            print("MRT = %f"%MRT)
            print("Vd(0-t) = %f"%Vd0t)
            print("Vd(0-inf) = %f"%Vd0inf)
            print("Vss = %f"%Vss)
            print("CL(0-t) = %f"%CL0t)
            print("CL(0-inf) = %f"%CL0inf)
            print("thalf = %f"%thalf)
        self.parameters = []
        self.parameters.append(AUC0t)
        self.parameters.append(AUC0inf)
        self.parameters.append(AUMC0t)
        self.parameters.append(AUMC0inf)
        self.parameters.append(MRT)
        self.parameters.append(Vd0t)
        self.parameters.append(Vd0inf)
        self.parameters.append(Vss)
        self.parameters.append(CL0t)
        self.parameters.append(CL0inf)
        self.parameters.append(thalf)
