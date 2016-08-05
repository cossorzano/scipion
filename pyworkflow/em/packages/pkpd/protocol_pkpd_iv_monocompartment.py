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
# *  All comments concerning this program package may be sent to thes
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import pyworkflow.protocol.params as params
from protocol_pkpd_ode_base import ProtPKPDODEBase
from pk_models import PK_Monocompartment
import biopharmaceutics


class ProtPKPDIVMonoCompartment(ProtPKPDODEBase):
    """ Fit a monocompartmental model to a set of measurements obtained by intravenous doses (any arbitrary dosing regimen is allowed)\n
        The differential equation is dC/dt = -Cl * C + 1/V * dD/dt\n
        where C is the concentration, Cl the clearance, V the distribution volume, and D the input dosing regime.
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'iv monocompartment'

    def __init__(self, **kwargs):
        ProtPKPDODEBase.__init__(self, **kwargs)
        self.ncaProt = None

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form)
        form.addParam('initType', params.EnumParam, choices=["From NCA","Specify bounds"], label="Initialize from", default=0,
                      help='You may initialize the model from NCA estimates or by giving the bounds for the different parameters')
        form.addParam('ncaProtocol', params.PointerParam, label="NCA protocol",
                      pointerClass='ProtPKPDNCAIVExp, ProtPKPDNCAIVObs', allowsNull = True, condition="initType==0",
                      help='NCA Protocol')
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], Cl, V)", default="",
                      help="Bounds for the tlag (if it must be estimated), clearance and volume. Example: (0,2);(0.1,0.2);(10,20). "\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'If tlag must be estimated, its bounds must always be specified')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default="t", condition="initType==1",
                      help='Y is predicted as an exponential function of X, Y=f(X)')
        form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default="Cp", condition="initType==1",
                      help='Y is predicted as an exponential function of X, Y=f(X)')

    def getListOfFormDependencies(self):
        retval = ProtPKPDODEBase.getListOfFormDependencies(self)
        return retval + [self.initType.get(), self.ncaProtocol.get()]

    def setupFromFormParameters(self):
        if self.initType==0:
            self.ncaProt = self.ncaProtocol.get()
            self.ncaExperiment = self.readExperiment(self.ncaProt.outputExperiment.fnPKPD)
        else:
            self.ncaProt = None

    def getXYvars(self):
        if self.initType == 1:
            self.varNameX = self.predictor.get()
            self.varNameY = self.predicted.get()
        else:
            if self.ncaProt != None:
                self.ncaProt.getXYvars()
                self.varNameX = self.ncaProt.varNameX
                self.varNameY = self.ncaProt.varNameY
            else:
                self.varNameX = None
                self.varNameY = None

    def configureSource(self, drugSource):
        drugSource.type = biopharmaceutics.DrugSource.IV

    def createModel(self):
        return PK_Monocompartment()

    def setBounds(self,sample):
        if self.initType == 0:
            boundsString = ""
            if sample.getSampleName() in self.ncaExperiment.samples:
                sampleNCA = self.ncaExperiment.samples[sample.getSampleName()]
                if "CL" in sampleNCA.descriptors:
                    Cl = float(sampleNCA.descriptors["CL"]) # NCA from exponential fitting
                    V  = float(sampleNCA.descriptors["Vd"])
                elif "CL_0inf" in sampleNCA.descriptors:
                    Cl = float(sampleNCA.descriptors["CL_0inf"]) # NCA from observations
                    V  = float(sampleNCA.descriptors["Vd_0inf"])
                boundsString = "(%f,%f);(%f,%f)"%(Cl/3,Cl*3,V/3,V*3)
                if self.findtlag:
                    boundsString = self.bounds.get()+";"+boundsString
            self.parseBounds(boundsString)
        else:
            self.parseBounds(self.bounds.get())
        ProtPKPDODEBase.setBounds(self)
        self.model.bounds = self.boundsPK
