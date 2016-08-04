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


class ProtPKPDEV1MonoCompartment(ProtPKPDODEBase):
    """ Fit a monocompartmental model with first order absorption to a set of measurements obtained by oral doses (any arbitrary dosing regimen is allowed)\n
        The differential equation is D(t)=Amax*(1-exp(-Ka*t)) and dC/dt = -Cl * C + 1/V * dD/dt\n
        where C is the concentration, Cl the clearance, V the distribution volume, and D the input dosing regime. Ka is the absorbption rate and
        Amax the dose (maximum amount).
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'ev1 monocompartment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form, True, "t", "Cp")
        form.addParam('initType', params.EnumParam, choices=["From absorption protocol","Specify bounds"], label="Initialize from", default=0,
                      help='You may initialize the model from estimates of a previous protocol or by giving the bounds for the different parameters')
        form.addParam('absorptionProtocol', params.PointerParam, label="Absorption rate protocol",
                      pointerClass='ProtPKPDAbsorptionRate', allowsNull = True, condition="initType==0",
                      help='Absorption rate Protocol')
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], Ka, Cl, V)", default="", condition="initType==1",
                      help="Bounds for the tlag (if it must be estimated), absorption, clearance and volume. Example: (0.01,0.04);(0.2,0.4);(10,20). "\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'Be careful that Cl bounds must be given here. If you have an estimate of the elimination rate, this is Ke=Cl/V. Consequently, Cl=Ke*V ')

    def setupFromFormParameters(self):
        if self.initType==0:
            self.absProt = self.absorptionProtocol.get()
            self.absExperiment = self.readExperiment(self.absProt.outputExperiment.fnPKPD)
        else:
            self.ncaProt = None

    def setBounds(self,sample):
        if self.initType == 0:
            boundsString = ""
            if sample.getSampleName() in self.absExperiment.samples:
                sampleAbs = self.absExperiment.samples[sample.getSampleName()]
                Ka = float(sampleAbs.descriptors["Ka"])
                Ke = float(sampleAbs.descriptors["Ke"])
                V = float(sampleAbs.descriptors["Vd"])
                Cl = Ke*V
                boundsString = "(%f,%f);(%f,%f);(%f,%f)"%(Ka/3,Ka*3,Cl/3,Cl*3,V/3,V*3)
                if self.findtlag:
                    tlag = float(sampleAbs.descriptors["tlag"])
                    boundsString = "(%f,%f);%s"%(tlag/3,tlag*3,boundsString)
            self.parseBounds(boundsString)
            self.setBoundsFromBoundsList()
        else:
            ProtPKPDODEBase.setBounds(self,sample)

    def configureSource(self):
        self.drugSource.type = biopharmaceutics.DrugSource.EV
        self.drugSource.evProfile = biopharmaceutics.BiopharmaceuticsModelOrder1()

    def createModel(self):
        return PK_Monocompartment()

    def _validate(self):
        self.getXYvars()
        errors=[]
        if self.varNameX!=None:
            experiment = self.readExperiment(self.getInputExperiment().fnPKPD, False)
            if not self.varNameX in experiment.variables:
                errors.append("Cannot find %s as variable"%self.varNameX)
            if not self.varNameY in experiment.variables:
                errors.append("Cannot find %s as variable"%self.varNameY)
        return errors
