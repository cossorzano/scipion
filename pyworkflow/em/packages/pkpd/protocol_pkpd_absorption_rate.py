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

import pyworkflow.protocol.params as params
from protocol_pkpd_fit_base import ProtPKPDFitBase
from pk_models import PKPDSimpleNonIVModel
from pyworkflow.em.pkpd_units import PKPDUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.object import Integer


class ProtPKPDAbsorptionRate(ProtPKPDFitBase):
    """ Estimation of the absorption rate for a non-intravenous route. The estimation is performed after estimating
        the elimination rate. The experiment is determined by the\n
        Protocol created by http://www.kinestatpharma.com\n.
        See the theory at http://www.pharmpress.com/files/docs/Basic%20Pharmacokinetics%20sample.pdf"""
    _label = 'absorption rate'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('protElimination', params.PointerParam, label="Elimination rate",
                      pointerClass='ProtPKPDEliminationRate',
                      help='Select an execution of a protocol estimating the elimination rate')
        form.addParam("absorptionF", params.FloatParam, label="Absorption fraction", default=1,
                      help="Between 0 (=no absorption) and 1 (=full absorption)")

        form.addParam('bounds', params.StringParam, label="Ka, Ke, F bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Bounds for Ka (absorption constant), Ke (elimination constant) and F (bioavailability).\nExample 1: (0,1e-3);(0,1e-2);(0.8,1) -> Ka in (0,1e-3), Ke in (0,1e-2), and F in (0.8,1)\n')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        self.fitType=Integer() # Logarithmic fit
        self.fitType.set(1)

    def getListOfFormDependencies(self):
        return [self.protElimination.get().getObjId(), self.absorptionF.get(), self.bounds.get(), self.confidenceInterval.get()]

    #--------------------------- STEPS functions --------------------------------------------
    def setupFromFormParameters(self):
        self.eliminationFitting = self.readFitting(self.protElimination.get().outputFitting.fnFitting.get())
        self.model.F = self.absorptionF.get()

    def getXYvars(self):
        self.varNameX = self.protElimination.get().predictor.get()
        self.varNameY = self.protElimination.get().predicted.get()

    def createModel(self):
        return PKPDSimpleNonIVModel()

    def prepareForSampleAnalysis(self, sampleName):
        # Keep only the ascending part of the curve
        idx = []
        for i in range(1,self.model.y.shape[0]):
            idx.append(i-1)
            if self.model.y[i-1]>self.model.y[i]:
                break
        self.model.setXYValues(self.model.x[idx],self.model.y[idx])

        sample = self.experiment.samples[sampleName]
        sampleFit = self.eliminationFitting.getSampleFit(sampleName)
        sample.interpretDose()
        self.model.D = sample.getDoseAt(0.0)
        self.model.Dunits = sample.getDoseUnits()
        if sampleFit == None:
            print("  Cannot process %s because its elimination rate cannot be found\n\n"%sampleName)
            return False
        self.model.C0 = sampleFit.parameters[0]
        self.model.C0units = PKPDUnit()
        self.model.C0units.unit = self.eliminationFitting.modelParameterUnits[0]
        self.model.Ke = sampleFit.parameters[1]
        self.model.KeUnits = PKPDUnit()
        self.model.KeUnits.unit = self.eliminationFitting.modelParameterUnits[1]
        print("Concentration at t=0 = %f [%s]"%(self.model.C0,self.model.C0units._toString()))
        print("Elimination rate = %f [%s]"%(self.model.Ke,self.model.KeUnits._toString()))

        return True

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protElimination.get().predicted.get())
        return msg

    def _validate(self):
        msg = []
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)!=0:
            msg.append("This protocol is meant only for bolus regimens. Check the doses for %s"%(','.join(incorrectList)))

        return msg

