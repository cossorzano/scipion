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
from protocol_pkpd_sa_base import ProtPKPDSABase
from sa_models import NCAObsIVModel


class ProtPKPDNCAIVObs(ProtPKPDSABase):
    """ Non-compartmental analysis based on observations.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca iv observations'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protElimination', params.PointerParam, label="Elimination rate",
                      pointerClass='ProtPKPDEliminationRate',
                      help='Select an execution of a protocol estimating the elimination rate')
        form.addParam("absorptionF", params.FloatParam, label="Absorption fraction", default=1,
                      help="Between 0 (=no absorption) and 1 (=full absorption)")

    def getInputExperiment(self):
        return self.protElimination.get().outputExperiment

    def getListOfFormDependencies(self):
        return [self.protElimination.get().getObjId()]

    #--------------------------- STEPS functions --------------------------------------------
    def setupFromFormParameters(self):
        self.fitting = self.readFitting(self.protElimination.get().outputFitting.fnFitting.get())

    def getXYvars(self):
        self.varNameX = self.protElimination.get().predictor.get()
        self.varNameY = self.protElimination.get().predicted.get()

    def createAnalysis(self):
        self.analysis = NCAObsIVModel()
        self.analysis.setExperiment(self.experiment)
        self.analysis.setXVar(self.varNameX)
        self.analysis.setYVar(self.varNameY)
        self.analysis.F = self.absorptionF.get()

    def prepareForSampleAnalysis(self, sampleName):
        sampleFit = self.fitting.getSampleFit(sampleName)
        sample = self.experiment.samples[sampleName]
        sample.interpretDose()
        self.analysis.D = sample.getDoseAt(0.0)
        if sampleFit == None:
            print("  Cannot process %s because its elimination rate cannot be found\n\n"%sampleName)
            return False
        self.analysis.lambdaz = sampleFit.parameters[1]
        print("Elimination rate = %f"%self.analysis.lambdaz)
        return True

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protElimination.get().predicted.get())
        return msg

    def _validate(self):
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = []
        for sampleName, sample in experiment.samples.iteritems():
            sample.interpretDose()
            if not sample.isDoseABolus():
                incorrectList.append(sampleName)
        if len(incorrectList)==0:
            return []
        else:
            return ["This protocol is meant only for intravenous bolus regimens. Check the doses for %s"%(','.join(incorrectList))]

