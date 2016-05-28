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

import sys

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDSampleSignalAnalysis, PKPDSignalAnalysis
from sa_models import NCAObsICModel


class ProtPKPDNCAIVObs(ProtPKPD):
    """ Non-compartmental analysis.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca iv observations'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('protElimination', params.PointerParam, label="Elimination rate",
                      pointerClass='ProtPKPDEliminationRate',
                      help='Select an execution of a protocol estimating the elimination rate')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runAnalysis',self.inputExperiment.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runAnalysis(self, objId):
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        fitting = self.readFitting(self.protElimination.get().outputFitting.fnFitting.get())

        varNameX = self.protElimination.get().predictor.get()
        varNameY = self.protElimination.get().predicted.get()

        self.analysis = NCAObsIVModel()
        self.analysis.setExperiment(experiment)
        self.analysis.setXVar(varNameX)
        self.analysis.setYVar(varNameY)

        self.signalAnalysis = PKPDSignalAnalysis()
        self.signalAnalysis.fnExperiment.set(self.inputExperiment.get().fnPKPD.get())
        self.signalAnalysis.predictor=experiment.variables[varNameX]
        self.signalAnalysis.predicted=experiment.variables[varNameY]
        self.signalAnalysis.analysisDescription=self.analysis.getDescription()
        self.signalAnalysis.analysisParameters = self.analysis.getParameterNames()

        self.printSection("Processing samples")
        for sampleName, sample in experiment.samples.iteritems():
            print("%s -------------------\n"%sampleName)
            sampleFit = fitting.getSampleFit(sampleName)
            if sampleFit == None:
                print("  Cannot process %s because its elimination rate cannot be found\n\n"%sampleName)
                continue
            lambdaz = sampleFit.parameters[1]
            print("Elimination rate = %f"%lambdaz)

            # Actually analyze
            x, y = sample.getXYValues(varNameX,varNameY)
            self.analysis.setXYValues(x, y)
            self.analysis.calculateParameters()

            # Keep this result
            sampleAnalysis = PKPDSampleSignalAnalysis()
            sampleAnalysis.sampleName = sampleName
            sampleAnalysis.x = x
            sampleAnalysis.y = y
            sampleAnalysis.analysisVariables = self.analysis.getParameterNames()
            sampleAnalysis.parameters = self.analysis.parameters
            self.signalAnalysis.sampleAnalyses.append(sampleAnalysis)
        self.signalAnalysis.write(self._getPath("analysis.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputAnalysis=self.signalAnalysis)
        self._defineSourceRelation(self.inputExperiment, self.signalAnalysis)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the variable %s"%self.protElimination.get().predicted.get())
        return msg