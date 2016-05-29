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
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDSampleSignalAnalysis, PKPDSignalAnalysis


class ProtPKPDSABase(ProtPKPD):
    """ Base signal analysis protocol"""

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams1(self, form, addPredictorPredicted=True):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        if addPredictorPredicted:
            form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default="t",
                          help='Y is predicted as an exponential function of X, Y=f(X)')
            form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default="Cp",
                          help='Y is predicted as an exponential function of X, Y=f(X)')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runAnalysis',self.inputExperiment.get().getObjId(), self.getListOfFormDependencies())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getListOfFormDependencies(self):
        return None

    def getXYvars(self):
        pass

    def createAnalysis(self):
        pass

    def setupFromFormParameters(self):
        pass

    def prepareForSampleAnalysis(self, sampleName):
        pass

    def runAnalysis(self, objId, otherDependencies):
        self.getXYvars()
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        self.setupFromFormParameters()
        self.createAnalysis()

        self.signalAnalysis = PKPDSignalAnalysis()
        self.signalAnalysis.fnExperiment.set(self.inputExperiment.get().fnPKPD.get())
        self.signalAnalysis.predictor=self.experiment.variables[self.varNameX]
        self.signalAnalysis.predicted=self.experiment.variables[self.varNameY]
        self.signalAnalysis.analysisDescription=self.analysis.getDescription()
        self.signalAnalysis.analysisParameters = self.analysis.getParameterNames()

        self.printSection("Processing samples")
        for sampleName, sample in self.experiment.samples.iteritems():
            print("%s -------------------"%sampleName)
            if not self.prepareForSampleAnalysis(sampleName):
                continue

            # Actually analyze
            x, y = sample.getXYValues(self.varNameX,self.varNameY)
            self.analysis.setXYValues(x, y)
            self.analysis.calculateParameters(show=True)

            # Keep this result
            sampleAnalysis = PKPDSampleSignalAnalysis()
            sampleAnalysis.sampleName = sampleName
            sampleAnalysis.x = x
            sampleAnalysis.y = y
            sampleAnalysis.analysisVariables = self.analysis.getParameterNames()
            sampleAnalysis.parameters = self.analysis.parameters
            self.signalAnalysis.sampleAnalyses.append(sampleAnalysis)

            print(" ")
        self.signalAnalysis.write(self._getPath("analysis.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputAnalysis=self.signalAnalysis)
        self._defineSourceRelation(self.inputExperiment, self.signalAnalysis)
