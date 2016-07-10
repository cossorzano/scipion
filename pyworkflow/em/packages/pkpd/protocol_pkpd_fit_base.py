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

import math
from itertools import izip

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDDEOptimizer, PKPDLSOptimizer, PKPDFitting, PKPDSampleFit, PKPDVariable
from utils import parseRange

class ProtPKPDFitBase(ProtPKPD):
    """ Base fit protocol"""

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams1(self, form, defaultPredictor, defaultPredicted):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default=defaultPredictor,
                      help='Y is predicted as an exponential function of X, Y=f(X)')
        form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default=defaultPredicted,
                      help='Y is predicted as an exponential function of X, Y=f(X)')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputExperiment.get().getObjId(), self.predictor.get(), \
                                 self.predicted.get(), self.fitType.get(), self.bounds.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def createModel(self):
        pass

    def setupFromFormParameters(self):
        pass

    def runFit(self, objId, X, Y, fitType, bounds):
        reportX = parseRange(self.reportX.get())
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        # Setup model
        self.printSection("Model setup")
        self.model = self.createModel()
        self.model.setExperiment(self.experiment)
        self.model.setXVar(self.predictor.get())
        self.model.setYVar(self.predicted.get())
        self.setupFromFormParameters()
        self.model.printSetup()

        # Create output object
        self.fitting = PKPDFitting()
        self.fitting.fnExperiment.set(self.inputExperiment.get().fnPKPD.get())
        self.fitting.predictor=self.experiment.variables[self.predictor.get()]
        self.fitting.predicted=self.experiment.variables[self.predicted.get()]
        self.fitting.modelDescription=self.model.getDescription()
        self.fitting.modelParameters = self.model.getParameterNames()
        self.fitting.modelParameterUnits = None

        # Actual fitting
        if self.fitType.get()==0:
            fitType = "linear"
        elif self.fitType.get()==1:
            fitType = "log"
        elif self.fitType.get()==2:
            fitType = "relative"

        for sampleName, sample in self.experiment.samples.iteritems():
            self.printSection("Fitting "+sampleName)
            x, y = sample.getXYValues(self.predictor.get(),self.predicted.get())
            self.model.calculateParameterUnits(sample)
            if self.fitting.modelParameterUnits==None:
                self.fitting.modelParameterUnits = self.model.parameterUnits
            print("X= "+str(x))
            print("Y= "+str(y))
            print(" ")
            self.model.setBounds(self.bounds.get())
            self.model.setXYValues(x, y)
            self.model.prepare()
            print(" ")

            optimizer1 = PKPDDEOptimizer(self.model,fitType)
            optimizer1.optimize()
            optimizer2 = PKPDLSOptimizer(self.model,fitType)
            optimizer2.optimize()
            optimizer2.setConfidenceInterval(self.confidenceInterval.get())
            if reportX!=None:
                print("Evaluation of the model at specified time points")
                yreportX = self.model.forwardModel(self.model.parameters, reportX)
                print("==========================================")
                print("X     Ypredicted     log10(Ypredicted)")
                print("==========================================")
                for n in range(0,reportX.shape[0]):
                    print("%f %f %f"%(reportX[n],yreportX[n],math.log10(yreportX[n])))
                print(' ')

            # Keep this result
            sampleFit = PKPDSampleFit()
            sampleFit.sampleName = sample.varName
            sampleFit.x = x
            sampleFit.y = y
            sampleFit.yp = self.model.yPredicted
            sampleFit.yl = self.model.yPredictedLower
            sampleFit.yu = self.model.yPredictedUpper
            sampleFit.parameters = self.model.parameters
            sampleFit.modelEquation = self.model.getEquation()
            sampleFit.copyFromOptimizer(optimizer2)
            self.fitting.sampleFits.append(sampleFit)

            # Add the parameters to the sample and experiment
            for varName, varUnits, description, varValue in izip(self.model.getParameterNames(), self.model.parameterUnits, self.model.getParameterDescriptions(), self.model.parameters):
                self.experiment.addParameterToSample(sampleName, varName, varUnits, description, varValue)

        self.fitting.write(self._getPath("fitting.pkpd"))
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.fitting)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=['Predicting %s from %s'%(self.predicted.get(),self.predictor.get())]
        return msg

    def _validate(self):
        errors=[]
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD, False)
        if not self.predictor.get() in experiment.variables:
            errors.append("Cannot find %s as variable"%self.predictor.get())
        if not self.predicted.get() in experiment.variables:
            errors.append("Cannot find %s as variable"%self.predicted.get())
        return errors

    def _citations(self):
        return ['Spiess2010']
