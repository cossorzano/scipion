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

import numpy as np
import math

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDExponentialModel, PKPDDEOptimizer, PKPDLSOptimizer
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtPKPDExponentialFit(ProtPKPD):
    """ Fit a set of exponentials. The observed measurement is modelled as Y=sum_{i=1}^N c_i exp(-lambda_i * X)"""
    _label = 'fit exponentials'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default="t",
                      help='Y is predicted as an exponential function of X, Y=sum_{i=1}^N c_i exp(-lambda_i * X)')
        form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default="Cp",
                      help='Y is predicted as an exponential function of X, Y=sum_{i=1}^N c_i exp(-lambda_i * X)')
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic"], label="Fit mode", default=1,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2')
        form.addParam('Nexp', params.IntParam, label="Number of exponentials", default=1,
                      help='Number of exponentials to fit')
        form.addParam('cBounds', params.StringParam, label="Amplitude bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Bounds for the c_i amplitudes.\nExample 1: (0,10)\nExample 2: (0,10);(1,5)')
        form.addParam('lambdaBounds', params.StringParam, label="Time constant bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Bounds for the lambda_i time constants.\nExample 1: (0,0.01)\nExample 2: (0,1e-2);(0,1e-1)')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('reportTime', params.StringParam, label="Evaluate at X=", default="", expertLevel=LEVEL_ADVANCED,
                      help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: -10:5:10, from -10 to 10 in steps of 5')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputExperiment.get().getObjId(), self.predictor.get(), \
                                 self.predicted.get(), self.fitType.get(), self.Nexp.get(), self.cBounds.get(), \
                                 self.lambdaBounds.get())
        # self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getReportTime(self):
        auxString = self.reportTime.get()
        if auxString=="":
            return None
        elif auxString.startswith('['):
            auxString=auxString.replace('[','')
            auxString=auxString.replace(']','')
            tokens=auxString.split(',')
            auxArray = np.array(tokens, dtype='|S4')
            return auxArray.astype(np.float)
        elif ':' in auxString:
            tokens=auxString.split(':')
            if len(tokens)!=3:
                raise Exception("The X evaluation string is not well formatted: %s"%auxString)
            fromValue = float(tokens[0])
            step= float(tokens[1])
            toValue = float(tokens[2])
            return np.arange(fromValue,toValue,step)

    def runFit(self, objId, X, Y, fitType, Nexp, cBounds, lambdaBounds):
        reportTime = self.getReportTime()
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        # Setup model
        self.printSection("Model setup")
        model = PKPDExponentialModel()
        model.setExperiment(experiment)
        model.setXVar(self.predictor.get())
        model.setYVar(self.predicted.get())
        model.setNexp(self.Nexp.get())
        model.printSetup()

        # Actual fitting
        if self.fitType.get()==0:
            fitType = "linear"
        elif self.fitType.get()==1:
            fitType = "log"

        for sampleName, sample in experiment.samples.iteritems():
            self.printSection("Fitting "+sampleName)
            x, y = sample.getXYValues(self.predictor.get(),self.predicted.get())
            print("X= "+str(x))
            print("Y= "+str(y))
            print(" ")
            model.setBounds(self.cBounds.get(),self.lambdaBounds.get())
            model.setXYValues(x, y)
            model.prepare()
            print(" ")

            optimizer1 = PKPDDEOptimizer(model,fitType)
            optimizer1.optimize()
            optimizer2 = PKPDLSOptimizer(model,fitType)
            optimizer2.optimize()
            optimizer2.setConfidenceInterval(self.confidenceInterval.get())
            if reportTime!=None:
                print("Evaluation of the model at specified time points")
                yReportTime = model.forwardModel(model.parameters, reportTime)
                print("==========================================")
                print("X     Ypredicted     log10(Ypredicted)")
                print("==========================================")
                for n in range(0,reportTime.shape[0]):
                    print("%f %f %f"%(reportTime[n],yReportTime[n],math.log10(yReportTime[n])))
                print(' ')

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        return msg

    def _validate(self):
        errors=[]
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD, False)
        if not self.predictor.get() in experiment.variables:
            errors.append("Cannot find %s as variable"%self.predictor.get())
        if not self.predicted.get() in experiment.variables:
            errors.append("Cannot find %s as variable"%self.predicted.get())
        if self.Nexp.get()<1:
            errors.append("The number of exponentials has to be larger than 0")
        return errors