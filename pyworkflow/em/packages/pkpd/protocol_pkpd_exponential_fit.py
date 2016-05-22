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

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDDEOptimizer, PKPDLSOptimizer, PKPDFitting, PKPDSampleFit
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.object import String, Integer
from utils import parseRange
from pk_models import PKPDExponentialModel


class ProtPKPDExponentialFit(ProtPKPD):
    """ Fit a set of exponentials. The observed measurement is modelled as Y=sum_{i=1}^N c_i exp(-lambda_i * X).\n
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit exponentials'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default="t",
                      help='Y is predicted as an exponential function of X, Y=sum_{i=1}^N c_i exp(-lambda_i * X)')
        form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default="Cp",
                      help='Y is predicted as an exponential function of X, Y=sum_{i=1}^N c_i exp(-lambda_i * X)')
        if fullForm:
            form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=1,
                          help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                               "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
            form.addParam('Nexp', params.IntParam, label="Number of exponentials", default=1,
                          help='Number of exponentials to fit')
        else:
            self.fitType=Integer()
            self.fitType.set(1)
            self.Nexp=Integer()
            self.Nexp.set(1)
        form.addParam('bounds', params.StringParam, label="Amplitude and time constant bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Bounds for the c_i amplitudes.\nExample 1: (0,10);(0,1e-2) -> c1 in (0,10), lambda1 in (0,1e-2)\n'\
                           'Example 2: (0,10);(0,1e-2);(0,1);(0,1e-1) -> c1 in (0,10), lambda1 in (0,1e-2), c2 in (0,1), lambda2 in (0,1e-1)')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        if fullForm:
            form.addParam('reportTime', params.StringParam, label="Evaluate at X=", default="", expertLevel=LEVEL_ADVANCED,
                          help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:0.55:10, from 0 to 10 in steps of 0.5')
        else:
            self.reportTime=String()
            self.reportTime.set("")

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputExperiment.get().getObjId(), self.predictor.get(), \
                                 self.predicted.get(), self.fitType.get(), self.Nexp.get(), self.bounds.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFit(self, objId, X, Y, fitType, Nexp, bounds):
        reportTime = parseRange(self.reportTime.get())
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        # Setup model
        self.printSection("Model setup")
        model = PKPDExponentialModel()
        model.setExperiment(experiment)
        model.setXVar(self.predictor.get())
        model.setYVar(self.predicted.get())
        model.Nexp=self.Nexp.get()
        model.printSetup()

        # Create output object
        self.fitting = PKPDFitting()
        self.fitting.fnExperiment.set(self.inputExperiment.get().fnPKPD.get())
        self.fitting.predictor=experiment.variables[self.predictor.get()]
        self.fitting.predicted=experiment.variables[self.predicted.get()]
        self.fitting.modelDescription=model.getDescription()
        self.fitting.modelParameters = model.getParameterNames()

        # Actual fitting
        if self.fitType.get()==0:
            fitType = "linear"
        elif self.fitType.get()==1:
            fitType = "log"
        elif self.fitType.get()==2:
            fitType = "relative"

        for sampleName, sample in experiment.samples.iteritems():
            self.printSection("Fitting "+sampleName)
            x, y = sample.getXYValues(self.predictor.get(),self.predicted.get())
            print("X= "+str(x))
            print("Y= "+str(y))
            print(" ")
            model.setBounds(self.bounds.get())
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

            # Keep this result
            sampleFit = PKPDSampleFit()
            sampleFit.sampleName = sample.varName
            sampleFit.x = x
            sampleFit.y = y
            sampleFit.yp = model.yPredicted
            sampleFit.yl = model.yPredictedLower
            sampleFit.yu = model.yPredictedUpper
            sampleFit.parameters = model.parameters
            sampleFit.modelEquation = model.getEquation()
            sampleFit.copyFromOptimizer(optimizer2)
            self.fitting.sampleFits.append(sampleFit)
            self.fitting.write(self._getPath("fitting.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineSourceRelation(self.inputExperiment, self.fitting)

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

    def _citations(self):
        return ['Spiess2010']
