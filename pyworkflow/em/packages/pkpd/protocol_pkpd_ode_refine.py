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
from itertools import izip

import pyworkflow.protocol.params as params
from pyworkflow.em.data import PKPDFitting, PKPDSampleFit, PKPDLSOptimizer
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.biopharmaceutics import DrugSource
from protocol_pkpd_ode_base import ProtPKPDODEBase

class ProtPKPDODERefine(ProtPKPDODEBase):
    """ Refinement of an ODE protocol. The parameters are reestimated with a finer sampling rate.
    """

    _label = 'ODE refinement'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDMonoCompartmentUrine, ProtPKPDTwoCompartments', help='Select a run of an ODE model')
        form.addParam('deltaT', params.FloatParam, default=0.5, label='Step (min)', expertLevel=LEVEL_ADVANCED)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputODE.get().getObjId(), self.deltaT.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def parseBounds(self, boundsString):
        if boundsString!="" and boundsString!=None:
            tokens=boundsString.split(';')
            if len(tokens)!=self.getNumberOfParameters():
                raise Exception("The number of bound intervals does not match the number of parameters")
            self.boundsList=[]
            for token in tokens:
                values = token.strip().split(',')
                self.boundsList.append((float(values[0][1:]),float(values[1][:-1])))

    def setBounds(self,sample):
        self.parseBounds(self.protODE.bounds.get())
        self.setBoundsFromBoundsList()

    def setBoundsFromBoundsList(self):
        Nbounds = len(self.boundsList)
        Nsource = self.drugSource.getNumberOfParameters()
        Nmodel = self.model.getNumberOfParameters()
        if Nbounds!=Nsource+Nmodel:
            raise "The number of parameters (%d) and bounds (%d) are different"%(Nsource+Nmodel,Nbounds)
        self.boundsSource = self.boundsList[0:Nsource]
        self.boundsPK = self.boundsList[Nsource:]
        self.model.bounds = self.boundsPK

    def getBounds(self):
        return self.boundsList

    def runFit(self, objId, deltaT):
        self.protODE = self.inputODE.get()
        self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting)

        # Create drug source
        self.drugSource = DrugSource()

        # Setup model
        self.model = self.protODE.createModel()
        self.model.deltaT = self.deltaT.get()
        self.model.setExperiment(self.experiment)
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)

        # Create output object
        self.fitting = PKPDFitting()
        self.fitting.fnExperiment.set(self.experiment.fnPKPD.get())
        self.fitting.predictor=self.experiment.variables[self.varNameX]
        self.fitting.predicted=self.experiment.variables[self.varNameY]
        self.fitting.modelParameterUnits = None

        # Actual fitting
        if self.protODE.fitType.get()==0:
            fitType = "linear"
        elif self.protODE.fitType.get()==1:
            fitType = "log"
        elif self.protODE.fitType.get()==2:
            fitType = "relative"

        parameterNames = None
        for sampleName, sample in self.experiment.samples.iteritems():
            self.printSection("Fitting "+sampleName)

            # Get the values to fit
            x, y = sample.getXYValues(self.varNameX,self.varNameY)
            print("X= "+str(x))
            print("Y= "+str(y))

            # Interpret the dose
            self.protODE.varNameX = self.varNameX
            self.protODE.varNameY = self.varNameY
            self.protODE.model = self.model
            self.protODE.setTimeRange(sample)
            sample.interpretDose()

            self.drugSource.setDoses(sample.parsedDoseList, self.model.t0, self.model.tF)
            self.protODE.configureSource(self.drugSource)
            self.model.drugSource = self.drugSource

            # Prepare the model
            self.model.setSample(sample)
            self.calculateParameterUnits(sample)
            if self.fitting.modelParameterUnits==None:
                self.fitting.modelParameterUnits = self.parameterUnits

            # Get the initial parameters
            if parameterNames==None:
                parameterNames = self.getParameterNames()
            parameters0 = []
            for parameterName in parameterNames:
                parameters0.append(float(sample.descriptors[parameterName]))
            print("Initial solution: %s"%str(parameters0))
            print(" ")

            # Set bounds
            self.setBounds(sample)
            self.setXYValues(x, y)
            self.parameters = parameters0

            optimizer2 = PKPDLSOptimizer(self,fitType)
            optimizer2.optimize()
            optimizer2.setConfidenceInterval(self.protODE.confidenceInterval.get())
            self.setParameters(optimizer2.optimum)

            # Keep this result
            sampleFit = PKPDSampleFit()
            sampleFit.sampleName = sample.varName
            sampleFit.x = x
            sampleFit.y = y
            sampleFit.yp = self.yPredicted
            sampleFit.yl = self.yPredictedLower
            sampleFit.yu = self.yPredictedUpper
            sampleFit.parameters = self.parameters
            sampleFit.modelEquation = self.getEquation()
            sampleFit.copyFromOptimizer(optimizer2)
            self.fitting.sampleFits.append(sampleFit)
            if type(sampleFit.y[0])!=list and type(sampleFit.y[0])!=np.ndarray and (type(sampleFit.yp[0])==list or type(sampleFit.yp[0])==np.ndarray):
                sampleFit.yp=sampleFit.yp[0]
            if type(sampleFit.y[0])!=list and type(sampleFit.y[0])!=np.ndarray and (type(sampleFit.yl[0])==list or type(sampleFit.yl[0])==np.ndarray):
                sampleFit.yl=sampleFit.yl[0]
            if type(sampleFit.y[0])!=list and type(sampleFit.y[0])!=np.ndarray and (type(sampleFit.yu[0])==list or type(sampleFit.yu[0])==np.ndarray):
                sampleFit.yu=sampleFit.yu[0]

            # Add the parameters to the sample and experiment
            for varName, varUnits, description, varValue in izip(self.getParameterNames(), self.parameterUnits, self.getParameterDescriptions(), self.parameters):
                self.experiment.addParameterToSample(sampleName, varName, varUnits, description, varValue)

        self.fitting.modelParameters = self.getParameterNames()
        self.fitting.modelDescription = self.getDescription()
        self.fitting.write(self._getPath("fitting.pkpd"))
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputODE.get(), self.fitting)
        self._defineSourceRelation(self.inputODE.get(), self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("New sampling rate: %f"%self.deltaT.get())
        return msg

    def _validate(self):
        return []
