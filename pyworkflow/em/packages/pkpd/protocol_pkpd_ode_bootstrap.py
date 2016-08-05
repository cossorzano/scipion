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
import numpy as np

import pyworkflow.protocol.params as params
from pyworkflow.em.data import PKPDFitting, PKPDSampleFitBootstrap, PKPDLSOptimizer
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from biopharmaceutics import DrugSource
from protocol_pkpd_ode_base import ProtPKPDODEBase
from pyworkflow.em.pkpd_units import PKPDUnit

class ProtPKPDODEBootstrap(ProtPKPDODEBase):
    """ Bootstrap of an ODE protocol"""

    _label = 'ODE bootstrap'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDIVMonoCompartment, ProtPKPDEV1MonoCompartment', help='Select a run of an ODE model')
        form.addParam('Nbootstrap', params.IntParam, label="Bootstrap samples", default=200, expertLevel=LEVEL_ADVANCED,
                      help='Number of bootstrap realizations for each sample')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputODE.get().getObjId(), self.Nbootstrap.get(), self.confidenceInterval.get())
        # self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFit(self, objId, Nbootstrap, confidenceInterval):
        self.protODE = self.inputODE.get()
        self.findtlag = self.protODE.findtlag
        self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting)

        # Create drug source
        self.drugSource = DrugSource()

        # Setup model
        self.model = self.protODE.createModel()
        self.model.setExperiment(self.experiment)
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)

        # Create output object
        self.fitting = PKPDFitting("PKPDSampleFitBootstrap")
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

            # Output object
            sampleFit = PKPDSampleFitBootstrap()
            sampleFit.sampleName = sample.varName
            sampleFit.parameters = np.zeros((self.Nbootstrap.get(),len(parameters0)),np.double)

            # Bootstrap samples
            idx = [k for k in range(0,len(x))]
            for n in range(0,self.Nbootstrap.get()):
                idxB = sorted(np.random.choice(idx,len(idx)))
                xB = [x[i] for i in idxB]
                yB = [y[i] for i in idxB]
                print("Bootstrap sample %d"%n)
                print("X= "+str(xB))
                print("Y= "+str(yB))
                self.setXYValues(xB, yB)
                self.parameters = parameters0

                optimizer2 = PKPDLSOptimizer(self,fitType)
                optimizer2.verbose = 0
                optimizer2.optimize()
                optimizer2.evaluateQuality()
                print(optimizer2.optimum)

                # Keep this result
                sampleFit.parameters[n,:] = optimizer2.optimum
                sampleFit.copyFromOptimizer(optimizer2)

            self.fitting.sampleFits.append(sampleFit)

        self.fitting.modelParameters = self.getParameterNames()
        self.fitting.modelDescription = self.getDescription()
        self.fitting.write(self._getPath("fitting.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineSourceRelation(self.getInputExperiment(), self.fitting)

    #--------------------------- INFO functions --------------------------------------------
    # def _summary(self):
    #     msg = []
    #     self.getXYvars()
    #     if self.varNameX!=None:
    #         msg.append('Predicting %s from %s'%(self.varNameX,self.varNameY))
    #     return msg
    #
    def _validate(self):
        return []
