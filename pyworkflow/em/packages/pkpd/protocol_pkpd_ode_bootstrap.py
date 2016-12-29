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

import pyworkflow.protocol.params as params
from pyworkflow.em.data import PKPDFitting, PKPDSampleFitBootstrap, PKPDLSOptimizer
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.biopharmaceutics import DrugSource
from protocol_pkpd_ode_base import ProtPKPDODEBase

# TESTED in test_workflow_gabrielsson_pk02.py

class ProtPKPDODEBootstrap(ProtPKPDODEBase):
    """ Bootstrap of an ODE protocol"""

    _label = 'ODE bootstrap'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDMonoCompartment',
                      # pointerClass='ProtPKPDMonoCompartment, ProtPKPDIVMonoCompartmentUrine, ProtPKPDEV0MonoCompartment, ProtPKPDEV01MonoCompartment, ProtPKPDEV1MonoCompartment',
                      help='Select a run of an ODE model')
        form.addParam('Nbootstrap', params.IntParam, label="Bootstrap samples", default=200, expertLevel=LEVEL_ADVANCED,
                      help='Number of bootstrap realizations for each sample')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('deltaT', params.FloatParam, default=2, label='Step (min)', expertLevel=LEVEL_ADVANCED)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputODE.get().getObjId(), self.Nbootstrap.get(), self.confidenceInterval.get())
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

    def runFit(self, objId, Nbootstrap, confidenceInterval):
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
        if type(self.fitting.predicted)==list:
            self.varNameY = [v.varName for v in self.fitting.predicted]
        else:
            self.varNameY = self.fitting.predicted.varName
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)

        # Create output object
        self.fitting = PKPDFitting("PKPDSampleFitBootstrap")
        self.fitting.fnExperiment.set(self.experiment.fnPKPD.get())
        self.fitting.predictor=self.experiment.variables[self.varNameX]
        if type(self.varNameY)==list:
            self.fitting.predicted=[self.experiment.variables[v] for v in self.varNameY]
        else:
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

            # Output object
            sampleFit = PKPDSampleFitBootstrap()
            sampleFit.sampleName = sample.varName
            sampleFit.parameters = np.zeros((self.Nbootstrap.get(),len(parameters0)),np.double)
            sampleFit.xB = []
            sampleFit.yB = []

            # Bootstrap samples
            if type(x[0])==list:
                idx = [[k for k in range(0,len(xj))] for xj in x]
            else:
                idx = [k for k in range(0,len(x))]
            for n in range(0,self.Nbootstrap.get()):
                if type(x[0])==list:
                    idxB = [sorted(np.random.choice(idxj,len(idxj))) for idxj in idx]
                    xB=[]
                    yB=[]
                    for j in range(len(idxB)):
                        idxBj = idxB[j]
                        xj=x[j]
                        yj=y[j]
                        xB.append([xj[i] for i in idxBj])
                        yB.append([yj[i] for i in idxBj])
                else:
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

                # Evaluate the quality on the whole data set
                self.setXYValues(x, y)
                optimizer2.evaluateQuality()
                print(optimizer2.optimum)
                print("   R2 = %f R2Adj=%f AIC=%f AICc=%f BIC=%f"%(optimizer2.R2,optimizer2.R2adj,optimizer2.AIC,\
                                                                   optimizer2.AICc,optimizer2.BIC))

                # Keep this result
                sampleFit.parameters[n,:] = optimizer2.optimum
                sampleFit.xB.append(str(xB))
                sampleFit.yB.append(str(yB))
                sampleFit.copyFromOptimizer(optimizer2)

            self.fitting.sampleFits.append(sampleFit)

        self.fitting.modelParameters = self.getParameterNames()
        self.fitting.modelDescription = self.getDescription()
        self.fitting.write(self._getPath("bootstrapPopulation.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputPopulation=self.fitting)
        self._defineSourceRelation(self.inputODE.get(), self.fitting)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("Number of bootstrap realizations: %d"%self.Nbootstrap.get())
        msg.append("Confidence interval: %f"%self.confidenceInterval.get())
        return msg

    def _validate(self):
        return []
