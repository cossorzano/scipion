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
import random

import pyworkflow.protocol.params as params
from pyworkflow.em.data import PKPDExperiment, PKPDDose, PKPDSample, PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from biopharmaceutics import DrugSource
from protocol_pkpd_ode_base import ProtPKPDODEBase
from pyworkflow.em.pkpd_units import createUnit, multiplyUnits, strUnit
from utils import find_nearest

class ProtPKPDODEBootstrapSimulate(ProtPKPDODEBase):
    """ Simulate a population of ODE parameters obtained by bootstrap"""

    _label = 'ODE bootstrap simulate'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputPopulation', params.PointerParam, label="Input population",
                      pointerClass='PKPDFitting', help='It must be a fitting coming from a bootstrap sample')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDIVMonoCompartment, ProtPKPDEV1MonoCompartment', help='Select a run of an ODE model')
        form.addParam('doses', params.TextParam, label="Doses", default="RepeatedBolus ; repeated_bolus t=0:24:120 d=60; h; mg",
                      help="Structure: [Dose Name] ; [Description] ; [Units] \n"\
                           "The dose name should have no space or special character\n"\
                           "Valid units are: h, mg, ug, ...\n"\
                           "The description is either a bolus or an infusion as shown in the examples\n"\
                           "\nIt is important that there are two semicolons.\n"\
                           "Examples:\n"\
                           "Infusion0 ; infusion t=0.500000...0.750000 d=60; mg\n"\
                           "Bolus0 ; bolus t=0.000000 d=60; min; mg\n"\
                           "RepeatedBolus ; repeated_bolus t=0:24:120 d=60; h; mg")
        form.addParam('t0', params.FloatParam, label="Initial time (h)", default=0)
        form.addParam('tF', params.FloatParam, label="Final time (h)", default=24*7)
        form.addParam('Nsimulations', params.IntParam, label="Simulation samples", default=200, expertLevel=LEVEL_ADVANCED,
                      help='Number of simulations')
        form.addParam('addStats', params.BooleanParam, label="Add simulation statistics", default=True, expertLevel=LEVEL_ADVANCED,
                      help="Mean, lower and upper confidence levels are added to the output")
        form.addParam('confidenceLevel', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters', condition="addStats")
        form.addParam('addIndividuals', params.BooleanParam, label="Add individual simulations", default=False, expertLevel=LEVEL_ADVANCED,
                      help="Individual simulations are added to the output")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate',self.inputODE.get().getObjId(), self.Nsimulations.get(),
                                 self.confidenceLevel.get(), self.doses.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, doseName, simulationsX, y):
        newSample = PKPDSample()
        newSample.varName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.doseDictPtr = self.outputExperiment.doses
        newSample.descriptors = {}
        newSample.doseList = [doseName]
        newSample.addMeasurementPattern([self.varNameY])
        newSample.addMeasurementColumn("t",simulationsX)
        newSample.addMeasurementColumn(self.varNameY,y)
        self.outputExperiment.samples[sampleName] = newSample

    def runSimulate(self, objId, Nsimulations, confidenceInterval, doses):
        self.protODE = self.inputODE.get()
        self.findtlag = self.protODE.findtlag
        self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.inputPopulation.get().fnFitting,cls="PKPDSampleFitBootstrap")
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName

        # Create drug source
        self.drugSource = DrugSource()

        # Create output object
        self.outputExperiment = PKPDExperiment()
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit("min")

        self.outputExperiment.variables[self.varNameX] = tvar
        self.outputExperiment.variables[self.varNameY] = self.experiment.variables[self.varNameY]
        self.outputExperiment.general["title"]="Simulated ODE response"
        self.outputExperiment.general["comment"]="Simulated ODE response"

        # Setup model
        self.model = self.protODE.createModel()
        self.model.setExperiment(self.outputExperiment)
        if hasattr(self.protODE,"deltaT"):
            self.model.deltaT = self.protODE.deltaT.get()
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)
        Nsamples = int(60*math.ceil((self.tF.get()-self.t0.get())/self.model.deltaT))+1
        self.model.x = [self.t0.get()+i*self.model.deltaT for i in range(0,Nsamples)]

        # Read the doses
        for line in self.doses.get().split(';;'):
            tokens = line.split(';')
            if len(tokens)!=4:
                print("Skipping dose: ",line)
                continue
            dosename = tokens[0].strip()
            self.outputExperiment.doses[dosename] = PKPDDose()
            self.outputExperiment.doses[dosename].parseTokens(tokens)

        # Check units
        # Dunits = self.outputExperiment.doses[dosename].dunits
        # Cunits = self.experiment.variables[self.varNameY].units

        # Simulate the different responses
        simulationsX = self.model.x
        simulationsY = np.zeros((self.Nsimulations.get(),len(simulationsX)))
        AUCarray = np.zeros(self.Nsimulations.get())
        AUMCarray = np.zeros(self.Nsimulations.get())
        MRTarray = np.zeros(self.Nsimulations.get())
        for i in range(0,self.Nsimulations.get()):
            # Create sample
            sampleName = "Simulation_%d"%i
            newSample = PKPDSample()
            newSample.varName = sampleName
            newSample.variableDictPtr = self.outputExperiment.variables
            newSample.doseDictPtr = self.outputExperiment.doses
            newSample.descriptors = {}
            newSample.doseList = [dosename]
            newSample.addMeasurementPattern([self.varNameY])
            self.setTimeRange(newSample)

            # Take parameters randomly from the population
            nfit = int(random.uniform(0,len(self.fitting.sampleFits)))
            sampleFit = self.fitting.sampleFits[nfit]
            nprm = int(random.uniform(0,sampleFit.parameters.shape[0]))
            parameters = sampleFit.parameters[nprm,:]
            print("Simulated sample %d: %s"%(i,str(parameters)))

            # Prepare source and this object
            newSample.interpretDose()
            self.drugSource.setDoses(newSample.parsedDoseList, self.model.t0, self.model.tF)
            self.protODE.configureSource(self.drugSource)
            self.model.drugSource = self.drugSource
            parameterNames = self.getParameterNames() # Necessary to count the number of source and PK parameters

            # Prepare the model
            self.setParameters(parameters)
            self.model.setSample(newSample)
            y = self.forwardModel(parameters, simulationsX)
            newSample.addMeasurementColumn("t",simulationsX)
            newSample.addMeasurementColumn(self.varNameY,y)

            # Evaluate AUC, AUMC and MRT in the last full period
            tperiod0 = self.drugSource.parsedDoseList[-2].t0
            tperiodF = self.drugSource.parsedDoseList[-1].t0-self.model.deltaT
            idx0 = find_nearest(simulationsX,tperiod0)
            idxF = find_nearest(simulationsX,tperiodF)
            AUC0t = 0
            AUMC0t = 0
            t = simulationsX
            C = y
            for idx in range(idx0,idxF+1):
                dt = (t[idx+1]-t[idx])
                if C[idx+1]>=C[idx]: # Trapezoidal in the raise
                    AUC0t  += 0.5*dt*(C[idx]+C[idx+1])
                    AUMC0t += 0.5*dt*(C[idx]*t[idx]+C[idx+1]*t[idx+1])
                else: # Log-trapezoidal in the decay
                    decrement = C[idx]/C[idx+1]
                    K = math.log(decrement)
                    B = K/dt
                    AUC0t  += dt*(C[idx]-C[idx+1])/K
                    AUMC0t += (C[idx]*(t[idx]-tperiod0)-C[idx+1]*(t[idx+1]-tperiod0))/B-(C[idx+1]-C[idx])/(B*B)
            MRT = AUMC0t/AUC0t
            Cunits = self.experiment.variables[self.varNameY].units
            AUCunits = multiplyUnits(tvar.units.unit,Cunits.unit)
            AUMCunits = multiplyUnits(tvar.units.unit,AUCunits)
            print("   AUC0t=%f [%s]"%(AUC0t,strUnit(AUCunits)))
            print("   AUMC0t=%f [%s]"%(AUMC0t,strUnit(AUMCunits)))
            print("   MRT=%f [min]"%MRT)

            # Keep result
            simulationsY[i,:] = y
            AUCarray[i] = AUC0t
            AUMCarray[i] = AUMC0t
            MRTarray[i] = MRT
            if self.addIndividuals:
                self.outputExperiment.samples[sampleName] = newSample
            else:
                newSample = None # Free memory

        # Report NCA statistics
        alpha_2 = (100-self.confidenceLevel.get())/2
        limits = np.percentile(AUCarray,[alpha_2,100-alpha_2])
        print("AUC %f%% confidence interval=[%f,%f] [%s]"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(AUCunits)))
        limits = np.percentile(AUMCarray,[alpha_2,100-alpha_2])
        print("AUMC %f%% confidence interval=[%f,%f] [%s]"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(AUMCunits)))
        limits = np.percentile(MRTarray,[alpha_2,100-alpha_2])
        print("MRT %f%% confidence interval=[%f,%f] [min]"%(self.confidenceLevel.get(),limits[0],limits[1]))

        # Calculate statistics
        if self.addStats:
            mu = np.mean(simulationsY,axis=0)
            limits = np.percentile(simulationsY,[alpha_2,100-alpha_2],axis=0)

            self.addSample("Mean", dosename, simulationsX, mu)
            self.addSample("LowerLimit", dosename, simulationsX, limits[0])
            self.addSample("UpperLimit", dosename, simulationsX, limits[1])

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputODE.get(), self.outputExperiment)
        self._defineSourceRelation(self.inputPopulation.get(), self.outputExperiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("Number of simulations: %d"%self.Nsimulations.get())
        msg.append("Confidence interval: %f"%self.confidenceInterval.get())
        msg.append("Dose: %s"%self.doses.get())
        return msg

    def _validate(self):
        msg=[]
        if not self.inputPopulation.get().fnFitting.get().endswith("bootstrapPopulation.pkpd"):
            msg.append("Population must be a bootstrap sample")
        return msg