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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************


import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pyworkflow.em.packages.pkpd import *
from test_workflow import TestWorkflow

class TestGabrielssonPK02Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK02')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK02Workflow(self):
        #First, import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t>50',
                                                  filterType=1, condition='$(t)>50 and $(t)<=300')
        protFilterTime.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Filter concentration variable
        print "Filter Cp"
        protFilterCp = self.newProtocol(ProtPKPDFilterMeasurements,
                                        objLabel='pkpd - filter measurements Cp>0',
                                        filterType=0, condition='$(Cp)<=0')
        protFilterCp.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protFilterCp)
        self.assertIsNotNone(protFilterCp.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('ProtFilterCp', protFilterCp)

        # Fit a single exponential to the input data
        print "Fitting an exponential..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               objLabel='pkpd - elimination rate',
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.assertIsNotNone(protEliminationRate.outputFitting.fnFitting, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEliminationRate.outputExperiment.fnPKPD)
        c1 = float(experiment.samples['Individual1'].descriptors['c1'])
        lambda1 = float(experiment.samples['Individual1'].descriptors['lambda1'])
        self.assertAlmostEqual(c1,4.1394,2)
        self.assertAlmostEqual(lambda1,0.0088,3)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9357,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9278,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-11.3863,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-8.3863,3)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-11.4945,3)

        # Estimate absorption rate
        print "Estimation of the absorption rate..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                              objLabel='pkpd - absorption rate')
        protAbsorptionRate.inputExperiment.set(protFilterCp.outputExperiment)
        protAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(protAbsorptionRate)
        self.assertIsNotNone(protAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the absorption rate estimation ")
        self.assertIsNotNone(protAbsorptionRate.outputFitting.fnFitting, "There was a problem with the absorption rate estimation ")
        self.validateFiles('protAbsorptionRate', protAbsorptionRate)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protAbsorptionRate.outputExperiment.fnPKPD)
        Cmax = float(experiment.samples['Individual1'].descriptors['Cmax'])
        Ka = float(experiment.samples['Individual1'].descriptors['Ka'])
        Ke = float(experiment.samples['Individual1'].descriptors['Ke'])
        Vd = float(experiment.samples['Individual1'].descriptors['Vd'])
        tlag = float(experiment.samples['Individual1'].descriptors['tlag'])
        tmax = float(experiment.samples['Individual1'].descriptors['tmax'])
        self.assertAlmostEqual(Cmax,1.8326,3)
        self.assertAlmostEqual(Ka,0.0356,3)
        self.assertAlmostEqual(Ke,0.0088,3)
        self.assertAlmostEqual(Vd,33.4011,2)
        self.assertAlmostEqual(tlag,12.3937,2)
        self.assertAlmostEqual(tmax,51.8842,1)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDEV1MonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  findtlag=True,  predictor='t',
                                                  predicted='Cp', initType=0, deltaT=0.25)
        protEV1MonoCompartment.inputExperiment.set(protFilterCp.outputExperiment)
        protEV1MonoCompartment.absorptionProtocol.set(protAbsorptionRate)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual1'].descriptors['Cl'])
        V = float(experiment.samples['Individual1'].descriptors['V'])
        Ka = float(experiment.samples['Individual1'].descriptors['Ka'])
        tlag = float(experiment.samples['Individual1'].descriptors['tlag'])
        self.assertAlmostEqual(Cl,0.2912,2)
        self.assertAlmostEqual(V,27.8,0)
        self.assertAlmostEqual(Ka,0.027,2)
        self.assertAlmostEqual(tlag,12.3725,0)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9712,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9615,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-24.9912,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-19.9912,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-23.0516,1)

if __name__ == "__main__":
    unittest.main()
