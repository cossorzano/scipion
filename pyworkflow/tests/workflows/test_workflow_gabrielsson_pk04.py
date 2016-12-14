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

class TestGabrielssonPK04Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK04')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK04Workflow(self):
        #First, import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model with first order..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDEV1MonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  findtlag=True,  predictor='t',
                                                  predicted='Cp', initType=1,
                                                  bounds='(0.0, 60.0); (0.0, 0.005); (0.05, 0.25); (30.0, 100.0)')
        protEV1MonoCompartment.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        self.assertAlmostEqual(Cl,0.1338,2)
        self.assertAlmostEqual(V,96.6806,None,None,50.0)
        self.assertAlmostEqual(Ka,0.003,2)
        self.assertAlmostEqual(tlag,40.0001,1)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.986,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9843,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-102.5974,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-101.5974,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-97.565,1)

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                          objLabel='pkpd - filter measurements (t<=24h)',
                                          filterType=1, condition='$(t)<=1440')
        protFilterTime.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model with first order..."
        protEV1MonoCompartment2 = self.newProtocol(ProtPKPDEV1MonoCompartment,
                                                   objLabel='pkpd - ev1 monocompartment 2',
                                                   findtlag=True,  predictor='t',
                                                   predicted='Cp', initType=1, deltaT=1.0,
                                                   bounds='(0,120); (0.0, 0.01); (0.0, 1.0); (0.0, 200.0)')
        protEV1MonoCompartment2.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment2)
        self.assertIsNotNone(protEV1MonoCompartment2.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment2', protEV1MonoCompartment2)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment2.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        self.assertAlmostEqual(Cl,0.1355,2)
        self.assertAlmostEqual(V,52.1712,1)
        self.assertAlmostEqual(Ka,0.0025,2)
        self.assertAlmostEqual(tlag,38.3506,1)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment2.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9992,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9989,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-76.4385,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-71.4385,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-74.4989,1)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model with first order..."
        protEV1MonoCompartment3 = self.newProtocol(ProtPKPDEV1MonoCompartment,
                                                   objLabel='pkpd - ev1 monocompartment 3',
                                                   findtlag=True,  predictor='t',
                                                   predicted='Cp', initType=1, globalSearch=False,
                                                   bounds='(0.0, 60.0); (0.0, 0.005); (0.05, 0.25); (30.0, 100.0)')
        protEV1MonoCompartment3.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment3)
        self.assertIsNotNone(protEV1MonoCompartment3.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment3', protEV1MonoCompartment3)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment3.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        self.assertAlmostEqual(Cl,0.1353,2)
        self.assertAlmostEqual(V,54.4305,1)
        self.assertAlmostEqual(Ka,0.0024,2)
        self.assertAlmostEqual(tlag,31.1877,1)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment3.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9733,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9705,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-85.8543,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-84.8543,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-80.822,1)

if __name__ == "__main__":
    unittest.main()
