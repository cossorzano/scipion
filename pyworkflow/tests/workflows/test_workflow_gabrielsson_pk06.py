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

class TestGabrielssonPK06Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK06')
        cls.exptIVFn = cls.dataset.getFile('experiment_iv')
        cls.exptEVFn = cls.dataset.getFile('experiment_ev')
    
    def testGabrielssonPK06Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses)"
        protImportIV = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment iv',
                                      inputFile=self.exptIVFn)
        self.launchProtocol(protImportIV)
        self.assertIsNotNone(protImportIV.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportIV', protImportIV)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImportIV.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses
        print "Fitting a monocompartmental model (intravenous)..."
        protIVMonoCompartment = self.newProtocol(ProtPKPDIVMonoCompartment,
                                                  objLabel='pkpd - iv monocompartment',
                                                  predictor='t', predicted='Cp',
                                                  initType=1, bounds='(0.0, 0.1); (0.0, 300.0)')
        protIVMonoCompartment.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVMonoCompartment)
        self.assertIsNotNone(protIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVMonoCompartment', protIVMonoCompartment)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protIVMonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertAlmostEqual(Cl,0.09653,5)
        self.assertAlmostEqual(V,274.62859,1)
        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protIVMonoCompartment.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.99456,5)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.994097,5)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-47.64876,5)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-46.64876,5)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-46.518865,5)

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses and urine
        print "Fitting a monocompartmental model (intravenous doses and urine )..."
        protIVMonoCompartmentUrine = self.newProtocol(ProtPKPDIVMonoCompartmentUrine,
                                                      objLabel='pkpd - iv monocompartment urine',
                                                      globalSearch=False,
                                                      bounds='(0.0,0.2);(150, 400.0);(0.0,1.0)')
        protIVMonoCompartmentUrine.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVMonoCompartmentUrine)
        self.assertIsNotNone(protIVMonoCompartmentUrine.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protIVMonoCompartmentUrine.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVMonoCompartmentUrine', protIVMonoCompartmentUrine)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protIVMonoCompartmentUrine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        fe = float(experiment.samples['Individual'].descriptors['fe'])
        self.assertAlmostEqual(Cl,0.096526,5)
        self.assertAlmostEqual(V,274.63698,1)
        self.assertAlmostEqual(fe,0.0688,5)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protIVMonoCompartmentUrine.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.99739,2)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.996972,2)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-55.73286,2)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-53.73286,2)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-53.608718,2)

        # Import an experiment (extravascular)
        print "Import Experiment (extravascular)"
        protImportEV = self.newProtocol(ProtImportExperiment,
                                        objLabel='pkpd - import experiment po',
                                        inputFile=self.exptEVFn)
        self.launchProtocol(protImportEV)
        self.assertIsNotNone(protImportEV.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportEV', protImportEV)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit2 = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit2.inputExperiment.set(protImportEV.outputExperiment)
        self.launchProtocol(protChangeTimeUnit2)
        self.assertIsNotNone(protChangeTimeUnit2.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeTimeUnit2', protChangeTimeUnit2)

        # Fit a monocompartmental model to a set of measurements obtained by extravascular doses and urine
        print "Fitting a monocompartmental model (extravascular and urine)..."
        protEV1MonoCompartmentUrine = self.newProtocol(ProtPKPDEV1MonoCompartmentUrine,
                                                       objLabel='pkpd - ev1 monocompartment urine',
                                                       globalSearch=False,
                                                       bounds='(0.0, 0.2); (0.0, 10.0); (170.0, 370.0); (0.0, 0.15)')
        protEV1MonoCompartmentUrine.inputExperiment.set(protChangeTimeUnit2.outputExperiment)
        self.launchProtocol(protEV1MonoCompartmentUrine)
        self.assertIsNotNone(protEV1MonoCompartmentUrine.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartmentUrine.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartmentUrine', protEV1MonoCompartmentUrine)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartmentUrine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        fe = float(experiment.samples['Individual'].descriptors['fe'])
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        self.assertAlmostEqual(Cl,0.0851,3)
        self.assertAlmostEqual(V,256.1935,1)
        self.assertAlmostEqual(fe,0.0883,3)
        self.assertAlmostEqual(Ka,0.003,3)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartmentUrine.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.89445,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.87984,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,7.67788,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,11.67788,3)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,10.51008,3)

if __name__ == "__main__":
    unittest.main()
