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

class TestGabrielssonPK05Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK05')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK05Workflow(self):
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

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses
        print "Fitting a monocompartmental model (intravenous)..."
        protIVMonoCompartment = self.newProtocol(ProtPKPDIVMonoCompartment,
                                                  objLabel='pkpd - iv monocompartment',
                                                  predictor='t', predicted='Cp',
                                                  initType=1, globalSearch=False,
                                                  bounds='(0.01, 0.1); (0.0, 20.0)')
        protIVMonoCompartment.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVMonoCompartment)
        self.assertIsNotNone(protIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVMonoCompartment', protIVMonoCompartment)

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses and urine
        print "Fitting a monocompartmental model (intravenous doses and urine )..."
        protIVMonoCompartmentUrine = self.newProtocol(ProtPKPDIVMonoCompartmentUrine,
                                                  objLabel='pkpd - iv monocompartment urine',
                                                  bounds='(0,0.04);(5,15);(0,1)')
        protIVMonoCompartmentUrine.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVMonoCompartmentUrine)
        self.assertIsNotNone(protIVMonoCompartmentUrine.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVMonoCompartmentUrine', protIVMonoCompartmentUrine)

        # # Filter time variable
        # print "Filter time"
        # protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
        #                                   objLabel='pkpd - filter measurements (t<=24h)',
        #                                   filterType=1, condition='$(t)<=1440')
        # protFilterTime.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        # self.launchProtocol(protFilterTime)
        # self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        # self.validateFiles('protFilterTime', protFilterTime)
        #
        # # Fit a monocompartmental model with first order absorption
        # print "Fitting monocompartmental model with first order..."
        # protEV1MonoCompartment2 = self.newProtocol(ProtPKPDEV1MonoCompartment,
        #                                            objLabel='pkpd - ev1 monocompartment 2',
        #                                            findtlag=True,  predictor='t',
        #                                            predicted='Cp', initType=1, deltaT=1.0,
        #                                            bounds='(0,120); (0.0, 0.01); (0.0, 1.0); (0.0, 200.0)')
        # protEV1MonoCompartment2.inputExperiment.set(protFilterTime.outputExperiment)
        # self.launchProtocol(protEV1MonoCompartment2)
        # self.assertIsNotNone(protEV1MonoCompartment2.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        # self.validateFiles('protEV1MonoCompartment2', protEV1MonoCompartment2)
        #
        # # Fit a monocompartmental model with first order absorption
        # print "Fitting monocompartmental model with first order..."
        # protEV1MonoCompartment3 = self.newProtocol(ProtPKPDEV1MonoCompartment,
        #                                            objLabel='pkpd - ev1 monocompartment',
        #                                            findtlag=True,  predictor='t',
        #                                            predicted='Cp', initType=1, globalSearch=False,
        #                                            bounds='(0.0, 60.0); (0.0, 0.005); (0.05, 0.25); (30.0, 100.0)')
        # protEV1MonoCompartment3.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        # self.launchProtocol(protEV1MonoCompartment3)
        # self.assertIsNotNone(protEV1MonoCompartment3.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        # self.validateFiles('protEV1MonoCompartment3', protEV1MonoCompartment3)

if __name__ == "__main__":
    unittest.main()
