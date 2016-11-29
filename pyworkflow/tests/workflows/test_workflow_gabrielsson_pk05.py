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

        # Calculate statistics of the labels
        print "Calculate statics of the labels"
        protStatisticsLabel = self.newProtocol(ProtPKPDIVMonoCompartmentUrine,
                                               objLabel='pkpd - statistics labels')
        protStatisticsLabel.inputExperiment.set(protIVMonoCompartmentUrine.outputExperiment)
        self.launchProtocol(protStatisticsLabel)
        self.assertIsNotNone(protStatisticsLabel.outputExperiment.fnPKPD, "There was a problem with the statics calculation ")
        self.validateFiles('protStatisticsLabel', protStatisticsLabel)

        # Simulate a generic pharmacodynamic response Y=f(X)
        print "Simulate a generic pharmacodynamic response"
        protSimulateGenericPD = self.newProtocol(ProtPKPDSimulateGenericPD,
                                                 objLabel='pkpd - simulate generic',
                                                 paramValues='7;3')
        protSimulateGenericPD.inputExperiment.set(protIVMonoCompartmentUrine.outputExperiment)
        self.launchProtocol(protSimulateGenericPD)
        self.assertIsNotNone(protSimulateGenericPD.outputExperiment.fnPKPD, "There was a problem with the simulation ")
        self.validateFiles('protSimulateGenericPD', protSimulateGenericPD)


if __name__ == "__main__":
    unittest.main()
