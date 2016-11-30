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

class TestGabrielssonPK07Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK07')
        cls.exptFn = cls.dataset.getFile('experiment')
    
    def testGabrielssonPK07Workflow(self):
        # Import an experiment

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

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        print "Fitting a two-compartmentx model (intravenous)..."
        protPKPDIVTwoCompartments = self.newProtocol(ProtPKPDIVTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     predictor='t', predicted='Cp',
                                                     globalSearch=False,
                                                     bounds='(0.0, 0.02); (0.0, 100.0); (0.0, 0.05); (0.0, 100.0)')
        protPKPDIVTwoCompartments.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protPKPDIVTwoCompartments)
        self.assertIsNotNone(protPKPDIVTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDIVTwoCompartments', protPKPDIVTwoCompartments)

if __name__ == "__main__":
    unittest.main()
