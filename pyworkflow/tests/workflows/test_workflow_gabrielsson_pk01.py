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

class TestGabrielssonPK01Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK01')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK01Workflow(self):
        #First, import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        # Change the concentration unit to mg/L
        print "Change Units"
        protChangeUnits = self.newProtocol(ProtPKPDChangeUnits,
                                           labelToChange='Cp', newUnitsCategory=4, newUnitsCategoryConc=1)
        protChangeUnits.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeUnits)
        self.assertIsNotNone(protChangeUnits.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeUnits)

        # Fit a single exponential to the input data
        print "Compute elimination rate..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protChangeUnits.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)

        # Non-compartmental analysis
        print "Performing Non-compartmental analysis..."
        protNCAIVObs = self.newProtocol(ProtPKPDNCAIVObs)
        protNCAIVObs.inputExperiment.set(protChangeUnits.outputExperiment)
        protNCAIVObs.protElimination.set(protEliminationRate)
        self.launchProtocol(protNCAIVObs)
        self.assertIsNotNone(protNCAIVObs.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protNCAIVObs', protNCAIVObs)

        # Fit a monocompartmental model
        print "Fitting monocompartmental model..."
        protIVMonoCompartment = self.newProtocol(ProtPKPDIVMonoCompartment,
                                                 findtlag=False, initType=0)
        protIVMonoCompartment.inputExperiment.set(protChangeUnits.outputExperiment)
        protIVMonoCompartment.ncaProtocol.set(protNCAIVObs)
        self.launchProtocol(protIVMonoCompartment)
        self.assertIsNotNone(protIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('ProtIVMonoCompartment', protIVMonoCompartment)

if __name__ == "__main__":
    unittest.main()
