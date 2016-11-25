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
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)


        # Non-compartmental analysis
        print "Performing Non-compartmental analysis..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate)
        protAbsorptionRate.inputExperiment.set(protFilterCp.outputExperiment)
        protAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(protAbsorptionRate)
        self.assertIsNotNone(protAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protAbsorptionRate', protAbsorptionRate)

        # Fit a monocompartmental model
        # print "Fitting monocompartmental model..."
        # protIVMonoCompartment = self.newProtocol(ProtPKPDIVMonoCompartment, findtlag=False, initType=0)
        # protIVMonoCompartment.inputExperiment.set(protFilterCp.outputExperiment.outputExperiment)
        # protIVMonoCompartment.ncaProtocol.set(protNCAIVObs.outputAnalysis)
        # self.launchProtocol(protIVMonoCompartment)
        # self.assertIsNotNone(protIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        # self.validateFiles('protIVMonoCompartment', protIVMonoCompartment)

if __name__ == "__main__":
    unittest.main()
