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
        ProtFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t>50',
                                                  filterType=1, condition='$(t)>50 and $(t)<=300')
        ProtFilterTime.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(ProtFilterTime)
        self.assertIsNotNone(ProtFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('ProtFilterTime', ProtFilterTime)

        # Filter concentration variable
        print "Filter Cp"
        ProtFilterCp = self.newProtocol(ProtPKPDFilterMeasurements,
                                        objLabel='pkpd - filter measurements Cp>0',
                                        filterType=0, condition='$(Cp)<=0')
        ProtFilterCp.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(ProtFilterCp)
        self.assertIsNotNone(ProtFilterCp.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('ProtFilterCp', ProtFilterCp)

        # Fit a single exponential to the input data
        print "Fitting an exponential..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(ProtFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)

        # Non-compartmental analysis
        print "Performing Non-compartmental analysis..."
        ProtAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                              absorptionF=0)
        ProtAbsorptionRate.inputExperiment.set(ProtFilterCp)
        #ProtAbsorptionRate.inputExperiment.setExtended('outputExperiment')
        ProtAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(ProtAbsorptionRate)
        self.assertIsNotNone(ProtAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('ProtAbsorptionRate', ProtAbsorptionRate)
        #
        # # Fit a monocompartmental model
        # print "Fitting monocompartmental model..."
        # ProtIVMonoCompartment = self.newProtocol(ProtPKPDIVMonoCompartment,
        #                                          findtlag=False, initType=0)
        # ProtIVMonoCompartment.inputExperiment.set(protChangeUnits.outputExperiment)
        # ProtIVMonoCompartment.ncaProtocol.set(protNCAIVObs.outputAnalysis)
        # self.launchProtocol(ProtIVMonoCompartment)
        # self.assertIsNotNone(ProtIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        # self.validateFiles('ProtIVMonoCompartment', ProtIVMonoCompartment)

if __name__ == "__main__":
    unittest.main()
