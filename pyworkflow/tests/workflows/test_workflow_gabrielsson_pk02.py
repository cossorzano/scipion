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
        self.validateFiles('protEliminationRate', protEliminationRate)


        # Estimate absorption rate
        print "Estimation of the absorption rate..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                              objLabel='pkpd - absorption rate')
        protAbsorptionRate.inputExperiment.set(protFilterCp.outputExperiment)
        protAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(protAbsorptionRate)
        self.assertIsNotNone(protAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the absorption rate estimation ")
        self.validateFiles('protAbsorptionRate', protAbsorptionRate)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDEV1MonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  findtlag=True,  predictor='t',
                                                  predicted='Cp', initType=0)
        protEV1MonoCompartment.inputExperiment.set(protFilterCp.outputExperiment)
        protEV1MonoCompartment.absorptionProtocol.set(protAbsorptionRate)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

if __name__ == "__main__":
    unittest.main()
