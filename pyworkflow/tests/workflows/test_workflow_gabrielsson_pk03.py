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

class TestGabrielssonPK03Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK03')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK03Workflow(self):
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

        # Change the concentration unit to mg/L
        print "Change Units"
        protChangeCpUnit = self.newProtocol(ProtPKPDChangeUnits,
                                           objLabel='pkpd - change units (Cp to mg/L)',
                                           labelToChange='Cp', newUnitsCategory=4, newUnitsCategoryConc=1)
        protChangeCpUnit.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protChangeCpUnit)
        self.assertIsNotNone(protChangeCpUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeCpUnit)

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements (t>=5h)',
                                                  filterType=1, condition='$(t)>=300')
        protFilterTime.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Fit a monocompartmental model with zero order absorption
        print "Fitting monocompartmental model with zero order..."
        protEV0MonoCompartment = self.newProtocol(ProtPKPDEV0MonoCompartment,
                                                  objLabel='pkpd - ev0 monocompartment',
                                                  findtlag=True, predictor='t', predicted='Cp',
                                                  bounds='(0.0, 100.0); (0.0, 0.4); (0.0, 2.0); (0.0, 100.0)')
        protEV0MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protEV0MonoCompartment)
        self.assertIsNotNone(protEV0MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV0MonoCompartment', protEV0MonoCompartment)

        # Fit a monocompartmental model with mixed zero and first order absorption
        print "Fitting monocompartmental model with zero and first order..."
        protEV01MonoCompartment = self.newProtocol(ProtPKPDEV01MonoCompartment,
                                                  objLabel='pkpd - ev01 monocompartment',
                                                  findtlag=True, predictor='t', predicted='Cp',
                                                  bounds='(0.0, 30.0); (0.0, 0.3); (0.0, 300.0);'
                                                         ' (0.0, 0.5); (0.0, 1.0); (0.0, 120.0)')
        protEV01MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protEV01MonoCompartment)
        self.assertIsNotNone(protEV01MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV01MonoCompartment', protEV01MonoCompartment)

        # Fit a single exponential to the input data
        print "Compute elimination rate..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)

        # Estimate absorption rate
        print "Estimation of the absorption rate..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                             objLabel='pkpd - absorption rate')
        protAbsorptionRate.inputExperiment.set(protChangeCpUnit.outputExperiment)
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
        protEV1MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        protEV1MonoCompartment.absorptionProtocol.set(protAbsorptionRate)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

if __name__ == "__main__":
    unittest.main()
