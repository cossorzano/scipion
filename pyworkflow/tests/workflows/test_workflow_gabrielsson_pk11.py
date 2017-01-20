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

class TestGabrielssonPK11Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('Gabrielsson_PK11')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK11Workflow(self):
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

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t<=1440',
                                                  filterType=1, condition='$(t)<=1440')
        protFilterTime.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        # print "Import Experiment PO"
        # protImportPO = self.newProtocol(ProtImportExperiment,
        #                               objLabel='pkpd - import experiment',
        #                               inputFile=self.exptFnPO)
        # self.launchProtocol(protImportPO)
        # self.assertIsNotNone(protImportPO.outputExperiment.fnPKPD, "There was a problem with the import")
        # self.validateFiles('protImport', protImportPO)
        #
        # # Fit a two-compartment model with oral absorption to a set of measurements
        # print "Fitting a two-compartment model (oral)..."
        # protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
        #                                              objLabel='pkpd - ev1 two-compartments',
        #                                              globalSearch=False,
        #                                              bounds='(10.0, 20.0); (0.25, 0.45); (0.02, 0.06); (0.9, 1.15); (40.0, 70.0); (1.0, 3.0); (40.0, 70.0)')
        # protPKPDPOTwoCompartments.inputExperiment.set(protImportPO.outputExperiment)
        # self.launchProtocol(protPKPDPOTwoCompartments)
        # self.assertIsNotNone(protPKPDPOTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        # self.assertIsNotNone(protPKPDPOTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        # self.validateFiles('protPKPDPOTwoCompartments', protPKPDPOTwoCompartments)
        # experiment = PKPDExperiment()
        # experiment.load(protPKPDPOTwoCompartments.outputExperiment.fnPKPD)
        # Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        # Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        # V = float(experiment.samples['Individual'].descriptors['V'])
        # Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        # Ka = float(experiment.samples['Individual'].descriptors['Bolus_Ka'])
        # F = float(experiment.samples['Individual'].descriptors['Bolus_bioavailability'])
        # tlag = float(experiment.samples['Individual'].descriptors['Bolus_tlag'])
        # self.assertTrue(Cl>0.95 and Cl<1.06)
        # self.assertTrue(Clp>1 and Clp<2.4)
        # self.assertTrue(V>45)
        # self.assertTrue(Vp>40 and Vp<70)
        # self.assertTrue(Ka>0.018 and Ka<0.05)
        # self.assertTrue(F>0.27 and F<0.38)
        # self.assertTrue(tlag>9 and tlag<20)
        # fitting = PKPDFitting()
        # fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        # self.assertTrue(fitting.sampleFits[0].R2>0.6)


if __name__ == "__main__":
    unittest.main()
