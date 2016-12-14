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
                                                  deltaT=1.0, findtlag=True,
                                                  predictor='t', predicted='Cp',
                                                  bounds='(0.0, 100.0); (0.0, 0.4); (0.0, 2.0); (0.0, 100.0)')
        protEV0MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protEV0MonoCompartment)
        self.assertIsNotNone(protEV0MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV0MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV0MonoCompartment', protEV0MonoCompartment)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV0MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Rin = float(experiment.samples['Individual'].descriptors['Rin'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        self.assertAlmostEqual(Cl,0.726,1)
        self.assertAlmostEqual(V,96.2128,None,None,2.0)
        self.assertAlmostEqual(Rin,0.0838,2)
        self.assertAlmostEqual(tlag,17.7708,None,None,2.0)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV0MonoCompartment.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9600,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9472,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-21.6269,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-16.6269,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-19.6873,1)

        # Fit a monocompartmental model with mixed zero and first order absorption
        print "Fitting monocompartmental model with zero and first order..."
        protEV01MonoCompartment = self.newProtocol(ProtPKPDEV01MonoCompartment,
                                                   objLabel='pkpd - ev01 monocompartment',
                                                   findtlag=True, predictor='t',
                                                   predicted='Cp', deltaT=1.0,
                                                   bounds='(0.0, 30.0); (0.0, 0.3); (0.0, 300.0); '
                                                          '(0.0, 0.5); (0.0, 1.0); (0.0, 120.0)')
        protEV01MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protEV01MonoCompartment)
        self.assertIsNotNone(protEV01MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV01MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV01MonoCompartment', protEV01MonoCompartment)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV01MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        #Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Rin = float(experiment.samples['Individual'].descriptors['Rin'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        t0 = float(experiment.samples['Individual'].descriptors['t0'])
        self.assertAlmostEqual(Cl,0.7251,1)
        #self.assertAlmostEqual(Ka,0.24314,2)
        self.assertAlmostEqual(V,96.0924,None,None,80.0)
        self.assertAlmostEqual(Rin,0.0837,1)
        self.assertAlmostEqual(tlag,17.4759,None,None,10.0)
        self.assertAlmostEqual(t0,292.7140,None,None,100.0)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV01MonoCompartment.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9600,1)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9296,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-17.6269,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-1.6269,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-14.7175,1)

        # Fit a single exponential to the input data
        print "Compute elimination rate..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.assertIsNotNone(protEliminationRate.outputFitting.fnFitting, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEliminationRate.outputExperiment.fnPKPD)
        c1 = float(experiment.samples['Individual'].descriptors['c1'])
        lambda1 = float(experiment.samples['Individual'].descriptors['lambda1'])
        self.assertAlmostEqual(c1,0.5978,2)
        self.assertAlmostEqual(lambda1,0.0073,3)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9731,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9673,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-13.992,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-9.992,3)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-14.4084,3)

        # Estimate absorption rate
        print "Estimation of the absorption rate..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                             objLabel='pkpd - absorption rate')
        protAbsorptionRate.inputExperiment.set(protChangeCpUnit.outputExperiment)
        protAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(protAbsorptionRate)
        self.assertIsNotNone(protAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the absorption rate estimation ")
        self.assertIsNotNone(protAbsorptionRate.outputFitting.fnFitting, "There was a problem with the absorption rate estimation ")
        self.validateFiles('protAbsorptionRate', protAbsorptionRate)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protAbsorptionRate.outputExperiment.fnPKPD)
        Cmax = float(experiment.samples['Individual'].descriptors['Cmax'])
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        Ke = float(experiment.samples['Individual'].descriptors['Ke'])
        Vd = float(experiment.samples['Individual'].descriptors['Vd'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        tmax = float(experiment.samples['Individual'].descriptors['tmax'])
        self.assertAlmostEqual(Cmax,0.0742,3)
        self.assertAlmostEqual(Ka,0.0091,3)
        self.assertAlmostEqual(Ke,0.0073,3)
        self.assertAlmostEqual(Vd,107.8489,1)
        self.assertAlmostEqual(tlag,23.8888,2)
        self.assertAlmostEqual(tmax,122.1951,1)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protAbsorptionRate.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9348,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9255,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-17.771,2)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-14.771,2)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-16.3162,2)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDEV1MonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  findtlag=True,  predictor='t',
                                                  predicted='Cp', initType=0, deltaT=1.0)
        protEV1MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        protEV1MonoCompartment.absorptionProtocol.set(protAbsorptionRate)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        # Compare model parameter values with Gold standard values
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        self.assertAlmostEqual(Cl,0.7872,2)
        self.assertAlmostEqual(V,96.2126,None,None,2.0)
        self.assertAlmostEqual(Ka,0.0081,2)
        self.assertAlmostEqual(tlag,17.4759,None,None,10.0)
        self.assertAlmostEqual(tlag,24.5894,0)

        # Compare fitting parameter values with Gold standard values
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertAlmostEqual(fitting.sampleFits[0].R2,0.9359,3)
        self.assertAlmostEqual(fitting.sampleFits[0].R2adj,0.9176,3)
        self.assertAlmostEqual(fitting.sampleFits[0].AIC,-15.985,1)
        self.assertAlmostEqual(fitting.sampleFits[0].AICc,-10.985,1)
        self.assertAlmostEqual(fitting.sampleFits[0].BIC,-14.0453,1)

if __name__ == "__main__":
    unittest.main()
