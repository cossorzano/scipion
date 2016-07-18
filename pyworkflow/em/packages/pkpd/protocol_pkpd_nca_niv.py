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
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import pyworkflow.protocol.params as params
from protocol_pkpd_sa_base import ProtPKPDSABase
from sa_models import NCANIVModel
from pyworkflow.em.pkpd_units import PKPDUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED

class ProtPKPDNCANIV(ProtPKPDSABase):
    """ Non-compartmental analysis of a non-intravenous bolus.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca non-iv'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protAbsorption', params.PointerParam, label="Input experiment",
                      pointerClass='ProtPKPDAbsorptionRate',
                      help='Select an execution of a protocol estimating the absorption rate.')
        form.addParam('areaCalc', params.EnumParam, choices=["Trapezoidal","Mixed"],
                      label="Method for AUC, AUMC calculation", default=1, expertLevel=LEVEL_ADVANCED,
                      help='The mixed integration uses the trapezoidal is used in the raising side and the log-trapezoidal in the decay side.\n'
                            'See explanation at http://learnpkpd.com/2011/04/02/calculating-auc-linear-and-log-linear\n')

    def getListOfFormDependencies(self):
        return [self.protAbsorption.get().getObjId()]

    #--------------------------- STEPS functions --------------------------------------------
    def getInputExperiment(self):
        return self.protAbsorption.get().outputExperiment

    def setupFromFormParameters(self):
        pass

    def getXYvars(self):
        self.varNameX = self.protAbsorption.get().protElimination.get().predictor.get()
        self.varNameY = self.protAbsorption.get().protElimination.get().predicted.get()

    def createAnalysis(self):
        self.analysis = NCANIVModel()
        self.analysis.setExperiment(self.experiment)
        self.analysis.setXVar(self.varNameX)
        self.analysis.setYVar(self.varNameY)
        self.analysis.F = self.protAbsorption.get().absorptionF.get()
        if self.areaCalc == 0:
            self.analysis.areaCalc = "Trapezoidal"
        elif self.areaCalc == 1:
            self.analysis.areaCalc = "Mixed"

    def prepareForSampleAnalysis(self, sampleName):
        sample = self.experiment.samples[sampleName]
        sample.interpretDose()
        self.analysis.D = sample.getDoseAt(0.0)
        self.analysis.Ke = float(sample.descriptors["Ke"])
        return True

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protAbsorption.get().protElimination.get().predicted.get())
        return msg

    def _validate(self):
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)==0:
            return []
        else:
            return ["This protocol is meant only for intravenous bolus regimens. Check the doses for %s"%(','.join(incorrectList))]

