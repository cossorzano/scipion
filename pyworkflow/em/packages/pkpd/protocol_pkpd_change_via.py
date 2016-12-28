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
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.pkpd_units import PKPDUnit

class ProtPKPDChangeVia(ProtPKPD):
    """ Change via of administration\n
        This protocol may also be used to change the bioavailability or the tlag
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'change via'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('doseName', params.StringParam, label="Dose name",
                      help='Name of the dose whose via you want to change')
        form.addParam('doseVia', params.StringParam, label="New via",
                      help='New via of the dose, leave it empty to keep current via. Valid vias are iv, ev0, ev1, ev01, evFractional')
        form.addParam('tlag', params.StringParam, label="New tlag",
                      help='New tlag of the dose, leave it empty to let it free so that it can be optimized by an ODE model')
        form.addParam('bioavailability', params.StringParam, label="New bioavailability",
                      help='New bioavailability of the dose, leave it empty to let it free so that it can be optimized by an ODE model')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runChange',self.inputExperiment.get().getObjId(),self.doseName.get(), self.doseVia.get(),
                                 self.tlag.get(),self.bioavailability.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runChange(self, objId, doseName, doseVia, tlag, bioavailability):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        dose = self.experiment.doses[doseName]
        if doseVia!="":
            dose.via = doseVia

        if not tlag:
            if not 'tlag' in dose.paramsToOptimize:
                dose.paramsToOptimize.append("tlag")
                dose.paramsUnitsToOptimize.append(PKPDUnit.UNIT_TIME_MIN)
        else:
            if 'tlag' in dose.paramsToOptimize:
                idx=dose.paramsToOptimize.index('tlag')
                del dose.paramsToOptimize[idx]
                del dose.paramsUnitsToOptimize[idx]
            dose.tlag=float(tlag)

        if not bioavailability:
            if not 'bioavailability' in dose.paramsToOptimize:
                dose.paramsToOptimize.append("bioavailability")
                dose.paramsUnitsToOptimize.append(PKPDUnit.UNIT_NONE)
        else:
            if 'bioavailability' in dose.paramsToOptimize:
                idx=dose.paramsToOptimize.index('bioavailability')
                del dose.paramsToOptimize[idx]
                del dose.paramsUnitsToOptimize[idx]
            dose.bioavailability=float(bioavailability)
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        if not self.doseName.get() in experiment.doses:
            errors.append("%s is not a dose of the experiment"%self.doseName.get())
        return errors

    def _summary(self):
        msg=[]
        return msg