# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************

import sys

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment


class ProtPKPDNCA(ProtPKPD):
    """ Non-compartmental analysis """
    _label = 'non-compartmental'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runNCA',self.inputExperiment.get().getObjId())
        # self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runNCA(self, objId):
        experiment = PKPDExperiment()
        experiment.load(self.inputExperiment.get().fnPKPD)
        print("Reading %s"%self.inputExperiment.get().fnPKPD)
        experiment._printToStream(sys.stdout)

    def createOutputStep(self):
        pass
        # self._defineSourceRelation(self.inputClasses, vol)

    #--------------------------- INFO functions --------------------------------------------
