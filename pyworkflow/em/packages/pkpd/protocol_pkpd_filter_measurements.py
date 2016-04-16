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

import sys

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment


class ProtPKPDFilterMeasurements(ProtPKPD):
    """ Filter measurements """
    _label = 'filter measurements'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('filterType', params.EnumParam, choices=["Exclude","Keep","Remove NA"], label="Filter mode", default=0,
                      help='Exclude or keep measurements meeting the following condition.\n"\
                           "NA values are excluded in keep filters, and kept in exclude filters')
        form.addParam('condition', params.TextParam, label="Condition",
                      help='Example: $(t)<200\n $(Cp)>=1000 and $(Cp)<=2000"')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runFilter',self.inputExperiment.get().getObjId(), self.filterType.get(), \
                                 self.condition.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFilter(self, objId, filterType, condition):
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        self.printSection("Filtering")
        if self.filterType.get()==0:
            filterType="exclude"
        else:
            filterType="keep"
        self.experiment = experiment.filterSamples(self.condition.get(),filterType)
        self.experiment._printToStream(sys.stdout)
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
