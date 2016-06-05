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
from pyworkflow.em.data import PKPDExperiment
from pyworkflow.em.pkpd_units import PKPDUnit


class ProtPKPDCreateLabel(ProtPKPD):
    """ Create label by performing calculations on already existing labels.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'create label'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelToAdd', params.StringParam, label="Label to add", default="",
                      help='Name of the variable to add')
        form.addParam('expression', params.StringParam, label="Expression to calculate", default="",
                      help='For example, to normalize the apparent volume of distribution by the animal weight use $(Vd)/$(weight)')
        form.addParam('units', params.StringParam, label="Units", default="None",
                      help='For example, L/kg')
        form.addParam('comment', params.StringParam, label="Label comment", default="",
                      help='For example, apparent volume of distribution per kilogram')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCreate',self.inputExperiment.get().getObjId(), self.labelToAdd.get(),
                                 self.expression.get(), self.comment.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runCreate(self, objId, labelToAdd, expression, comment):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        labelToAdd = self.labelToAdd.get().replace(' ',"_")
        units = PKPDUnit(self.units.get())
        for sampleName, sample in self.experiment.samples.iteritems():
            varValue = sample.evaluateExpression(self.expression.get())
            self.experiment.addParameterToSample(sampleName, labelToAdd, units.unit, self.comment.get(), varValue)

        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["%s created as %s"%(self.labelToAdd.get(),self.expression.get())]
        return msg
