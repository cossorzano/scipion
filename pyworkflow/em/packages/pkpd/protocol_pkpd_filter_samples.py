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
from pyworkflow.em.data import PKPDExperiment, PKPDVariable


class ProtPKPDFilterSamples(ProtPKPD):
    """ Filter samles """
    _label = 'filter samples'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('filterType', params.EnumParam, choices=["Exclude","Keep","Remove NA"], label="Filter mode", default=0,
                      help='Exclude or keep samples meeting the following condition')
        form.addParam('condition', params.TextParam, label="Condition",
                      help='Example: $(weight)<200 and $(sex)=="female"')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runFilter',self.inputExperiment.get().getObjId(), self.filterType.get(), \
                                 self.condition.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFilter(self, objId, filterType, condition):
        import copy
        experiment = PKPDExperiment()
        experiment.load(self.inputExperiment.get().fnPKPD)
        print("**********************************************************************************************")
        print("Reading %s"%self.inputExperiment.get().fnPKPD)
        print("**********************************************************************************************")
        experiment._printToStream(sys.stdout)
        print("**********************************************************************************************")
        print("Filtering")
        print("**********************************************************************************************")
        if self.filterType.get()==0:
            filterType="exclude"
        else:
            filterType="keep"

        filteredExperiment = PKPDExperiment()
        filteredExperiment.general = copy.copy(experiment.general)
        filteredExperiment.variables = copy.copy(experiment.variables)
        filteredExperiment.samples = {}
        filteredExperiment.doses = {}

        # http://stackoverflow.com/questions/701802/how-do-i-execute-a-string-containing-python-code-in-python
        safe_list = ['descriptors']
        safe_dict = dict([ (k, locals().get(k, None)) for k in safe_list ])
        usedDoses = []
        for sampleKey, sample in experiment.samples.iteritems():
            ok = False
            try:
                conditionPython = copy.copy(condition)
                for key, variable in experiment.variables.iteritems():
                    if key in sample.descriptors:
                        if variable.varType == PKPDVariable.TYPE_NUMERIC:
                            conditionPython = conditionPython.replace("$(%s)"%key,"%f"%float(sample.descriptors[key]))
                        else:
                            conditionPython = conditionPython.replace("$(%s)"%key,"'%s'"%sample.descriptors[key])
                ok=eval(conditionPython, {"__builtins__" : None }, {})
            except:
                pass
            if (ok and filterType=="keep") or (not ok and filterType=="exclude"):
                filteredExperiment.samples[sampleKey] = copy.copy(sample)
                usedDoses.append(sample.doseName)

        if len(usedDoses)>0:
            for doseName in usedDoses:
                filteredExperiment.doses[doseName] = copy.copy(experiment.doses[doseName])

        self.experiment = filteredExperiment
        self.experiment._printToStream(sys.stdout)
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
