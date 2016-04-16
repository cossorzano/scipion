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
from pyworkflow.em.data import PKPDExperiment, PKPDSample


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
        form.addParam('condition', params.TextParam, label="Condition", condition="filterType!=2",
                      help='Example: $(t)<200\n $(Cp)>=1000 and $(Cp)<=2000"')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runFilter',self.inputExperiment.get().getObjId(), self.filterType.get(), \
                                 self.condition.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFilter(self, objId, filterType, condition):
        import copy
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        self.printSection("Filtering")
        if self.filterType.get()==0:
            filterType="exclude"
        elif self.filterType.get()==1:
            filterType="keep"
        else:
            filterType="rmNA"

        filteredExperiment = PKPDExperiment()
        filteredExperiment.general = copy.copy(experiment.general)
        filteredExperiment.variables = copy.copy(experiment.variables)
        filteredExperiment.samples = {}
        filteredExperiment.doses = {}

        usedDoses = []
        for sampleKey, sample in experiment.samples.iteritems():
            candidateSample = PKPDSample()
            candidateSample.variableDictPtr    = copy.copy(sample.variableDictPtr)
            candidateSample.doseDictPtr        = copy.copy(sample.doseDictPtr)
            candidateSample.varName            = copy.copy(sample.varName)
            candidateSample.doseName           = copy.copy(sample.doseName)
            candidateSample.dose               = copy.copy(sample.dose)
            candidateSample.descriptors        = copy.copy(sample.descriptors)
            candidateSample.measurementPattern = copy.copy(sample.measurementPattern)

            N = 0 # Number of initial measurements
            if len(sample.measurementPattern)>0:
                aux=getattr(sample,"measurement_%s"%sample.measurementPattern[0])
                N = len(aux)
            if N==0:
                continue

            # Create empty output variables
            Nvar = len(sample.measurementPattern)
            for i in range(0,Nvar):
                exec("candidateSample.measurement_%s = []"%sample.measurementPattern[i])

            for n in range(0,N):
                toAdd = []
                okToAddTimePoint = True
                for i in range(0,Nvar):
                    exec("aux=sample.measurement_%s[%d]"%(sample.measurementPattern[i],n))
                    # aux=getattr(sample,"measurement_%s"%sample.measurementPattern[i])
                    if filterType=="rmNA":
                        if aux=="NA":
                            okToAddTimePoint = False
                        else:
                            toAdd.append(aux)
                if okToAddTimePoint:
                    for i in range(0,Nvar):
                        exec("candidateSample.measurement_%s.append('%s')"%(sample.measurementPattern[i],toAdd[i]))

            N = len(getattr(sample,"measurement_%s"%sample.measurementPattern[0])) # Number of final measurements
            if N!=0:
                filteredExperiment.samples[candidateSample.varName] = candidateSample
                usedDoses.append(candidateSample.doseName)

        if len(usedDoses)>0:
            for doseName in usedDoses:
                filteredExperiment.doses[doseName] = copy.copy(experiment.doses[doseName])

        self.writeExperiment(filteredExperiment,self._getPath("experiment.pkpd"))
        self.experiment = filteredExperiment

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------