# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Carlos Oscar Sorzano (info@kinestat.com)
# *              
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os
from os.path import basename, exists
import Tkinter as tk
import ttk

import pyworkflow.object as pwobj
from pyworkflow.wizard import Wizard
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from pyworkflow.gui.tree import BoundTree, TreeProvider

from protocol_pkpd_regression_labels import ProtPKPDRegressionLabel
from protocol_pkpd_drop_measurements import ProtPKPDDropMeasurements
from protocol_pkpd_change_units import ProtPKPDChangeUnits
from protocol_pkpd_scale_to_common_dose import ProtPKPDScaleToCommonDose
from protocol_pkpd_stats_oneExperiment_twoSubgroups_mean import ProtPKPDStatsExp1Subgroups2Mean
from protocol_pkpd_exponential_fit import ProtPKPDExponentialFit
from protocol_pkpd_elimination_rate import ProtPKPDEliminationRate
from protocol_pkpd_ev0_monocompartment import ProtPKPDEV0MonoCompartment
from protocol_pkpd_simulate_generic_pd import ProtPKPDSimulateGenericPD
from protocol_pkpd_stats_twoExperiments_twoSubgroups_mean import ProtPKPDStatsExp2Subgroups2Mean
from protocol_pkpd_import_from_csv import ProtPKPDImportFromText, getSampleNamesFromCSVfile
from protocol_pkpd_bootstrap_simulate import ProtPKPDODESimulate

class FilterVariablesTreeProvider(TreeProvider):
    """ Simplified view of VariablesTreeProvider with less columns.
    Additionally, we can filter by a given function. """

    def __init__(self, experiment, filter=None):
        self.experiment = experiment
        self.filter = filter or self._noFilter

    def _noFilter(self, v):
        return True

    def getColumns(self):
        return [('Name', 100), ('Unit', 100), ('Comment', 300)]

    def getObjects(self):
        sortedVars = []
        for key in sorted(self.experiment.variables.keys()):
            v = self.experiment.variables[key]
            if self.filter(v):
                sortedVars.append(v)
        return sortedVars

    def getObjectInfo(self, obj):
        key = obj.varName
        return {'key': key, 'text': key,
                'values': (obj.getUnitsString(),
                           obj.comment)}


class PKPDChooseVariableWizard(Wizard):
    _targets = [(ProtPKPDRegressionLabel, ['labelX']),
                (ProtPKPDRegressionLabel, ['labelY']),
                (ProtPKPDChangeUnits, ['labelToChange']),
                (ProtPKPDStatsExp1Subgroups2Mean, ['labelToCompare']),
                (ProtPKPDExponentialFit, ['predictor']),
                (ProtPKPDExponentialFit, ['predicted']),
                (ProtPKPDEliminationRate, ['predictor']),
                (ProtPKPDEliminationRate, ['predicted']),
                (ProtPKPDEV0MonoCompartment, ['predictor']),
                (ProtPKPDEV0MonoCompartment, ['predicted']),
                (ProtPKPDSimulateGenericPD, ['predictor']),
                (ProtPKPDStatsExp2Subgroups2Mean, ['label1', 'inputExperiment1']),
                (ProtPKPDStatsExp2Subgroups2Mean, ['label2', 'inputExperiment2'])
                ]

    def show(self, form, *params):
        protocol = form.protocol
        label = params[0]

        fullParams = self._getAllParams(protocol)
        experiment = None
        if len(fullParams)>1:
            experiment = protocol.getAttributeValue(fullParams[1], None)
        else:
            experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select input experiment first.")
        else:
            experiment.load()
            filterFunc = getattr(protocol, 'filterVarForWizard', None)
            tp = FilterVariablesTreeProvider(experiment, filter=filterFunc)
            dlg = dialog.ListDialog(form.root, self.getTitle(), tp,
                             selectmode=self.getSelectMode())
            if dlg.resultYes():
                self.setFormValues(form, label, dlg.values)

    def getTitle(self):
        return "Choose variable"

    def setFormValues(self, form, label, values):
        form.setVar(label, values[0].varName)

    def getSelectMode(self):
        return "browse"


class PKPDChooseSeveralVariableWizard(PKPDChooseVariableWizard):
    _targets = [(ProtPKPDDropMeasurements, ['varsToDrop']),
                (ProtPKPDScaleToCommonDose,['measurementsToChange'])]

    def getTitle(self):
        return "Choose variable(s)"

    def setFormValues(self, form, label, values):
        form.setVar(label, ', '.join([var.varName for var in values]))

    def getSelectMode(self):
        return "extended"


class PKPDVariableTemplateWizard(Wizard):
    _targets = [(ProtPKPDImportFromText, ['variables'])
                ]

    def show(self, form, *params):
        label = params[0]
        protocol = form.protocol
        currentValue = protocol.getAttributeValue(label, "")
        form.setVar(label, currentValue+"\n[Variable Name] ; [Units] ; [numeric/text] ; [time/label/measurement] ; [Comment]")

class PKPDDoseTemplateWizard(Wizard):
    _targets = [(ProtPKPDImportFromText, ['doses']),
                (ProtPKPDODESimulate, ['doses'])
                ]

    def show(self, form, *params):
        label = params[0]
        protocol = form.protocol
        currentValue = protocol.getAttributeValue(label, "")
        template = "\nInfusion0 ; infusion t=0.5...0.75 d=60*weight/1000; h; mg\n"\
                   "Bolus1 ; bolus t=2 d=100; h; mg\n"\
                   "Bolus0 ; bolus t=0 d=60*weight/1000; min; mg"
        form.setVar(label, currentValue+template)


class SimpleListTreeProvider(TreeProvider):
    """ A simple TreeProvider over the elements of a string list """

    def __init__(self, strList, name='Name', width=100):
        self.objects = [pwobj.Object(value=s) for s in strList]
        self.columns = [(name, width)]

    def getColumns(self):
        return self.columns

    def getObjects(self):
        return self.objects

    def getObjectInfo(self, obj):
        key = obj.get()
        return {'key': key, 'text': key,
                'values': ()}


class PKPDDosesToSamplesTemplateWizard(Wizard):
    _targets = [(ProtPKPDImportFromText, ['dosesToSamples'])
                ]
    def show(self, form, *params):
        label = params[0]
        protocol = form.protocol
        fnCSV = protocol.getAttributeValue('inputFile', "")
        if not os.path.exists(fnCSV):
            form.showError("Select a valid CSV input file first.")
        else:
            doseNames = []
            for line in protocol.doses.get().replace('\n',';;').split(';;'):
                tokens = line.split(';')
                if len(tokens)==4:
                    doseNames.append(tokens[0].strip())
            print(doseNames)

            sampleNames = getSampleNamesFromCSVfile(fnCSV)
            print(sampleNames)

            currentValue = protocol.getAttributeValue(label, "")
            tp = SimpleListTreeProvider(doseNames, name="Doses")
            dlg = dialog.ListDialog(form.root, "Test", tp,
                             selectmode='extended')
            if dlg.resultYes():
                self.setFormValues(form, label, dlg.values)
            #form.setVar(label, currentValue+"\n[Sample Name] ; [DoseName1,DoseName2,...]\n")
