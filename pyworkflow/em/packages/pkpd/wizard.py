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

from pyworkflow.wizard import Wizard
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from pyworkflow.gui.tree import BoundTree, TreeProvider

from protocol_pkpd_regression_labels import ProtPKPDRegressionLabel


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
                (ProtPKPDRegressionLabel, ['labelY'])]

    def show(self, form, *params):
        protocol = form.protocol
        label = params[0]

        value = protocol.getAttributeValue(label)
        experiment = protocol.loadInputExperiment()

        if experiment is None:
            dialog.showError("Select input experiment first.")
        else:
            filterFunc = getattr(protocol, 'filterVarForWizard', None)
            tp = FilterVariablesTreeProvider(experiment, filter=filterFunc)
            dlg = dialog.ListDialog(form.root, "Choose variable", tp,
                             selectmode='browse')
            if dlg.resultYes():
                form.setVar(label, dlg.values[0].varName)


