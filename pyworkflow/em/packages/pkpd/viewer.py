# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
from os.path import basename, join, exists
import numpy as np

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.data import PKPDExperiment

from protocol_batch_create_experiment import BatchProtCreateExperiment
from protocol_pkpd_export_to_csv import ProtPKPDExportToCSV
from protocol_pkpd_statistics_labels import ProtPKPDStatisticsLabel
from pyworkflow.gui.text import openTextFileEditor

from tk_experiment import ExperimentWindow



class PKPDExperimentViewer(Viewer):
    """ Visualization of a given PKPDExperiment
    """
    _targets = [PKPDExperiment]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        obj.load()
        self.experimentWindow = self.tkWindow(ExperimentWindow,
                                           title='Experiment Viewer',
                                           experiment=obj,
                                           callback=self._createExperiment)
        self.experimentWindow.show()

    def _createExperiment(self):
        """ Create a new experiment after manipulation of
        the currently displayed experiment and register the action.
        """
        sampleKeys = self.experimentWindow.samplesTree.selection()
        samples = ';'.join([self.experimentWindow.experiment.samples[k].varName
                            for k in sampleKeys])

        prot = self.protocol
        project = prot.getProject()
        newProt = project.newProtocol(BatchProtCreateExperiment)
        newProt.inputExperiment.set(self.experimentWindow.experiment)
        newProt.listOfSamples.set(samples)
        newProt.newTitle.set(self.experimentWindow._titleVar.get())
        newProt.newComment.set(self.experimentWindow._commentText.getText())
        project.launchProtocol(newProt)


class PKPDCSVViewer(Viewer):
    """ Wrapper to visualize CSV files
    """
    _label = 'viewer csv'
    _targets = [ProtPKPDExportToCSV]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fnCSV = self.protocol.getFilenameOut()

        if exists(fnCSV):
            openTextFileEditor(fnCSV)


class PKPDStatisticsLabelViewer(Viewer):
    """ Wrapper to visualize statistics
    """
    _label = 'viewer statistics'
    _targets = [ProtPKPDStatisticsLabel]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fnStatistics = self.protocol._getPath("statistics.txt")

        if exists(fnStatistics):
            openTextFileEditor(fnStatistics)
