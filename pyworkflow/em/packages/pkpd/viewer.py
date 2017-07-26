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
from pyworkflow.gui.text import openTextFileEditor
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em.data import PKPDFitting, PKPDSignalAnalysis
from pyworkflow.em.pkpd_units import strUnit

from protocol_batch_create_experiment import BatchProtCreateExperiment
from protocol_pkpd_export_to_csv import ProtPKPDExportToCSV
from protocol_pkpd_statistics_labels import ProtPKPDStatisticsLabel
from protocol_pkpd_regression_labels import ProtPKPDRegressionLabel
from protocol_pkpd_ode_bootstrap import ProtPKPDODEBootstrap
from protocol_pkpd_filter_population import ProtPKPDFilterPopulation
from protocol_pkpd_merge_populations import ProtPKPDMergePopulations

from tk_experiment import ExperimentWindow
from tk_populations import PopulationWindow



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
        samples = ';'.join([self.experimentWindow.experiment.samples[k].sampleName
                            for k in sampleKeys])

        prot = self.protocol
        project = prot.getProject()
        newProt = project.newProtocol(BatchProtCreateExperiment)
        newProt.inputExperiment.set(self.experimentWindow.experiment)
        newProt.listOfSamples.set(samples)
        newProt.newTitle.set(self.experimentWindow._titleVar.get())
        newProt.newComment.set(self.experimentWindow._commentText.getText())
        project.launchProtocol(newProt)


class PKPDFittingViewer(Viewer):
    """ Visualization of a given PKPDFitting
    """
    _targets = [PKPDFitting]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fitting = obj
        if fitting.isPopulation():
            fitting.sampleFittingClass="PKPDSampleFitBootstrap"
            fitting.load()
            self.populationWindow = self.tkWindow(PopulationWindow,
                                                  title='Population Viewer',
                                                  population=fitting)
            self.populationWindow.show()
        else:
            fitting.load()
            experiment = fitting.loadExperiment()
            self.fittingWindow = self.tkWindow(ExperimentWindow,
                                               title='Fitting Viewer',
                                               experiment=experiment,
                                               fitting=fitting)
            self.fittingWindow.show()


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


class PKPDAnalysisViewer(Viewer):
    """ Wrapper to visualize analysis results
    """
    _label = 'viewer analysis'
    _targets = [PKPDSignalAnalysis]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fnAnalysis = obj.fnAnalysis.get()

        if exists(fnAnalysis):
            openTextFileEditor(fnAnalysis)


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


class PKPDRegressionLabelsViewer(Viewer):
    """ Wrapper to visualize regression
    """
    _label = 'viewer regression'
    _targets = [ProtPKPDRegressionLabel]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **kwargs):
        fnResults = self.protocol._getPath("results.txt")

        if exists(fnResults):
            X, Y, _, _ = self.protocol.getXYValues(False)

            minX = min(X)
            maxX = max(X)
            step = (maxX - minX) / 50
            xValues = np.arange(minX, maxX+step, step)
            yValues = self.protocol.evalFunction(xValues)

            plotter = EmPlotter()
            varNameX = self.protocol.labelX.get()
            varNameY = self.protocol.labelY.get()
            ax = plotter.createSubPlot("Regression Plot",
                                       "%s [%s]"%(varNameX,strUnit(self.protocol.experiment.variables[varNameX].units.unit)),
                                       "%s [%s]"%(varNameY,strUnit(self.protocol.experiment.variables[varNameY].units.unit)))
            ax.plot(xValues, yValues)
            ax.plot(X, Y, 'o')

            return [plotter]
        else:
            return [self.errorMessage("Result file '%s' not produced yet. "
                                      % fnResults)]


class PKPDPopulationViewer(Viewer):
    """ Visualization of a given Population
    """
    _targets = [ProtPKPDODEBootstrap,
                ProtPKPDFilterPopulation,
                ProtPKPDMergePopulations]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        if hasattr(obj,"outputPopulation"):
            population = PKPDFitting("PKPDSampleFitBootstrap")
            population.load(obj.outputPopulation.fnFitting)

            self.populationWindow = self.tkWindow(PopulationWindow,
                                                  title='Population Viewer',
                                                  population=population)
            self.populationWindow.show()

