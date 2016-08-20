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

import math
import numpy as np
from itertools import izip
import Tkinter as tk
import ttk

import pyworkflow.gui as gui
from pyworkflow.gui.widgets import Button, HotButton, ComboBox
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.em.plotter import EmPlotter


class VariablesTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Unit', 60)]

    def getObjects(self):
        sortedVars = []
        for key in sorted(self.experiment.variables.keys()):
            sortedVars.append(self.experiment.variables[key])
        return sortedVars

    def getObjectInfo(self, obj):
        key = obj.varName
        return {'key': key, 'text': key,
                'values': (obj.getUnitsString())}


class PopulationWindow(gui.Window):
    """
    """
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        self.experiment = kwargs.get('experiment')
        self.callback = kwargs.get('callback', None)
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)
        content.rowconfigure(0, weight=1)
        self.plotter = None

    def _createContent(self, content):
        # Create and fill the frame containing the Experiment
        # info, variables and doses
        self._createTopFrame(content)
        # Create the middle frame containing the Samples Box
        #self._createSamplesFrame(p)
        # Create the last frame with the buttons
        #self._createButtonsFrame(content)
        #self._updateSelectionLabel()

    def _createTopFrame(self, content):
        frame = tk.Frame(content)
        lfSamples = tk.LabelFrame(frame, text='Plot')
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.tree = self._addBoundTree(lfSamples, VariablesTreeProvider, 10)
        self.samplesTree.itemDoubleClick = self._onSampleDoubleClick
        self.samplesTree.itemClick = self._onSampleClick

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=(10, 5))

    def _addBoundTree(self, parent, ProviderClass, height):
        bt = BoundTree(parent, ProviderClass(self.experiment), height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

    def getUnits(self, varName):
        return self.experiment.variables[varName].getUnitsString()

    def getLabel(self, varName, useLog):
        varLabel = '%s [%s]' % (varName, self.getUnits(varName))
        if useLog:
            varLabel = "log10(%s)" % varLabel
        return varLabel

    def getTimeVarName(self):
        return self.timeWidget[0].getText()

    def useTimeLog(self):
        return self.timeWidget[2].get()

    def getTimeLabel(self):
        return self.getLabel(self.getTimeVarName(), self.useTimeLog())

    def getMeasureVarName(self):
        return self.measureWidget[0].getText()

    def useMeasureLog(self):
        return self.measureWidget[2].get()

    def getMeasureLabel(self):
        return self.getLabel(self.getMeasureVarName(), self.useMeasureLog())

    def getPlotValues(self, sample):
        xValues, yValues = sample.getXYValues(self.getTimeVarName(),
                                  self.getMeasureVarName())

        useMeasureLog = self.useMeasureLog()
        useTimeLog = self.useTimeLog()

        if not (useMeasureLog or useTimeLog):
            return xValues, yValues

        # If log will be used either for time or measure var
        # we need to filter elements larger than 0
        newXValues = []
        newYValues = []

        def _value(v, useLog):
            if useLog:
                return math.log10(v) if v > 0 else None
            return v

        for x, y in izip(xValues, yValues):
            x = _value(x, useTimeLog)
            y = _value(y, useMeasureLog)

            if x is not None and y is not None:
                newXValues.append(x)
                newYValues.append(y)

        return newXValues, newYValues

    def _onPlotClick(self, e=None):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            if self.plotter is None or self.plotter.isClosed():
                self.plotter = EmPlotter()
                doShow = True
                ax = self.plotter.createSubPlot("Plot", self.getTimeLabel(),
                                                self.getMeasureLabel())
                self.plotDict = {}
            else:
                doShow = False
                ax = self.plotter.getLastSubPlot()

            samples = [self.experiment.samples[k] for k in sampleKeys]
            for s in samples:
                if not s.varName in self.plotDict:
                    x, y = self.getPlotValues(s)
                    ax.plot(x, y, label=s.varName)
                    self.plotDict[s.varName] = True
            ax.legend()

            if doShow:
                self.plotter.show()
            else:
                self.plotter.draw()
        else:
            self.showInfo("Please select some sample(s) to plot.")

    def _onPlotSummaryClick(self, e=None):
        sampleKeys = self.samplesTree.selection()
        n = len(sampleKeys)

        if n == 1:
            self.showInfo("Please select several samples to plot.")
        else:
            if n > 1:
                samples = [self.experiment.samples[k] for k in sampleKeys]
            else:
                samples = self.experiment.samples.values()

            dataDict = {} # key will be time values
            for s in samples:
                xValues, yValues = self.getPlotValues(s)
                for x, y in izip(xValues, yValues):
                    if x in dataDict:
                        dataDict[x].append(y)
                    else:
                        dataDict[x] = [y]

            sortedTime = sorted(dataDict.keys())
            # We will store five values (min, 25%, 50%, 75%, max)
            # for each of the time entries computed
            percentileList = [0, 25, 50, 75, 100]
            Y = np.zeros((len(sortedTime), 5))
            for i, t in enumerate(sortedTime):
                Y[i,:] = np.percentile(dataDict[t], percentileList)

            plotter = EmPlotter()
            ax = plotter.createSubPlot("Summary Plot", self.getTimeLabel(),
                                       self.getMeasureLabel())
            ax.plot(sortedTime, Y[:, 0], 'r--', label="Minimum")
            ax.plot(sortedTime, Y[:, 1], 'b--', label="25%")
            ax.plot(sortedTime, Y[:, 2], 'g', label="50% (Median)")
            ax.plot(sortedTime, Y[:, 3], 'b--', label="75%")
            ax.plot(sortedTime, Y[:, 4], 'r--', label="Maximum")

            ax.legend()
            plotter.show()

    def _onCreateClick(self, e=None):
        sampleKeys = self.samplesTree.selection()
        if sampleKeys and self.callback:
            self.callback()
        else:
            self.showInfo("Please select some sample(s) to create a "
                          "new experiment.")

    def _onSampleClick(self, obj):
        if not (self.plotter is None or self.plotter.isClosed()):
            self._onPlotClick()

