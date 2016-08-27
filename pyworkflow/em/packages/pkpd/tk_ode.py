# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com)
# *              Carlos Oscar Sorzano (info@kinestat.com)
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

import math
import numpy as np
from itertools import izip
import Tkinter as tk
import ttk

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
import pyworkflow.gui.dialog as dialog
import pyworkflow.gui as gui
from pyworkflow.gui.widgets import Button, HotButton, ComboBox
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.em.plotter import EmPlotter


class SamplesTreeProvider(TreeProvider):
    def __init__(self, experiment, fitting=None):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Dose', 60)]

    def getObjects(self):
        sortedSamples = []
        for key in sorted(self.experiment.samples.keys()):
            sample = self.experiment.samples[key]
            sortedSamples.append(sample)
        return sortedSamples

    def getObjectInfo(self, obj):
        key = obj.varName
        values = [','.join(obj.doseList)]

        return {'key': key, 'text': key,
                'values': tuple(values)
                }


class MinMaxSlider(tk.Frame):
    """
    Create a personalized frame that contains:
        label, min entry, slider and max entry
    It also keeps a variable with the value
    """
    def __init__(self, master, label, from_=0, to=100, value=50, callback=None,
                 step=0.01):
        self.var = tk.DoubleVar()
        self.var.set(float(value))
        self.varMin = tk.DoubleVar()
        self.varMin.set(float(from_))
        self.varMax = tk.DoubleVar()
        self.varMax.set(float(to))

        tk.Frame.__init__(self, master)
        tk.Label(self, text=label).pack(side=tk.LEFT, padx=2, pady=2,
                                        anchor='s')

        def _entry(var):
            entry = tk.Entry(self, textvariable=var, bg='white', width=10)
            entry.pack(side=tk.LEFT, padx=2, pady=2, anchor='s')
            entry.bind('<Return>', self._onBoundChanged)

        _entry(self.varMin)

        self.slider = tk.Scale(self, from_=from_, to=to, variable=self.var,
                               bigincrement=step, resolution=step,
                               orient=tk.HORIZONTAL)
        self.slider.pack(side=tk.LEFT, padx=2)

        _entry(self.varMax)

        if callback:
            self.var.trace('w', callback)

    def getMinMax(self):
        return (self.varMin.get(), self.varMax.get())

    def getValue(self):
        return self.var.get()

    def _onBoundChanged(self, e=None):
        v = self.getValue()
        minV, maxV = self.getMinMax()
        self.slider.config(from_=minV, to=maxV)
        if v > maxV or v < minV:
            self.var.set(0.5 * (minV + maxV))


class PKPDODEDialog(dialog.Dialog):
    def __init__(self, parent, title, **kwargs):
        """ From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback: a callback function to validate selected items.
        """
        self.values = []
        self.protODE = kwargs['protODE']
        self.experiment = self.protODE.experiment
        self.varNameX = kwargs['varNameX']
        self.varNameY = kwargs['varNameY']
        self.provider = SamplesTreeProvider(self.experiment)
        self.model = self.loadModel()
        self.validateSelectionCallback = kwargs.get('validateSelectionCallback', None)

        dialog.Dialog.__init__(self, parent, title,
                        buttons=[('Select', dialog.RESULT_YES),
                                 ('Cancel', dialog.RESULT_CANCEL)])

    def loadModel(self):
        model = self.protODE.createModel()
        model.setExperiment(self.experiment)
        if hasattr(self.protODE, "deltaT"):
            model.deltaT = self.protODE.deltaT.get()
        model.setXVar(self.varNameX)
        model.setYVar(self.varNameY)
        return model

    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        gui.configureWeigths(bodyFrame)
        self._createSamplesFrame(bodyFrame)
        self._createSlidersFrame(bodyFrame)

    def _createSamplesFrame(self, content):
        frame = tk.Frame(content, bg='white')
        #frame = tk.LabelFrame(content, text='General')
        lfSamples = tk.LabelFrame(frame, text="Samples")
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.samplesTree = self._addBoundTree(lfSamples,
                                                  self.provider, 10)
        self.samplesTree.itemClick = self._onPlotFitClick

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=5)

    def _createSlidersFrame(self, content):
        frame = tk.Frame(content, bg='white')
        lfBounds = tk.LabelFrame(frame, text="Parameter Bounds")
        gui.configureWeigths(frame)

        i = 0
        self.sliders = {}

        for paramName, bounds in self.protODE.getParameterBounds().iteritems():
            bounds = bounds or (0, 1)
            value = 0.5 * (bounds[0] + bounds[1])
            slider = MinMaxSlider(lfBounds, paramName,
                                  bounds[0], bounds[1], value,
                                  callback=self._onVarChanged)
            slider.grid(row=i, column=0, padx=5, pady=5)
            self.sliders[paramName] = slider
            i += 1

        lfBounds.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        frame.grid(row=0, column=1, sticky='news', padx=5, pady=5)

    def _addBoundTree(self, parent, provider, height):
        bt = BoundTree(parent, provider, height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

    def apply(self):
        self.values = []

    def _onVarChanged(self, *args):
        for paramName, slider in self.sliders.iteritems():
            print "paramName: ", paramName, "values: ", slider.getMinMax(), slider.getValue()

    def getPlotValues(self, sample):
        xValues, yValues = sample.getXYValues(self.varNameX, self.varNameY)

        useMeasureLog = False
        useTimeLog = False

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

    def _onPlotFitClick(self, e=None):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            # Get first selected element
            #fit = self.fitting.getSampleFit(sampleKeys[0])
            plotter = EmPlotter()
            ax = plotter.createSubPlot("Plot", self.varNameX, self.varNameY)

            sample = self.experiment.samples[sampleKeys[0]]
            x, y = self.getPlotValues(sample)
            ax.plot(x, y, 'x', label="Observations")
            ax.legend()

            plotter.show()
        else:
            self.showInfo("Please select some sample(s) to plot.")
