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

from os.path import basename
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
        return [('Name', 60), ('Unit', 60), ('Format', 60),
                ('Type', 60), ('Role', 60), ('Comment', 100)]

    def getObjects(self):
        return self.experiment.variables.values()

    def getObjectInfo(self, obj):
        key = obj.varName
        return {'key': key, 'text': key,
                'values': (obj.getUnitsString(),
                           obj.displayString,
                           obj.getTypeString(),
                           obj.getRoleString(),
                           obj.comment)}


class DosesTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Dose', 60), ('Amount', 60),
                ('Unit', 60), ('Norm.', 60)]

    def getObjects(self):
        return self.experiment.doses.values()

    def getObjectInfo(self, obj):
        key = obj.varName
        return {'key': key, 'text': key,
                'values': (obj.getDoseString(),
                           obj.doseAmount,
                           obj.getUnitsString(),
                           obj.normalization)
                }


class SamplesTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment
        sample = self.experiment.samples.values()[0]
        self.columns = [(key, 60) for key, _ in sample.descriptors.iteritems()]

    def getColumns(self):
        return [('Name', 60), ('Dose', 60)] + self.columns

    def getObjects(self):
        return self.experiment.samples.values()

    def getObjectInfo(self, obj):
        key = obj.varName
        values = [obj.doseName] + [obj.descriptors[k] for k, _ in self.columns]
        return {'key': key, 'text': key,
                'values': tuple(values)
                }


class ExperimentWindow(gui.Window):
    """ This class creates a Window that will display some Point's
    contained in a Data object.
    It will allow to launch 1D, 2D and 3D plots by selecting any
    combination of the x1, x2...xn from the Point dimension.
    Points can be selected by either Click and Drag in the Scatter plot or..
    by creating an Expression.
    Finally, there is a button 'Create Cluster' that will call a callback
    fuction to take care of it.
    """
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        self.experiment = kwargs.get('experiment')
        self.callback = kwargs.get('callback', None)
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)
        self.plotter = None

    def _createContent(self, content):
        # Create and fill the frame containing the Experiment
        # info, variables and doses
        self._createTopFrame(content)
        # Create the middle frame containing the Samples Box
        self._createSamplesFrame(content)
        # Create the last frame with the buttons
        self._createButtonsFrame(content)
        #self._updateSelectionLabel()

    def _createTopFrame(self, content):
        frame = tk.Frame(content)
        frame.columnconfigure(0, weight=1)
        tab = ttk.Notebook(frame)
        tab.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        #frame = tk.LabelFrame(content, text='General')
        def addLabelFrame(label, row, col, rowspan=1):
            lf = tk.Frame(tab)
            #frame.columnconfigure(col, minsize=250)
            #frame.columnconfigure(col, weight=1)#, minsize=30)
            tab.add(lf, text=label)
            #lf.grid(row=row, column=col, sticky='news', padx=5, pady=5,
            #        rowspan=rowspan)
            return lf

        lfGeneral = addLabelFrame('General', 0, 0, rowspan=2)
        self._addLabel(lfGeneral, 'Title', 0, 0)
        self._titleVar = tk.StringVar()
        self._titleVar.set(self.experiment.general['title'])
        titleEntry = tk.Entry(lfGeneral, width=26, textvariable=self._titleVar,
                              bg='white')
        titleEntry.grid(row=0, column=1, sticky='nw', padx=5, pady=(5, 0))
        self._addLabel(lfGeneral, 'Comment', 1, 0)
        commentText = gui.text.Text(lfGeneral, width=30, height=3, bg='white')
        commentText.setText(self.experiment.general['comment'])
        commentText.grid(row=1, column=1, sticky='nw', padx=5, pady=(5, 0))
        self._commentText = commentText

        lfVars = addLabelFrame('Variables', 0, 1)
        self.varsTree = self._addBoundTree(lfVars, VariablesTreeProvider, 5)
        lfDoses = addLabelFrame('Doses', 1, 1)
        self.dosesTree = self._addBoundTree(lfDoses, DosesTreeProvider, 5)

        frame.grid(row=0, column=0, sticky='new', padx=5, pady=(10, 5))

    def _createSamplesFrame(self, content):
        frame = tk.Frame(content)
        #frame = tk.LabelFrame(content, text='General')
        lfSamples = tk.LabelFrame(frame, text='Samples')
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.samplesTree = self._addBoundTree(lfSamples, SamplesTreeProvider, 10)
        self.samplesTree.itemDoubleClick = self._onSampleDoubleClick
        self.samplesTree.itemClick = self._onSampleClick

        plotFrame = tk.Frame(lfSamples)
        plotFrame.grid(row=1, column=0, sticky='ws', padx=5, pady=5)

        # Add a combobox with the variable for time
        timeVars = [v.varName for v in self.experiment.variables.values()
                    if v.role == v.ROLE_TIME]

        measureVars = [v.varName for v in self.experiment.variables.values()
                    if v.role == v.ROLE_MEASUREMENT]

        def addVar(text, col, choices):
            varFrame = tk.Frame(plotFrame)
            varFrame.grid(row=0, column=col, sticky='new')
            label = tk.Label(varFrame, text=text, font=self.fontBold)
            label.grid(row=0, column=0, padx=5, pady=2, sticky='nw')
            combo = ComboBox(varFrame, choices, width=10)
            combo.grid(row=0, column=1, sticky='nw', padx=5, pady=5)
            radioVar = tk.IntVar()
            radio = tk.Checkbutton(varFrame, text='Log', variable=radioVar)
            radio.grid(row=0, column=2, sticky='nw', padx=5, pady=5)
            return combo, radio, radioVar

        self.timeWidget = addVar('Time variable', 0, timeVars)
        self.measureWidget = addVar('Measure variable', 1, measureVars)

        self.plotButton = Button(plotFrame, '   Plot   ', font=self.fontBold,
                                 command=self._onPlotClick,
                                 tooltip='Select one or more samples to plot '
                                         'theirs measures of the selected '
                                         'variables (optionally in log).')

        self.plotButton.grid(row=0, column=2, sticky='ne', padx=5)

        frame.grid(row=1, column=0, sticky='news', padx=5, pady=5)

    def _createButtonsFrame(self, content):
        frame = tk.Frame(content)
        gui.configureWeigths(frame)
        buttonsFrame = tk.Frame(frame)
        buttonsFrame.grid(row=0, column=0, sticky='ne')
        closeButton = Button(buttonsFrame, 'Close', command=self.close,
                             imagePath='fa-times.png')
        closeButton.grid(row=0, column=0, sticky='ne', padx=5)
        self.newButton = HotButton(buttonsFrame, '   New Experiment   ',
                                   command=self._onCreateClick,
                                   tooltip='Create a new experiment with the '
                                            'selected samples. You can also edit'
                                            'title and comment.')

        self.newButton.grid(row=0, column=1, sticky='ne', padx=5)

        frame.grid(row=2, column=0, sticky='news', padx=5, pady=5)

    def _addLabel(self, parent, text, r, c):
        label = tk.Label(parent, text=text, font=self.fontBold)
        label.grid(row=r, column=c, padx=5, pady=5, sticky='ne')
        return label

    def _addBoundTree(self, parent, ProviderClass, height):
        bt = BoundTree(parent, ProviderClass(self.experiment), height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

    def _onPlotClick(self, e=None):
        sampleKeys = self.samplesTree.selection()
        if sampleKeys:
            samples = [self.experiment.samples[k] for k in sampleKeys]
            timeVarName = self.timeWidget[0].getText()
            measureVarName = self.measureWidget[0].getText()

            if self.plotter is None or self.plotter.isClosed():
                self.plotter = EmPlotter()
                doShow = True
                ax = self.plotter.createSubPlot("Plot", timeVarName, measureVarName)
                self.plotDict = {}
            else:
                doShow = False
                ax = self.plotter.getLastSubPlot()


            for s in samples:
                if not s.varName in self.plotDict:
                    x, y = s.getXYValues(timeVarName, measureVarName)
                    ax.plot(x, y, label=s.varName)
                    self.plotDict[s.varName] = True
            ax.legend()

            if doShow:
                self.plotter.show()
            else:
                self.plotter.draw()
        else:
            self.showInfo("Please select some sample(s) to plot.")

    def _onCreateClick(self, e=None):
        sampleKeys = self.samplesTree.selection()
        if sampleKeys:
            print "samples: ", len(sampleKeys)
            print "title: ", self._titleVar.get()
            print "comment: ", self._commentText.getText()
        else:
            self.showInfo("Please select some sample(s) to create a "
                          "new experiment.")

    def _onSampleClick(self, obj):
        if not (self.plotter is None or self.plotter.isClosed()):
            self._onPlotClick()

    def _onSampleDoubleClick(self, obj):
        MeasureWindow(masterWindow=self, experiment=self.experiment).show()

    def _onClosing(self):
        if self.plotter:
            self.plotter.close()
        gui.Window._onClosing(self)


class MeasureWindow(gui.Window):
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        self.experiment = kwargs.get('experiment')
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)

    def _createContent(self, content):
        self._createMeasurementFrame(content)
        self._createButtonsFrame(content)

    def _createMeasurementFrame(self, content):
        frame = tk.Frame(content)
        #frame = tk.LabelFrame(content, text='General')
        lfSamples = tk.LabelFrame(frame, text='Measurement')
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.samplesTree = self._addBoundTree(lfSamples, SamplesTreeProvider, 10)

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=5)

    def _createButtonsFrame(self, content):
        frame = tk.Frame(content)
        gui.configureWeigths(frame)
        buttonsFrame = tk.Frame(frame)
        buttonsFrame.grid(row=0, column=0, sticky='ne')
        closeButton = Button(buttonsFrame, 'Close', command=self.close,
                             imagePath='fa-times.png')

        frame.grid(row=2, column=0, sticky='news', padx=5, pady=5)

    def _addBoundTree(self, parent, ProviderClass, height):
        bt = BoundTree(parent, ProviderClass(self.experiment), height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

