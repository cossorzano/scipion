# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

"""
Import experiment
"""
from base import ProtImportFiles
import pyworkflow.protocol.params as params
from pyworkflow.em.data import PKPDExperiment
from os.path import exists, basename
import sys
from pyworkflow.utils.path import copyFile

# TESTED in test_workflow_gabrielsson_pk01.py
# TESTED in test_workflow_gabrielsson_pk02.py
# TESTED in test_workflow_gabrielsson_pk03.py
# TESTED in test_workflow_gabrielsson_pk04.py
# TESTED in test_workflow_gabrielsson_pk05.py

class ProtImportExperiment(ProtImportFiles):
    """ Protocol to import an PKPD experiment\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'import experiment'
    IMPORT_FROM_FILES = 1 
    
    def __init__(self, **args):
        ProtImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputFile', params.PathParam,
                      label="File path", allowsNull=False,
                      help='Specify a path to desired experiment file.')
         
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep', self.inputFile.get())
        
    def createOutputStep(self, inputPath):
        localPath = self._getPath(basename(inputPath))
        experiment = PKPDExperiment()
        experiment.load(inputPath, verifyIntegrity=False)
        experiment._printToStream(sys.stdout)
        experiment.write(localPath)
        self._defineOutputs(outputExperiment=experiment)

    def _summary(self):
        summary = ['Experiment imported from file: *%s*' % self.inputFile]
        return summary
    
    def _validate(self):
        errors = []
        if (not exists(self.inputFile.get())):
            errors.append("Experiment not found at *%s*" % self.inputFile.get())
        return errors
