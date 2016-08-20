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
In this module are protocol base classes related to PKPD

"""
import sys
import os
from pyworkflow.em.protocol import *
from pyworkflow.em.data import PKPDExperiment, PKPDFitting

class ProtPKPD(EMProtocol):
    def printSection(self, msg):
        print("**********************************************************************************************")
        print("Section: %s"%msg)
        print("**********************************************************************************************")

    def readExperiment(self,fnIn, show=True):
        experiment = PKPDExperiment()
        experiment.load(fnIn)
        if show:
            self.printSection("Reading %s"%fnIn)
            experiment._printToStream(sys.stdout)
        return experiment

    def writeExperiment(self, experiment, fnOut):
        self.printSection("Writing %s"%fnOut)
        experiment._printToStream(sys.stdout)
        experiment.write(fnOut)

    def readFitting(self, fnIn, show=True, cls=""):
        fitting = PKPDFitting(cls)
        fitting.load(fnIn)
        if show:
            self.printSection("Reading %s"%fnIn)
            fitting._printToStream(sys.stdout)
        return fitting

    def doublePrint(self,fh,msg):
        fh.write(msg+"\n")
        print(msg)

    def addFileContentToMessage(self,msg,fn):
        if os.path.exists(fn):
            fh = open(fn)
            for line in fh.readlines():
                msg.append(line.strip())
            fh.close()

    def loadInputExperiment(self):
        """ If the protocol has an attribute 'inputExperiment',
        load that experiment from file. If not, return None. """
        experiment = self.getAttributeValue('inputExperiment', None)

        if experiment:
            experiment.load()
            return experiment

        return None
