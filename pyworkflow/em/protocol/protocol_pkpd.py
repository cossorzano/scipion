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
import pyworkflow.protocol.params as params

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

    def setInputExperiment(self):
        """ Set as self.experiment the experiment that is referenced
        in the attribute self.inputExperiment.
        """
        self.setExperiment(self.loadInputExperiment())

def addDoseToForm(form):
    form.addParam('doses', params.TextParam, height=5, width=70, label="Doses", default="",
                  help="Structure: [Dose Name] ; [Via] ; [Description] ; [Units] ; [Optional]\n"\
                       "The dose name should have no space or special character\n"\
                       "Valid vias are: iv (intravenous), ev0 (extra-vascular order 0), ev1 (extra-vascular order 1), \n"\
                       "     ev01 (extra-vascular first order 0 and then order 1), evFractional (extra-vascular fractional order)\n"\
                       "Valid units are: h, mg, ug, ug/mL, ...\n"\
                       "Optional parameters are tlag (e.g. tlag=0)\n"\
                       "   and bioavailability (e.g. bioavailability=0.8)\n"\
                       "The description is either a bolus or an infusion as shown in the examples\n"\
                       "\nIt is important that there are two semicolons.\n"\
                       "Examples:\n"\
                       "Infusion0 ; infusion t=0.500000...0.750000 d=60*weight/1000; h; mg\n"\
                       "Bolus1 ; bolus t=2.000000 d=100; h; mg\n"\
                       "Treatment ; repeated_bolus t=0:8:48 d=100; h; mg")
